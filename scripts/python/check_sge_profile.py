#!/usr/bin/env python3
"""
check_sge_profile.py
--------------------
Validate an SGE workflow profile before any job is submitted. Three things are
checked, each of which otherwise fails silently or late.

  1. Wall clock. A rule whose `runtime` exceeds the `h_rt` ceiling of the queue
     it is routed to is either rejected at submission or, worse, admitted and
     killed at the ceiling hours later. Snakemake's `restart-times` then retries
     it with the same doomed request.

  2. Memory complex. `-l h_vmem=...` only reserves memory if h_vmem is
     CONSUMABLE. A complex that is requestable but not consumable (mem_free is
     usually configured this way) is checked once at dispatch and never tracked,
     so jobs oversubscribe the node and the OOM killer arbitrates.

  3. Existence. Queues and the parallel environment must exist, or every qsub
     fails identically and the run dies at job one.

WHERE THE FACTS COME FROM, and why this is not as simple as calling qconf.

`qconf -sc`, `qconf -spl` and `qconf -sq` are permitted only from an ADMIN HOST.
From an ordinary submit host they fail with

    denied: host "morty.eb.local" is not an admin host

That is a refusal to answer, not the answer "no". An earlier version of this
script reported `queue 'standard.q' not found` on exactly that stderr, which is
the same undefined-treated-as-negative mistake the pipeline audit was written to
eliminate. So this script now distinguishes three outcomes per check --
VERIFIED, ASSUMED (from the recorded facts in cluster.yaml) and UNVERIFIABLE --
and only ever reports a problem from a query that actually answered.

Three sources, tried in order:

  qconf         authoritative, admin hosts only. When it answers, its values are
                cross-checked against cluster.yaml and drift is reported.
  cluster.yaml  facts recorded earlier on an admin host, with provenance. Used
                when qconf is denied. Regenerate with record_cluster_facts.sh.
  qsub -w v     a dry scheduling run on an empty cluster (see qsub(1) -w v). It
                does not submit, works from any submit host, and rejects a
                request whose h_rt exceeds the queue ceiling, whose complex is
                unknown, or whose PE or queue does not exist. It is the scheduler
                itself answering, so it settles points 1 and 3 without qconf --
                but it cannot see whether a complex is consumable, which is why
                point 2 still needs qconf or a recorded fact.

Usage:
    python scripts/python/check_sge_profile.py config/sge
    python scripts/python/check_sge_profile.py config/sge --no-qsub-verify
    python scripts/python/check_sge_profile.py config/sge --facts other/cluster.yaml

Exit 0 = nothing contradicted. Exit 1 = a check that answered, answered badly.
Unverifiable checks are printed and do not fail the run; they are listed at the
end with the command to run on an admin host.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

import yaml

INF = float("inf")
DENIED = re.compile(r"denied:|not an admin host|must be (a )?(manager|operator)",
                    re.I)


class Probe:
    """Result of running a scheduler query. `denied` is not `absent`.

    `text` holds stdout+stderr whatever the exit status, because a command can
    fail *and* have answered: `qsub -w v` exits 1 to say "no suitable queues",
    which is a real answer, not a refusal to answer.
    """

    def __init__(self, out=None, denied=False, missing=False, err="", text="", rc=0):
        self.out, self.denied, self.missing = out, denied, missing
        self.err, self.text, self.rc = err, text, rc

    @property
    def answered(self):
        return self.out is not None


def probe(cmd):
    if not shutil.which(cmd[0]):
        return Probe(missing=True, err=f"{cmd[0]} not on PATH")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as e:
        return Probe(err=str(e))
    text = (r.stderr or "") + (r.stdout or "")
    if DENIED.search(text):
        return Probe(denied=True, text=text, rc=r.returncode,
                     err=text.strip().splitlines()[0] if text.strip() else "denied")
    if r.returncode != 0:
        return Probe(err=text.strip(), text=text, rc=r.returncode)
    return Probe(out=r.stdout, text=text, rc=0)


def parse_hms(v):
    """'24:00:00' -> seconds. 'INFINITY' -> inf."""
    v = str(v).strip()
    if v.upper() in ("INFINITY", "INFINITE", "NONE", ""):
        return INF
    parts = [int(x) for x in v.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    h, m, s = parts[-3:]
    return h * 3600 + m * 60 + s


def fmt_h(sec):
    return "unlimited" if sec is INF else f"{sec / 3600:.0f}h"


# --------------------------------------------------------------------------
class Facts:
    """Cluster facts, each tagged with how we came to believe it."""

    def __init__(self, recorded):
        self.recorded = recorded or {}
        self.pes = None
        self.pes_src = None
        self.complexes = {}
        self.complex_src = {}
        self.h_rt = {}
        self.h_rt_src = {}
        self.drift = []
        self.unverifiable = []

    def rec_queue(self, q, key):
        return ((self.recorded.get("queues") or {}).get(q) or {}).get(key)

    def load(self, queues, mem_complex):
        # --- parallel environments -------------------------------------------
        p = probe(["qconf", "-spl"])
        if p.answered:
            self.pes, self.pes_src = set(p.out.split()), "qconf"
            rec = set(self.recorded.get("parallel_environments") or [])
            if rec and rec != self.pes:
                self.drift.append(f"parallel environments: cluster.yaml says "
                                  f"{sorted(rec)}, qconf says {sorted(self.pes)}")
        elif self.recorded.get("parallel_environments"):
            self.pes, self.pes_src = set(self.recorded["parallel_environments"]), "cluster.yaml"
        else:
            self.unverifiable.append(("parallel environment", "qconf -spl", p.err))

        # --- memory complex ---------------------------------------------------
        p = probe(["qconf", "-sc"])
        if p.answered:
            for line in p.out.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                if len(f) >= 6:
                    self.complexes[f[0]] = (f[4].upper(), f[5].upper())
                    self.complex_src[f[0]] = "qconf"
            rc = (self.recorded.get("complexes") or {}).get(mem_complex)
            got = self.complexes.get(mem_complex)
            if rc and got:
                want = {"per_slot": "YES", "per_job": "JOB", False: "NO"}.get(rc.get("consumable"))
                if want and got[1] != want:
                    self.drift.append(f"{mem_complex}: cluster.yaml says consumable="
                                      f"{rc.get('consumable')}, qconf says {got[1]}")
        else:
            rc = (self.recorded.get("complexes") or {}).get(mem_complex)
            if rc:
                cons = {"per_slot": "YES", "per_job": "JOB", False: "NO", None: "NO"}[rc.get("consumable")]
                self.complexes[mem_complex] = ("YES" if rc.get("requestable") else "NO", cons)
                self.complex_src[mem_complex] = "cluster.yaml"
            else:
                self.unverifiable.append(("memory complex consumability",
                                          f"qconf -sc | grep {mem_complex}", p.err))

        # --- queue h_rt ceilings ----------------------------------------------
        for q in queues:
            p = probe(["qconf", "-sq", q])
            if p.answered:
                m = re.search(r"^h_rt\s+(\S+)", p.out, re.M)
                self.h_rt[q] = parse_hms(m.group(1)) if m else INF
                self.h_rt_src[q] = "qconf"
                rec = self.rec_queue(q, "h_rt")
                if rec and parse_hms(rec) != self.h_rt[q]:
                    self.drift.append(f"{q} h_rt: cluster.yaml says {rec}, qconf says "
                                      f"{fmt_h(self.h_rt[q])}")
            elif p.denied or p.missing:
                rec = self.rec_queue(q, "h_rt")
                if rec:
                    self.h_rt[q] = parse_hms(rec)
                    self.h_rt_src[q] = "cluster.yaml"
            else:
                # qconf answered and said no such queue: that IS an absence
                self.h_rt[q] = None
                self.h_rt_src[q] = "qconf"


# --------------------------------------------------------------------------
def qsub_verify(queue, pe, threads, mem_complex, mem_mb_per_slot, runtime_min):
    """`qsub -w v`: dry scheduling run on an empty cluster. Does not submit.

    Returns (ok, detail). ok is None when qsub could not be asked.
    """
    h = runtime_min // 60
    m = runtime_min % 60
    args = ["qsub", "-w", "v", "-b", "y", "-q", queue,
            "-l", f"h_rt={h:02d}:{m:02d}:00"]
    if mem_complex != "none":
        args += ["-l", f"{mem_complex}={mem_mb_per_slot}M"]
    if threads > 1:
        args += ["-pe", pe, str(threads)]
    args += ["true"]
    p = probe(args)
    if p.denied or p.missing:
        return None, p.err
    # qsub -w v exits 1 to say "no suitable queues". That is an ANSWER, so read
    # the text before the exit status; only a refusal or a crash is unverifiable.
    text = (p.text or "").lower()
    if not text:
        return None, p.err or "no output"
    if "no suitable queue" in text or "error:" in text or "unknown resource" in text:
        return False, p.text.strip().splitlines()[-1]
    if "suitable queue" in text or "verification" in text:
        return True, ""
    return None, p.text.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("profile", help="directory containing config.yaml")
    ap.add_argument("--facts", default=None,
                    help="recorded cluster facts (default: <profile>/cluster.yaml, "
                         "then config/sge/cluster.yaml)")
    ap.add_argument("--mem-complex", default="h_vmem")
    ap.add_argument("--pe", default=None)
    ap.add_argument("--margin", type=float, default=0.10,
                    help="required headroom below a queue's h_rt (default 10%%)")
    ap.add_argument("--no-qsub-verify", action="store_true",
                    help="skip the `qsub -w v` dry scheduling run")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(f"{args.profile}/config.yaml")) or {}
    dflt = cfg.get("default-resources") or {}
    setres = cfg.get("set-resources") or {}
    default_queue = dflt.get("sge_queue")
    default_runtime = dflt.get("runtime", 0)
    default_mem = dflt.get("mem_mb", 0)
    pe = args.pe or cfg.get("sge-pe") or "parallel"

    factfile = args.facts
    for cand in (factfile, f"{args.profile}/cluster.yaml", "config/sge/cluster.yaml"):
        if cand and os.path.exists(cand):
            factfile = cand
            break
    else:
        factfile = None
    recorded = yaml.safe_load(open(factfile)) if factfile else {}

    # Threads decide two things at once: `-pe parallel N`, and the divisor that
    # turns a job's total mem_mb into the per-slot request. Taking them from
    # set-threads means the verified request is the request that gets submitted.
    # A rule with no entry falls back to 1, which is what Snakemake gives a rule
    # that declares no `threads:` -- and is flagged, because the alternative is
    # to guess a number and verify a request nobody will ever make.
    setthreads = cfg.get("set-threads") or {}
    assumed = [r for r in setres if r not in setthreads]

    rules = [(r, (v or {}).get("sge_queue", default_queue),
              (v or {}).get("runtime", default_runtime),
              (v or {}).get("mem_mb", default_mem),
              int(setthreads.get(r, 1)))
             for r, v in setres.items()]
    rules.append(("<default>", default_queue, default_runtime, default_mem, 1))
    queues = sorted({q for _, q, _, _, _ in rules if q})

    F = Facts(recorded)
    F.load(queues, args.mem_complex)

    problems, notes = [], []
    if factfile and any(s == "cluster.yaml" for s in
                        list(F.h_rt_src.values()) + [F.pes_src] + list(F.complex_src.values())):
        notes.append(f"qconf is not answering from this host, so some facts come from "
                     f"{factfile} (recorded {recorded.get('recorded_on', '?')} on "
                     f"{recorded.get('recorded_from', '?')}). Re-record on an admin host "
                     f"if the cluster has changed.")

    # --- parallel environment -------------------------------------------------
    if F.pes is not None and pe not in F.pes:
        problems.append(f"parallel environment {pe!r} does not exist "
                        f"(known: {' '.join(sorted(F.pes))}) [{F.pes_src}]")

    # --- memory complex -------------------------------------------------------
    flags = F.complexes.get(args.mem_complex)
    per_job = False
    if flags is None and args.mem_complex not in ("none",):
        if F.complexes and F.complex_src.get(args.mem_complex) is None and \
                any(s == "qconf" for s in F.complex_src.values()):
            problems.append(f"complex {args.mem_complex!r} is not defined [qconf]")
    elif flags:
        requestable, consumable = flags
        src = F.complex_src.get(args.mem_complex, "?")
        if requestable == "NO":
            problems.append(f"complex {args.mem_complex!r} is not requestable [{src}]")
        if consumable == "NO":
            problems.append(
                f"complex {args.mem_complex!r} is requestable but NOT consumable [{src}]: "
                f"requesting it reserves nothing, so jobs will oversubscribe the node. "
                f"Use a consumable complex.")
        elif consumable == "JOB":
            per_job = True
            notes.append(f"{args.mem_complex} is consumable=JOB (per job, not per slot) "
                         f"[{src}]: do NOT divide mem_mb by threads. Set SGE_MEM_PER_JOB=1 "
                         f"for config/sge-generic.")
        else:
            notes.append(f"{args.mem_complex} is consumable per slot [{src}]: mem_mb is "
                         f"divided by threads before it is requested.")

    # --- queues exist ---------------------------------------------------------
    p = probe(["qstat", "-g", "c"])
    if p.answered:
        known = set(re.findall(r"^(\S+\.q)", p.out, re.M))
        for q in queues:
            if known and q not in known:
                problems.append(f"queue {q!r} is not in `qstat -g c`")
    else:
        F.unverifiable.append(("queue existence", "qstat -g c", p.err))

    for q in queues:
        if F.h_rt.get(q) is None and F.h_rt_src.get(q) == "qconf":
            problems.append(f"queue {q!r} does not exist [qconf]")

    if assumed:
        notes.append(f"no set-threads entry for {', '.join(sorted(assumed))}; assuming 1 "
                     f"slot each (which is what Snakemake gives a rule with no `threads:` "
                     f"directive). State them in the profile so the verified request is "
                     f"the request that gets submitted.")

    # --- runtime vs ceiling ---------------------------------------------------
    checked = 0
    for rule, q, runtime, _, _ in sorted(rules):
        ceiling = F.h_rt.get(q)
        if not q or ceiling is None or ceiling is INF:
            continue
        checked += 1
        src = F.h_rt_src.get(q, "?")
        want = float(runtime) * 60
        if want > ceiling:
            problems.append(f"{rule}: runtime {runtime} min exceeds {q} h_rt "
                            f"{fmt_h(ceiling)} [{src}] -- the job will be rejected or killed")
        elif want > ceiling * (1 - args.margin):
            problems.append(f"{rule}: runtime {runtime} min is within {args.margin:.0%} of "
                            f"{q}'s h_rt ceiling ({fmt_h(ceiling)}) [{src}]. Leave headroom "
                            f"or move it to a longer queue")
    if not checked:
        F.unverifiable.append(("runtime vs queue h_rt ceiling",
                               "qconf -sq <queue> | grep h_rt", "no ceilings known"))

    # --- ask the scheduler itself (works from any submit host) ----------------
    verified = 0
    if not args.no_qsub_verify and shutil.which("qsub"):
        seen = set()
        for rule, q, runtime, mem_mb, threads in sorted(rules):
            if not q:
                continue
            per_slot = mem_mb if per_job else max(1, -(-int(mem_mb) // threads))
            key = (q, runtime, per_slot, threads)
            if key in seen:
                continue
            seen.add(key)
            ok, detail = qsub_verify(q, pe, threads, args.mem_complex, per_slot, int(runtime))
            if ok is False:
                problems.append(f"{rule}: `qsub -w v` rejects the request "
                                f"(-q {q} -pe {pe} {threads} -l h_rt={runtime}min,"
                                f"{args.mem_complex}={per_slot}M): {detail}")
            elif ok is True:
                verified += 1
            elif detail:
                F.unverifiable.append((f"qsub -w v for {rule}", "qsub -w v ...", detail))
    elif not args.no_qsub_verify:
        F.unverifiable.append(("scheduler dry run", "qsub -w v ...", "qsub not on PATH"))

    # --- report ---------------------------------------------------------------
    for n in notes:
        print(f"note:  {n}")
    for dft in F.drift:
        print(f"DRIFT: {dft}  <- cluster.yaml is stale; re-record it")
    for what, cmd, why in F.unverifiable:
        print(f"unverifiable from this host: {what}\n"
              f"       run on an admin host: {cmd}\n"
              f"       ({why})")

    if problems:
        print(f"\n{len(problems)} problem(s) in {args.profile}:", file=sys.stderr)
        for p_ in problems:
            print(f"  - {p_}", file=sys.stderr)
        return 1

    src_summary = ", ".join(f"{q} h_rt {fmt_h(c)} [{F.h_rt_src.get(q)}]"
                            for q, c in F.h_rt.items() if c is not None)
    print(f"\n{args.profile}: OK. {len(rules)} rule routings checked."
          + (f" Queues: {src_summary}." if src_summary else "")
          + (f" {verified} distinct request(s) accepted by `qsub -w v`." if verified else ""))
    if F.drift:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
