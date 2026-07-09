#!/usr/bin/env bash
# Verification suite. Run from the repository root:
#
#     bash tests/run_tests.sh
#
# 1. compile   every Python script parses
# 2. config    config.yaml parses and every key the scripts read exists
# 3. rules     every script referenced by a Snakemake rule exists on disk
# 4. stats     hand-derived numerical checks of each statistic
# 5. smoke     end-to-end runs on synthetic data with planted signals
#
# Needs: python3, pandas, numpy, scipy, statsmodels, pyyaml. pyarrow optional
# (tests/_shim.py backs parquet with pickle when it is absent). R, Snakemake and
# the external tools are NOT needed: those components are covered by the static
# checks in steps 1-3 and by review, as recorded in docs/AUDIT.md.
set -uo pipefail
cd "$(dirname "$0")/.."
FAIL=0
hdr() { printf '\n\033[1m== %s ==\033[0m\n' "$1"; }

hdr "1. compile all Python scripts"
python3 -m compileall -q scripts/python tests && echo "OK" || FAIL=1

hdr "2. config parses; no CHANGE_ME left in required inputs"
python3 - <<'PY' || FAIL=1
import sys, yaml
cfg = yaml.safe_load(open("config/config.yaml"))
need = [("inputs","metadata"), ("inputs","focal_niche"), ("stats","contrasts"),
        ("stats","consensus"), ("stats","fdr_alpha")]
missing = [".".join(k) for k in need
           if cfg.get(k[0], {}).get(k[1]) is None]
if missing:
    sys.exit("missing config keys: " + ", ".join(missing))
c = cfg["stats"]["consensus"]
for k in ("require_methods","require_same_direction","min_applicable_methods",
          "anchor_method","strict_all_methods"):
    if k not in c: sys.exit(f"stats.consensus.{k} missing")
if c["anchor_method"] not in c["require_methods"]:
    sys.exit("stats.consensus.anchor_method must be in require_methods")
print("OK")
PY

hdr "3. every script referenced by a rule exists; no {genome} templating in rules"
python3 - <<'PY' || FAIL=1
import glob, os, re, sys
missing, templated = [], []
for f in glob.glob("workflow/**/*.smk", recursive=True) + ["workflow/Snakefile"]:
    text = open(f).read()
    for m in re.finditer(r'scripts/(python|R|sh)/([\w.-]+\.(?:py|R|sh))', text):
        p = os.path.join("scripts", m.group(1), m.group(2))
        if not os.path.exists(p):
            missing.append(f"{f}: {p}")
    # annotation paths are resolved once by resolve_annotation_paths.py; a rule
    # expanding {genome} itself would bypass the manifest
    if re.search(r'sed\s+"s/\{genome\}', text):
        templated.append(f)
if missing:
    sys.exit("missing scripts:\n  " + "\n  ".join(sorted(set(missing))))
if templated:
    sys.exit("rules expanding {genome} directly (use annot_lookup): " + ", ".join(templated))
print("OK")
PY

hdr "3c. annot_lookup projects the requested columns for the requested genomes"
python3 - <<'PY' || FAIL=1
import os, subprocess, sys, tempfile
src = open("workflow/Snakefile").read()
i = src.index("def annot_lookup(")
j = src.index("\n\n", src.index("{ids} {paths or ANNOT_PATHS}"))
ns = {"ANNOT_PATHS": "ap.tsv"}
exec(src[i:j], ns)
al = ns["annot_lookup"]
d = tempfile.mkdtemp(); here = os.getcwd(); os.chdir(d)
open("ap.tsv", "w").write("genome\tprokka_faa\tdbcan_overview\n"
                          "A\t/x/a.faa\t/o/A.txt\nB\t/x/b.faa\t/o/B.txt\nC\t/x/c.faa\t/o/C.txt\n")
open("ids.txt", "w").write("A\nC\n")
def sh(c): return subprocess.run(c, shell=True, capture_output=True, text=True)
r = sh(al("ids.txt", ["prokka_faa", "dbcan_overview"], "ap.tsv"))
assert r.returncode == 0, r.stderr
assert r.stdout == "A\t/x/a.faa\t/o/A.txt\nC\t/x/c.faa\t/o/C.txt\n", repr(r.stdout)
assert sh(al("ids.txt", ["nope"], "ap.tsv")).returncode == 2, "missing column must fail"
open("ap2.tsv", "w").write("genome\tprokka_faa\nA\t\nC\t/x/c.faa\n")
assert sh(al("ids.txt", "prokka_faa", "ap2.tsv")).returncode == 2, "empty cell must fail"
os.chdir(here)
print("OK")
PY

hdr "3a. no conda environment can resolve from defaults/anaconda"
python3 scripts/python/check_env_channels.py envs || FAIL=1

hdr "3b. every R package attached with library() is declared in envs/r.yaml"
python3 - <<'PY' || FAIL=1
import glob, re, sys, yaml
declared = set()
for f in glob.glob("envs/*.yaml"):
    for dep in (yaml.safe_load(open(f)).get("dependencies") or []):
        if isinstance(dep, str):
            name = re.split(r"[=<>]", dep.strip())[0]
            declared.add(name)
            if name.startswith("r-"):
                declared.add(name[2:])
            if name.startswith("bioconductor-"):
                declared.add(name[len("bioconductor-"):])
# packages that ship with R itself
BASE = {"stats","utils","methods","grDevices","graphics","tools","parallel","grid","base"}
missing = []
for f in sorted(glob.glob("scripts/R/*.R")):
    for m in re.finditer(r'(?:library|require)\(\s*["\']?([\w.]+)', open(f).read()):
        p = m.group(1)
        if p not in BASE and p.lower() not in {d.lower() for d in declared}:
            missing.append(f"{f}: {p}")
if missing:
    sys.exit("R packages attached but not in any envs/*.yaml:\n  " + "\n  ".join(missing))
print("OK")
PY

hdr "3d. cluster profiles: valid YAML, and every per-rule key names a real rule"
python3 - <<'PY' || FAIL=1
import glob, re, sys, yaml
rules = set()
for f in glob.glob("workflow/**/*.smk", recursive=True) + ["workflow/Snakefile"]:
    rules |= set(re.findall(r'^(?:rule|checkpoint)\s+([\w]+):', open(f).read(), re.M))
if not rules:
    sys.exit("no rules found; the rule-name regex is wrong")
bad = []
for prof in sorted(glob.glob("config/*/config.yaml")):
    try:
        cfg = yaml.safe_load(open(prof)) or {}
    except yaml.YAMLError as e:
        sys.exit(f"{prof} is not valid YAML: {e}")
    # Snakemake silently ignores a set-resources/set-threads key that is not a
    # rule name, so the job runs with the defaults and you find out from an OOM.
    for section in ("set-resources", "set-threads"):
        for name in (cfg.get(section) or {}):
            if name not in rules:
                bad.append(f"{prof}: {section}: {name!r} is not a rule")
    if cfg.get("executor") == "cluster-generic":
        for key in ("cluster-generic-submit-cmd", "cluster-generic-status-cmd"):
            if key not in cfg:
                bad.append(f"{prof}: missing {key}")
        for s in ("qsub_submit.sh", "qsub_status.sh"):
            import os
            if not os.path.exists(f"{os.path.dirname(prof)}/{s}"):
                bad.append(f"{prof}: missing {s}")

    # Threads drive `-pe parallel N` and the per-slot memory divisor, and on SLURM
    # they must match cpus_per_task or the job idles cores it reserved. Every rule
    # given resources must state its threads; a rule that declares no `threads:`
    # silently gets 1, which is fine but must be said out loud.
    st = cfg.get("set-threads") or {}
    for name, res in (cfg.get("set-resources") or {}).items():
        if name not in st:
            bad.append(f"{prof}: set-resources has {name!r} but set-threads does not")
        elif (res or {}).get("cpus_per_task") not in (None, st[name]):
            bad.append(f"{prof}: {name}: cpus_per_task={res['cpus_per_task']} but "
                       f"set-threads={st[name]}")

# every profile must agree on how many threads a rule gets
by_rule = {}
for prof in sorted(glob.glob("config/*/config.yaml")):
    for name, t in ((yaml.safe_load(open(prof)) or {}).get("set-threads") or {}).items():
        by_rule.setdefault(name, {})[prof] = t
for name, seen in by_rule.items():
    if len(set(seen.values())) > 1:
        bad.append(f"profiles disagree on threads for {name!r}: {seen}")

if bad:
    sys.exit("cluster profile problems:\n  " + "\n  ".join(bad))
print(f"OK ({len(glob.glob('config/*/config.yaml'))} profiles checked against {len(rules)} rules)")
PY

hdr "3e. SGE submit/status scripts convert resources correctly"
python3 - <<'PY' || FAIL=1
import os, subprocess, sys, tempfile, textwrap
d = tempfile.mkdtemp(); bin = f"{d}/bin"; os.makedirs(bin)
open(f"{d}/js.sh", "w").write("#!/bin/sh\n")
def fake(name, body):
    p = f"{bin}/{name}"
    open(p, "w").write("#!/bin/sh\n" + textwrap.dedent(body))
    os.chmod(p, 0o755)
fake("qsub", '''echo "$*" >&2\necho "  12345"\n''')
env = {**os.environ, "PATH": f"{bin}:{os.environ['PATH']}", "TMPDIR": d}
sub = os.path.abspath("config/sge-generic/qsub_submit.sh")
sta = os.path.abspath("config/sge-generic/qsub_status.sh")
chk = os.path.abspath("scripts/python/check_sge_profile.py")
profiles = [os.path.abspath(p) for p in ("config/sge", "config/sge-generic")]

def submit(*a, **kw):
    r = subprocess.run(["bash", sub, *map(str, a), f"{d}/js.sh"], capture_output=True,
                       text=True, env={**env, **kw}, cwd=d)
    assert r.returncode == 0, r.stderr
    return r.stdout, r.stderr

out, err = submit("phyloglm_chunk", 8, 32000, 360, "long.q")
assert out == "12345", f"submit must print a bare job id, got {out!r}"
assert "-l h_vmem=4000M" in err, err          # 32000 MB total / 8 slots
assert "-l h_rt=06:00:00" in err, err         # 360 min
assert "-pe parallel 8" in err, err
assert "-q long.q" in err, err

_, err = submit("ingest", 1, 8000, 120, "standard.q")
assert "-pe" not in err, "single-threaded jobs must not request a PE"
assert "-l h_vmem=8000M" in err, err

_, err = submit("scaffold", 8, 320000, 1440, "long.q")
assert "-l h_vmem=40000M" in err and "-l h_rt=24:00:00" in err, err

_, err = submit("x", 4, 10000, 90, "standard.q", SGE_MEM_COMPLEX="mem_free")
assert "-l mem_free=2500M" in err, err
_, err = submit("x", 2, 10000, 45, "test.q", SGE_MEM_COMPLEX="none")
assert "h_vmem" not in err and "mem_free" not in err, err
# mem_free is not consumable on this cluster: co-requesting it reserves nothing,
# it only keeps the job off a node whose RAM is already committed.
_, err = submit("x", 4, 10000, 90, "standard.q", SGE_MEM_FREE_TOO="1")
assert "-l h_vmem=2500M" in err and "-l mem_free=2500M" in err, err
# a consumable=JOB complex must NOT be divided by the slot count
_, err = submit("x", 8, 32000, 90, "standard.q", SGE_MEM_PER_JOB="1")
assert "-l h_vmem=32000M" in err, err

def status(jid, **kw):
    r = subprocess.run(["bash", sta, str(jid)], capture_output=True, text=True,
                       env={**env, "SGE_STATUS_SLEEP": "0", **kw}, cwd=d)
    return r.stdout.strip()

fake("qstat", "exit 0\n"); fake("qacct", "exit 1\n")
assert status(1) == "running"
fake("qstat", "exit 1\n")
fake("qacct", "printf 'exit_status  0\\nfailed       0\\n'\n")
assert status(2) == "success"
fake("qacct", "printf 'exit_status  1\\nfailed       0\\n'\n")
assert status(3) == "failed"
fake("qacct", "printf 'exit_status  137\\nfailed       37 : h_rt\\n'\n")
assert status(4) == "failed"
# the qstat/qacct race: gone from one, not yet in the other. Reporting "failed"
# here would kill a healthy run.
fake("qacct", "exit 1\n")
assert status(5, SGE_STATUS_ATTEMPTS="3") == "running"
assert status(5, SGE_STATUS_ATTEMPTS="3") == "running"
assert status(5, SGE_STATUS_ATTEMPTS="3") == "failed"

# --- check_sge_profile.py against this cluster's real qconf output ------------
fake("qconf", '''
case "$1" in
  -sc) printf '#name shortcut type relop requestable consumable default urgency\\n'
       printf 'h_vmem h_vmem MEMORY <= YES YES 0 0\\n'
       printf 'mem_free mf MEMORY <= YES NO 0 0\\n' ;;
  -spl) printf 'openmpi\\nparallel\\n' ;;
  -sq) case "$2" in
         standard.q) echo "h_rt 24:00:00" ;;
         long.q) echo "h_rt 672:00:00" ;;
         *) exit 1 ;;
       esac ;;
esac
''')
fake("qstat", "printf 'cryo-em.q a\\nlong.q b\\nstandard.q c\\ntest.q d\\n'\n")

def check(profile, *extra, **kw):
    return subprocess.run([sys.executable, chk, profile, *extra], capture_output=True,
                          text=True, env={**env, **kw})

for p in profiles:
    r = check(p, "--no-qsub-verify")
    assert r.returncode == 0, f"{p} should validate:\n{r.stdout}{r.stderr}"
    assert "consumable per slot" in r.stdout, r.stdout

# requesting a non-consumable complex reserves nothing -> must be flagged
r = check(profiles[0], "--mem-complex", "mem_free", "--no-qsub-verify")
assert r.returncode == 1 and "NOT consumable" in r.stderr, r.stderr

# a runtime past, or flush against, a queue's h_rt ceiling must be flagged
import textwrap as _tw
bad = f"{d}/badprof"; os.makedirs(bad, exist_ok=True)
open(f"{bad}/config.yaml", "w").write(_tw.dedent("""
    sge-pe: smp
    default-resources: {runtime: 120, sge_queue: gpu.q}
    set-resources:
      iqtree_species: {runtime: 1500, sge_queue: standard.q}
      scaffold: {runtime: 1440, sge_queue: standard.q}
"""))
r = check(bad, "--no-qsub-verify")
assert r.returncode == 1, r.stdout
for expect in ("exceeds standard.q h_rt", "within 10% of standard.q",
               "'smp' does not exist", "gpu.q"):
    assert expect in r.stderr, f"{expect!r} not flagged:\n{r.stderr}"

# --- qconf is admin-host only. A REFUSAL to answer is not the answer "no". ----
# `qconf -sq standard.q` -> denied. The checker must fall back to cluster.yaml
# and must NOT report the queue as missing.
fake("qconf", '''echo 'denied: host "morty.eb.local" is not an admin host' >&2\nexit 1\n''')
r = check(profiles[0], "--no-qsub-verify")
assert r.returncode == 0, f"denied qconf must not fail the check:\n{r.stdout}{r.stderr}"
assert "not found" not in r.stderr and "does not exist" not in r.stderr, r.stderr
assert "cluster.yaml" in r.stdout, r.stdout

# ...but a qconf that ANSWERS "no such queue" is a real absence
fake("qconf", '''
case "$1" in
  -sc) printf 'h_vmem h_vmem MEMORY <= YES YES 0 0\\n' ;;
  -spl) printf 'parallel\\n' ;;
  -sq) [ "$2" = standard.q ] && { echo "h_rt 24:00:00"; exit 0; }
       [ "$2" = long.q ] && { echo "h_rt 672:00:00"; exit 0; }
       echo "error: queue \\"$2\\" does not exist" >&2; exit 1 ;;
esac
''')
gone = f"{d}/goneq"; os.makedirs(gone, exist_ok=True)
open(f"{gone}/config.yaml", "w").write(
    "sge-pe: parallel\ndefault-resources: {runtime: 60, mem_mb: 4000, sge_queue: nosuch.q}\n")
r = check(gone, "--facts", "/dev/null", "--no-qsub-verify")
assert r.returncode == 1 and "does not exist" in r.stderr, r.stderr

# stale recorded facts must be reported as drift, not silently believed
fake("qconf", '''
case "$1" in
  -sc) printf 'h_vmem h_vmem MEMORY <= YES JOB 0 0\\n' ;;
  -spl) printf 'openmpi\\n' ;;
  -sq) [ "$2" = standard.q ] && { echo "h_rt 12:00:00"; exit 0; }
       [ "$2" = long.q ] && { echo "h_rt 672:00:00"; exit 0; }
       exit 1 ;;
esac
''')
r = check(profiles[0], "--no-qsub-verify")
assert "DRIFT" in r.stdout and "consumable=JOB" in r.stdout, r.stdout

# --- `qsub -w v` works from any submit host. It exits 1 to say "no suitable
# queues", which is an ANSWER, not a refusal, and must be read as a rejection.
fake("qconf", '''echo 'denied: host "x" is not an admin host' >&2\nexit 1\n''')
fake("qsub", '''
q=""; hrt=""
while [ $# -gt 0 ]; do
  case "$1" in -q) q=$2; shift;; -l) case "$2" in h_rt=*) hrt=${2#h_rt=};; esac; shift;; esac
  shift
done
hh=${hrt%%:*}
if [ "$q" = standard.q ] && [ "${hh:-0}" -gt 24 ] 2>/dev/null; then
  echo "verification: no suitable queues"; exit 1
fi
echo "verification: found suitable queue(s)"
''')
r = check(profiles[0])
assert r.returncode == 0 and "accepted by `qsub -w v`" in r.stdout, r.stdout
over = f"{d}/overq"; os.makedirs(over, exist_ok=True)
open(f"{over}/config.yaml", "w").write(
    "sge-pe: parallel\n"
    "default-resources: {runtime: 60, mem_mb: 4000, sge_queue: standard.q}\n"
    "set-resources:\n  run_dbcan: {runtime: 1500, mem_mb: 16000, sge_queue: standard.q}\n")
r = check(over, "--facts", "/dev/null")     # no recorded ceilings: qsub alone must catch it
assert r.returncode == 1 and "qsub -w v` rejects" in r.stderr, f"{r.stdout}{r.stderr}"

# no scheduler on PATH at all: skip cleanly rather than fail the suite
r = subprocess.run([sys.executable, chk, profiles[0]], capture_output=True, text=True,
                   env={**os.environ, "PATH": "/usr/bin:/bin"})
assert r.returncode == 0, r.stdout + r.stderr
print("OK")
PY

hdr "4. statistical unit checks"
python3 tests/test_statistics.py || FAIL=1

hdr "5. end-to-end smoke tests"
python3 tests/test_pipeline_smoke.py || FAIL=1

printf '\n'
if [ "$FAIL" -ne 0 ]; then echo "SUITE FAILED"; exit 1; fi
echo "SUITE PASSED"
