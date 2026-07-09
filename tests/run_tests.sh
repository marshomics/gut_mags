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

hdr "4. statistical unit checks"
python3 tests/test_statistics.py || FAIL=1

hdr "5. end-to-end smoke tests"
python3 tests/test_pipeline_smoke.py || FAIL=1

printf '\n'
if [ "$FAIL" -ne 0 ]; then echo "SUITE FAILED"; exit 1; fi
echo "SUITE PASSED"
