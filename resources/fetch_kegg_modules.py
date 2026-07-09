#!/usr/bin/env python3
"""
fetch_kegg_modules.py
---------------------
Download KEGG module definitions via the public KEGG REST API and write the
TSV (module_id, name, definition) consumed by kegg_module_completeness.py.

Run ONCE on a login node with internet access, then point the workflow at the
output. KEGG REST is rate-limited; this script sleeps between calls and is
resumable (skips modules already written).

Usage:
    python fetch_kegg_modules.py --out resources/kegg_modules.tsv
"""
import argparse
import time
import os
import urllib.request

LIST_URL = "https://rest.kegg.jp/list/module"
GET_URL = "https://rest.kegg.jp/get/{mod}"


def fetch(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode()


def parse_definition(entry_text):
    name, definition = "", ""
    for line in entry_text.splitlines():
        if line.startswith("NAME"):
            name = line.replace("NAME", "").strip()
        elif line.startswith("DEFINITION"):
            definition = line.replace("DEFINITION", "").strip()
    return name, definition


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sleep", type=float, default=0.34)
    args = ap.parse_args()

    done = set()
    if os.path.exists(args.out):
        with open(args.out) as fh:
            done = {l.split("\t")[0] for l in fh if l.strip() and not l.startswith("module_id")}

    modules = [l.split("\t")[0].replace("md:", "")
               for l in fetch(LIST_URL).splitlines() if l.strip()]
    mode = "a" if done else "w"
    with open(args.out, mode) as out:
        if not done:
            out.write("module_id\tname\tdefinition\n")
        for i, m in enumerate(modules, 1):
            if m in done:
                continue
            try:
                name, defn = parse_definition(fetch(GET_URL.format(mod=m)))
                out.write(f"{m}\t{name}\t{defn}\n")
                out.flush()
            except Exception as e:
                print(f"WARN {m}: {e}")
            time.sleep(args.sleep)
            if i % 50 == 0:
                print(f"... {i}/{len(modules)}")
    print(f"Wrote module definitions to {args.out}")


if __name__ == "__main__":
    main()
