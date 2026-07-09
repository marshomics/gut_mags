#!/usr/bin/env python3
"""
fetch_kegg_genesets.py
----------------------
Fetch the canonical KEGG KO -> pathway membership used by the functional
enrichment step, from the public KEGG REST API. Run once on a node with internet
access. Writes resources/genesets/kegg_pathway.tsv (feature=KO, category_id,
category_name).

KO -> GO and BRITE are not produced here (KEGG REST does not expose them cleanly);
to test those systems, drop a map file (feature<TAB>category_id<TAB>name) into the
same directory as ko_go.tsv / brite.tsv and add the system to enrichment.systems.

Usage: python fetch_kegg_genesets.py --out-dir resources/genesets
"""
import argparse
import os
import urllib.request


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="resources/genesets")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # pathway names (map ids)
    names = {}
    for line in fetch("https://rest.kegg.jp/list/pathway").splitlines():
        if "\t" in line:
            pid, name = line.split("\t", 1)
            names[pid.replace("map", "")] = name      # store numeric suffix

    # KO -> pathway links (keep ko maps, i.e. path:koXXXXX)
    rows = []
    for line in fetch("https://rest.kegg.jp/link/pathway/ko").splitlines():
        if "\t" not in line:
            continue
        ko, path = line.split("\t")
        ko = ko.replace("ko:", "")
        if path.startswith("path:ko"):
            num = path.replace("path:ko", "")
            rows.append((ko, f"map{num}", names.get(num, f"map{num}")))

    out = os.path.join(args.out_dir, "kegg_pathway.tsv")
    with open(out, "w") as fh:
        fh.write("feature\tcategory_id\tcategory_name\n")
        for ko, cid, cname in rows:
            fh.write(f"{ko}\t{cid}\t{cname}\n")
    print(f"Wrote {len(rows)} KO-pathway memberships to {out}")


if __name__ == "__main__":
    main()
