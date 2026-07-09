#!/usr/bin/env python3
"""
compare_sensitivity.py
----------------------
Compare consensus signature calls between the main run and each sensitivity run.
For every layer x contrast it reports the overlap (Jaccard) of consensus
features and how many calls are gained/lost, so the report can state how stable
the signature is. High Jaccard against prev01/prev09/hqonly/nomouse = robust;
a large jump in the nophylo run quantifies how much phylogeny correction matters.
"""
import argparse
import glob
import os

import pandas as pd


def consensus_set(run_dir, layer, contrast):
    path = f"{run_dir}/05_diff/{layer}/{contrast}_signatures.tsv"
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, sep="\t")
    d["consensus_signature"] = d["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    return set(d.loc[d["consensus_signature"], "feature"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # discover layer/contrast pairs from base
    pairs = []
    for p in glob.glob(f"{a.base}/05_diff/*/*_signatures.tsv"):
        layer = os.path.basename(os.path.dirname(p))
        contrast = os.path.basename(p).replace("_signatures.tsv", "")
        pairs.append((layer, contrast))
    rows = []
    for run in a.runs:
        tag = os.path.basename(run).replace("results_", "")
        for layer, contrast in sorted(set(pairs)):
            b = consensus_set(a.base, layer, contrast)
            r = consensus_set(run, layer, contrast)
            if b is None or r is None:
                continue
            inter = len(b & r); union = len(b | r)
            rows.append({"run": tag, "layer": layer, "contrast": contrast,
                         "n_base": len(b), "n_run": len(r), "shared": inter,
                         "lost": len(b - r), "gained": len(r - b),
                         "jaccard": round(inter / union, 4) if union else 1.0})
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"Wrote {a.out}")


if __name__ == "__main__":
    main()
