#!/usr/bin/env python3
"""
transition_meta.py
------------------
Aggregate the per-species verdicts into the cross-species result, which is the
strong claim: repeated, independent niche acquisitions in a consistent direction
across many species are far more convincing than any single genome.

For each ordered niche pair (source -> derived) it counts the species showing a
strong/moderate recent acquisition in that direction, and tests whether the
direction is consistent with a two-sided binomial (sign) test against 1:1. It
also summarises which niche is most often the recent acquisition and the
distribution of transition depths.

Inputs: all per-species transition_verdict.tsv files.
Outputs: meta_directionality.tsv, meta_summary.json
"""
import argparse
import glob
import json
import os
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from hgn_utils import load_config, get_logger

log = get_logger("trans-meta")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--verdict-glob", required=True,
                    help="glob for per-species transition_verdict.tsv")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    files = glob.glob(args.verdict_glob)
    frames = [pd.read_csv(f, sep="\t") for f in files if os.path.getsize(f)]
    frames = [f for f in frames if "derived_niche" in f.columns and len(f)]
    if not frames:
        json.dump({"n_species": 0}, open(f"{args.out_prefix}_summary.json", "w"))
        pd.DataFrame().to_csv(f"{args.out_prefix}_directionality.tsv", sep="\t"); return
    allv = pd.concat(frames, ignore_index=True)
    allv.to_csv(f"{args.out_prefix}_all_calls.tsv", sep="\t", index=False)

    # keep supported calls
    strong = allv[allv["tier"].isin(["strong", "moderate"])].copy()
    strong["direction"] = strong["ancestral_niche"] + "->" + strong["derived_niche"]

    # counts per ordered direction, and binomial test per unordered pair
    dir_counts = strong["direction"].value_counts().to_dict()
    rows = []
    pairs = set()
    for d in dir_counts:
        a, b = d.split("->")
        pairs.add(frozenset((a, b)))
    for pr in pairs:
        a, b = sorted(pr)
        ab = dir_counts.get(f"{a}->{b}", 0)
        ba = dir_counts.get(f"{b}->{a}", 0)
        n = ab + ba
        p = binomtest(ab, n, 0.5).pvalue if n > 0 else np.nan
        major = f"{a}->{b}" if ab >= ba else f"{b}->{a}"
        rows.append({"pair": f"{a}|{b}", f"{a}->{b}": ab, f"{b}->{a}": ba,
                     "n_species": n, "majority_direction": major,
                     "binomial_p": p})
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}_directionality.tsv", sep="\t", index=False)

    summary = {
        "n_species_analysed": int(allv["species_id"].nunique()),
        "n_supported_acquisitions": int(len(strong)),
        "tier_counts": allv["tier"].value_counts().to_dict(),
        "recent_acquisition_niche_counts": Counter(strong["derived_niche"]).most_common(),
        "direction_counts": dir_counts,
        "median_transition_depth": float(pd.to_numeric(strong["transition_depth"],
                                                       errors="coerce").median()),
    }
    json.dump(summary, open(f"{args.out_prefix}_summary.json", "w"), indent=2, default=str)
    log.info("Meta: %d species, supported acquisitions %d, directions %s",
             summary["n_species_analysed"], summary["n_supported_acquisitions"], dir_counts)


if __name__ == "__main__":
    main()
