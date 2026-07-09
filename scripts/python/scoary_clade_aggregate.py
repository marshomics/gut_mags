#!/usr/bin/env python3
"""
scoary_clade_aggregate.py
-------------------------
Combine clade-stratified Scoary2 results into one table and summary. Because
ortholog-family ids are clade-specific, recurrence is summarised at the level of
how many clades show niche-associated families per contrast (and per niche
direction), rather than by shared family id.

Inputs: parsed per-clade Scoary tables at <root>/clades/<clade_id>/scoary/<contrast>.tsv
Outputs: clade_scoary_all.tsv, clade_scoary_summary.tsv
"""
import argparse
import glob
import os

import pandas as pd

from hgn_utils import get_logger

log = get_logger("scoary-clade-agg")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="clade results root (…/clade)")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    files = glob.glob(f"{args.root}/clades/*/scoary/*.tsv")
    frames = []
    for f in files:
        try:
            d = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if d.empty:
            continue
        d["clade_id"] = os.path.basename(os.path.dirname(os.path.dirname(f)))
        d["contrast"] = os.path.splitext(os.path.basename(f))[0]
        frames.append(d)
    if not frames:
        pd.DataFrame().to_csv(f"{args.out_prefix}_all.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False)
        log.warning("No clade Scoary results found."); return
    alld = pd.concat(frames, ignore_index=True)
    alld.to_csv(f"{args.out_prefix}_all.tsv", sep="\t", index=False)

    sig = alld[alld["scoary_sig"].astype(str).isin(["True", "TRUE", "1"])]
    summ = (sig.groupby(["contrast", "clade_id"])
            .agg(n_sig_families=("feature", "nunique")).reset_index())
    per_contrast = (summ.groupby("contrast")
                    .agg(n_clades_with_hits=("clade_id", "nunique"),
                         total_sig_families=("n_sig_families", "sum")).reset_index())
    per_contrast.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False)
    log.info("Clade Scoary aggregate: %d clade-contrast results, %s",
             len(summ), per_contrast.to_dict("records"))


if __name__ == "__main__":
    main()
