#!/usr/bin/env python3
"""
enrichment_aggregate.py
-----------------------
Concatenate the per-layer/per-contrast combined enrichment tables into one, and
summarise the functional themes that define each niche (categories enriched by
both methods), for the report and figures.

Output: enrichment_all.tsv, enrichment_top_by_contrast.tsv
"""
import argparse
import glob

import pandas as pd

from hgn_utils import get_logger

log = get_logger("enrich-agg")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    frames = []
    for f in args.inputs:
        try:
            d = pd.read_csv(f, sep="\t")
            if len(d):
                frames.append(d)
        except Exception:
            continue
    if not frames:
        pd.DataFrame().to_csv(f"{args.out_prefix}_all.tsv", sep="\t", index=False)
        pd.DataFrame().to_csv(f"{args.out_prefix}_top_by_contrast.tsv", sep="\t", index=False)
        log.warning("No enrichment inputs."); return
    alld = pd.concat(frames, ignore_index=True)
    alld.to_csv(f"{args.out_prefix}_all.tsv", sep="\t", index=False)

    both = alld[alld["confidence"] == "both"].copy()
    # select only the columns that exist (GSEA/ORA columns can be absent if one
    # method produced nothing for a layer x contrast)
    want = ["contrast", "system", "category_id", "category_name", "layer",
            "direction", "fold_enrichment", "NES", "ora_q", "gsea_padj"]
    if len(both):
        sort_col = "ora_q" if "ora_q" in both.columns else "gsea_padj"
        grp = [c for c in ["contrast", "system", "category_id", "category_name"]
               if c in both.columns]
        top = (both.sort_values(sort_col).groupby(grp, as_index=False).first())
        top = top[[c for c in want if c in top.columns]]
        top.sort_values([c for c in ["contrast", "system", sort_col] if c in top.columns]).to_csv(
            f"{args.out_prefix}_top_by_contrast.tsv", sep="\t", index=False)
    else:
        pd.DataFrame(columns=want).to_csv(f"{args.out_prefix}_top_by_contrast.tsv",
                                          sep="\t", index=False)
    log.info("Enrichment aggregate: %d rows, %d both-significant across %d contrasts",
             len(alld), len(both), alld["contrast"].nunique())


if __name__ == "__main__":
    main()
