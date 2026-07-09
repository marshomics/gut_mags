#!/usr/bin/env python3
"""
enrichment_combine.py
---------------------
Merge the over-representation (ORA) and preranked-GSEA results for one layer x
contrast into one table, and flag agreement. A category called enriched by both
an over-representation test and a threshold-free GSEA is the defensible call;
single-method hits are reported but tiered lower.

Output: enrichment_combined_<layer>_<contrast>.tsv
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("enrich-combine")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--ora", required=True)
    ap.add_argument("--gsea", required=True)
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    alpha = cfg["enrichment"]["fdr_alpha"]

    ora = pd.read_csv(args.ora, sep="\t")
    gsea = pd.read_csv(args.gsea, sep="\t")

    # ORA: best (most significant) direction per category
    if len(ora):
        ora_best = (ora.sort_values("q").groupby(["system", "category_id"], as_index=False)
                    .first()[["system", "category_id", "category_name", "direction",
                              "fold_enrichment", "k", "K", "q", "overlap"]]
                    .rename(columns={"q": "ora_q", "direction": "ora_direction"}))
    else:
        ora_best = pd.DataFrame(columns=["system", "category_id", "category_name",
                                         "ora_direction", "fold_enrichment", "k", "K",
                                         "ora_q", "overlap"])

    g = gsea[["system", "category_id", "category_name", "NES", "padj", "size",
              "leadingEdge"]].rename(columns={"padj": "gsea_padj"}) if len(gsea) else \
        pd.DataFrame(columns=["system", "category_id", "category_name", "NES",
                              "gsea_padj", "size", "leadingEdge"])

    m = ora_best.merge(g, on=["system", "category_id"], how="outer",
                       suffixes=("", "_g"))
    m["category_name"] = m["category_name"].fillna(m.get("category_name_g"))
    m = m.drop(columns=[c for c in m.columns if c.endswith("_g")], errors="ignore")
    m["layer"] = args.layer
    m["contrast"] = args.contrast
    m["ora_sig"] = m["ora_q"] < alpha
    m["gsea_sig"] = m["gsea_padj"] < alpha
    m["confidence"] = np.where(m["ora_sig"].fillna(False) & m["gsea_sig"].fillna(False), "both",
                       np.where(m["ora_sig"].fillna(False), "ora_only",
                       np.where(m["gsea_sig"].fillna(False), "gsea_only", "ns")))
    # direction: ORA direction if present else sign(NES)
    nes_dir = pd.Series(np.where(m["NES"] > 0, "up", "down"), index=m.index)
    m["direction"] = m["ora_direction"].fillna(nes_dir)
    m.sort_values(["confidence", "ora_q"]).to_csv(args.out, sep="\t", index=False)
    log.info("%s/%s combined: %d categories (%d both, %d any)",
             args.layer, args.contrast, len(m),
             int((m["confidence"] == "both").sum()),
             int((m["confidence"] != "ns").sum()))


if __name__ == "__main__":
    main()
