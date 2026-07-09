#!/usr/bin/env python3
"""
fig_enrichment_dotplot.py
-------------------------
Enrichment dot plot for one contrast: the functional categories over-represented
among the niche-signature features. x = GSEA NES (sign shows enrichment among
positively- vs negatively-associated features), dot size = number of signature
features in the category, colour = -log10 q. Categories are grouped by system
(KEGG pathway, CAZyme class, ...). Both-method-significant categories are drawn
with a bold outline. PNG + editable-text SVG.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--enrichment", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    d = pd.read_csv(args.enrichment, sep="\t")
    d = d[d["contrast"] == args.contrast].copy()
    if d.empty:
        fig, ax = plt.subplots(figsize=(mm(120), mm(60)))
        ax.text(0.5, 0.5, f"no enrichment for {args.contrast}", ha="center")
        save(fig, args.out, cfg); return

    d["q"] = d[["ora_q", "gsea_padj"]].min(axis=1)
    d["nes"] = pd.to_numeric(d.get("NES"), errors="coerce")
    d["k"] = pd.to_numeric(d.get("k"), errors="coerce").fillna(
        pd.to_numeric(d.get("size"), errors="coerce")).fillna(5)
    # rank: both-significant first, then by q
    d["both"] = d["confidence"] == "both"
    d = d.sort_values(["both", "q"], ascending=[False, True]).head(args.top)
    d = d.sort_values(["system", "nes"])
    d["label"] = "[" + d["system"].str.replace("kegg_", "", regex=False) + "] " + \
                 d["category_name"].astype(str).str.slice(0, 42)

    y = np.arange(len(d))
    neglogq = -np.log10(d["q"].clip(lower=1e-300))
    sizes = 20 + 120 * (d["k"] / max(d["k"].max(), 1))
    fig, ax = plt.subplots(figsize=(mm(150), mm(max(60, 4 + 4 * len(d)))))
    sc = ax.scatter(d["nes"].fillna(0), y, s=sizes, c=neglogq, cmap="viridis",
                    edgecolors=np.where(d["both"], "black", "none"), linewidths=0.8)
    ax.axvline(0, lw=0.5, color="#888")
    ax.set_yticks(y); ax.set_yticklabels(d["label"], fontsize=5)
    ax.set_xlabel("GSEA NES (+ enriched among positively-associated features)")
    ax.set_title(f"Functional enrichment — {args.contrast}", loc="left", fontweight="bold")
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02); cb.set_label("-log10 q")
    # size legend
    for kk in sorted(set([int(d['k'].min()), int(d['k'].median()), int(d['k'].max())])):
        ax.scatter([], [], s=20 + 120 * (kk / max(d['k'].max(), 1)), c="#888",
                   label=f"{kk} genes")
    ax.legend(frameon=False, fontsize=5, loc="lower right", title="overlap")
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
