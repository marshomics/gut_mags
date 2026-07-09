#!/usr/bin/env python3
"""
fig_enrichment_heatmap.py
-------------------------
Functional themes across niches: a heatmap of category (rows) by contrast
(columns), coloured by signed enrichment (GSEA NES, or signed -log10 q where NES
is absent). Shows at a glance which functional categories distinguish the human
gut and how they compare with the animal and free-living contrasts. Top
both-method-significant categories are shown. PNG + editable-text SVG.
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
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=40)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    d = pd.read_csv(args.enrichment, sep="\t")
    if d.empty:
        fig, ax = plt.subplots(); ax.text(0.5, 0.5, "no enrichment", ha="center")
        save(fig, args.out, cfg); return

    d["q"] = d[["ora_q", "gsea_padj"]].min(axis=1)
    d["nes"] = pd.to_numeric(d.get("NES"), errors="coerce")
    d["score"] = d["nes"].fillna(
        np.sign((d["direction"] == "up").map({True: 1, False: -1})) *
        -np.log10(d["q"].clip(lower=1e-300)))
    d["cat"] = "[" + d["system"].str.replace("kegg_", "", regex=False) + "] " + \
               d["category_name"].astype(str).str.slice(0, 40)
    # pick top categories by best significance across contrasts
    best = d.groupby("cat")["q"].min().sort_values().head(args.top).index
    sub = d[d["cat"].isin(best)]
    mat = sub.pivot_table(index="cat", columns="contrast", values="score", aggfunc="mean")
    contrasts = [c for c in cfg["enrichment"]["contrasts"] if c in mat.columns]
    mat = mat[contrasts] if contrasts else mat

    vmax = np.nanpercentile(np.abs(mat.values), 98) or 1
    fig, ax = plt.subplots(figsize=(mm(40 + 16 * mat.shape[1]),
                                    mm(max(60, 4 + 4.2 * mat.shape[0]))))
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=35, ha="right")
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("signed enrichment (NES)")
    ax.set_title("Functional themes by niche contrast", loc="left", fontweight="bold")
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
