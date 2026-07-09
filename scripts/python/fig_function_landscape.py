#!/usr/bin/env python3
"""
fig_function_landscape.py
-------------------------
Figure: specialised functional load across niches, per species. Violin/box of
CAZyme, BGC and AMR richness by niche (each species one observation). The
phylogenetic regression (PGLS) result for the niche term is annotated, so the
reader sees both the raw distribution and the phylogeny-controlled test in one
place rather than trusting the boxplot alone.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, niche_palette, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--traits", required=True)
    ap.add_argument("--pgls", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]

    tr = pd.read_csv(args.traits, sep="\t")
    try:
        pgls = pd.read_csv(args.pgls, sep="\t")
    except Exception:
        pgls = pd.DataFrame(columns=["trait", "term", "p"])

    layers = [("cazyme_richness", "CAZyme families"),
              ("bgc_richness", "BGC regions"),
              ("amr_richness", "AMR genes")]
    fig, ax = plt.subplots(1, len(layers), figsize=(mm(180), mm(70)))

    for j, (col, title) in enumerate(layers):
        if col not in tr.columns:
            ax[j].text(0.5, 0.5, f"{col}\nnot available", ha="center"); continue
        data = [tr.loc[tr["niche_primary"] == n, col].dropna().values for n in niches]
        parts = ax[j].violinplot(data, showmedians=True, widths=0.8)
        for i, b in enumerate(parts["bodies"]):
            b.set_facecolor(pal[niches[i]]); b.set_alpha(0.6); b.set_edgecolor("none")
        for key in ("cmedians", "cbars", "cmins", "cmaxes"):
            if key in parts:
                parts[key].set_color("#333333"); parts[key].set_linewidth(0.6)
        ax[j].set_xticks(range(1, len(niches) + 1)); ax[j].set_xticklabels(niches)
        ax[j].set_ylabel("count per species"); ax[j].set_title(title, loc="left", fontweight="bold")
        # annotate PGLS niche-term p (human contrast)
        sub = pgls[(pgls["trait"] == col) & (pgls["term"].astype(str).str.contains("niche"))]
        if len(sub):
            p = pd.to_numeric(sub["p"], errors="coerce").min()
            star = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
            ax[j].text(0.97, 0.95, f"PGLS niche {star}\n(p={p:.2g})",
                       transform=ax[j].transAxes, ha="right", va="top", fontsize=5)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
