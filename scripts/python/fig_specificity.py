#!/usr/bin/env python3
"""
fig_specificity.py
------------------
Figure: how niche-specific taxa are.
  A  number of FDR-significant niche-specific taxa per rank, coloured by the
     niche they are specific to (from the permutation specificity test).
  B  standardised niche-breadth distribution at the family rank, coloured by
     dominant niche (most taxa near 0 = specialists).
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
    ap.add_argument("--spec-dir", required=True)
    ap.add_argument("--rank", default="family")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]

    summ = pd.read_csv(f"{args.spec_dir}/specialist_summary.tsv", sep="\t")
    rankf = f"{args.spec_dir}/taxon_specificity_{args.rank}.tsv"
    fig, ax = plt.subplots(1, 2, figsize=(mm(180), mm(75)))

    # A
    ranks = list(dict.fromkeys(summ["rank"]))
    x = np.arange(len(ranks)); w = 0.25
    for i, n in enumerate(niches):
        vals = [summ[(summ["rank"] == r) & (summ["niche"] == n)]
                ["n_significant_specific_taxa"].sum() for r in ranks]
        ax[0].bar(x + (i - 1) * w, vals, w, color=pal[n], label=n)
    ax[0].set_xticks(x); ax[0].set_xticklabels(ranks, rotation=30, ha="right")
    ax[0].set_ylabel("FDR-significant niche-specific taxa")
    ax[0].set_title("A  Niche-specific taxa per rank", loc="left", fontweight="bold")
    ax[0].legend(frameon=False)

    # B
    try:
        d = pd.read_csv(rankf, sep="\t")
        for n in niches:
            sub = d[d["dominant_niche"] == n]["levins_B_std"].dropna()
            ax[1].hist(sub, bins=30, alpha=0.6, color=pal[n], label=n)
        ax[1].axvline(cfg["taxonomy"]["specialist_breadth_cutoff"], ls="--", lw=0.6,
                      color="#333333")
        ax[1].set_xlabel("standardised niche breadth")
        ax[1].set_ylabel(f"{args.rank} taxa")
        ax[1].set_title(f"B  Breadth ({args.rank})", loc="left", fontweight="bold")
        ax[1].legend(frameon=False)
    except Exception as e:
        ax[1].text(0.5, 0.5, f"n/a\n{e}", ha="center")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
