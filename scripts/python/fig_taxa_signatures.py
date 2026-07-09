#!/usr/bin/env python3
"""
fig_taxa_signatures.py
----------------------
Figure: which taxa characterise each niche.
  A  enrichment heatmap at a chosen rank (default family): log2 odds ratio of
     each top taxon across niches, FDR-significant cells marked. Shows the taxa
     that distinguish the human gut.
  B  top indicator taxa (IndVal.g) for the focal niche, ranked by indicator
     statistic, from the indicspecies output.
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
    ap.add_argument("--enrichment", required=True)     # enrichment_<rank>.tsv
    ap.add_argument("--indicator", required=True)       # indicator_<rank>.tsv
    ap.add_argument("--rank", default="family")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]
    focal = cfg["inputs"]["focal_niche"]

    en = pd.read_csv(args.enrichment, sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(mm(185), mm(95)),
                           gridspec_kw={"width_ratios": [1.2, 1]})

    # A enrichment heatmap (top taxa by max |log2 OR| among significant)
    sig = en[en["q"] < 0.05]
    top_taxa = (sig.assign(absor=sig["log2_or"].abs())
                .sort_values("absor", ascending=False)["taxon"].drop_duplicates().head(25).tolist())
    piv = en[en["taxon"].isin(top_taxa)].pivot_table(index="taxon", columns="niche",
                                                     values="log2_or")[niches]
    qpiv = en[en["taxon"].isin(top_taxa)].pivot_table(index="taxon", columns="niche",
                                                      values="q")[niches]
    piv = piv.reindex(top_taxa)
    vmax = np.nanpercentile(np.abs(piv.values), 98) or 1
    im = ax[0].imshow(piv.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax[0].set_xticks(range(len(niches))); ax[0].set_xticklabels(niches)
    ax[0].set_yticks(range(len(piv))); ax[0].set_yticklabels(
        [t.replace("f__", "") for t in piv.index], fontsize=5)
    for i in range(len(piv)):
        for j, n in enumerate(niches):
            if qpiv.iloc[i][n] < 0.05:
                ax[0].text(j, i, "*", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax[0], fraction=0.04, pad=0.02, label="log2 OR")
    ax[0].set_title(f"A  {args.rank} enrichment (* q<0.05)", loc="left", fontweight="bold")

    # B top indicators for focal niche
    try:
        ind = pd.read_csv(args.indicator, sep="\t")
        f = ind[ind["niche_combination"] == focal].sort_values("stat", ascending=False).head(15)
        ax[1].barh(range(len(f))[::-1], f["stat"], color=pal[focal])
        ax[1].set_yticks(range(len(f))[::-1])
        ax[1].set_yticklabels([t.replace("f__", "").replace("g__", "") for t in f["taxon"]],
                              fontsize=5)
        ax[1].set_xlabel("IndVal.g statistic")
        ax[1].set_title(f"B  Top {focal} indicators", loc="left", fontweight="bold")
    except Exception as e:
        ax[1].text(0.5, 0.5, f"indicator n/a\n{e}", ha="center", fontsize=6)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
