#!/usr/bin/env python3
"""
fig_community.py
----------------
Community-level metabolic interaction across niches (rarefied to equal species).
  A  auxotrophy dependency: mean auxotrophies per species by niche (CI).
  B  cross-feeding potential and number of active byproduct exchanges by niche.
  C  KEGG-module complementarity (collective coverage / mean per species).
  D  community trait composition (anaerobe / spore / bile / motile fractions).
A second figure shows the per-amino-acid auxotrophy and per-metabolite cross-
feeding heatmaps. PNG + editable-text SVG.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, niche_palette, save, mm


def load(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def err(df, col):
    return [np.clip(df[col] - df[f"{col}_lo"], 0, None),
            np.clip(df[f"{col}_hi"] - df[col], 0, None)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--community-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    d = args.community_dir
    summ = load(f"{d}/community_summary.tsv")
    aa = load(f"{d}/auxotrophy_by_aa.tsv")
    cf = load(f"{d}/crossfeeding_by_metabolite.tsv")
    tr = load(f"{d}/trait_composition.tsv")
    niches = [n for n in cfg["inputs"]["niche_levels"] if n in set(summ.get("niche", []))]
    summ = summ.set_index("niche").reindex(niches).reset_index() if len(summ) else summ

    fig, ax = plt.subplots(2, 2, figsize=(mm(185), mm(150)))
    x = np.arange(len(niches))

    if len(summ):
        ax[0, 0].bar(x, summ["mean_aux"], color=[pal[n] for n in niches],
                     yerr=err(summ, "mean_aux"), error_kw=dict(lw=0.6))
        ax[0, 0].set_xticks(x); ax[0, 0].set_xticklabels(niches)
        ax[0, 0].set_ylabel("mean auxotrophies / species")
    ax[0, 0].set_title("A  Auxotrophy dependency", loc="left", fontweight="bold")

    if len(summ):
        ax[0, 1].bar(x, summ["crossfeed_potential"], color=[pal[n] for n in niches],
                     yerr=err(summ, "crossfeed_potential"), error_kw=dict(lw=0.6))
        ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(niches)
        ax[0, 1].set_ylabel("cross-feeding potential (norm)")
    ax[0, 1].set_title("B  Byproduct cross-feeding", loc="left", fontweight="bold")

    if len(summ):
        ax[1, 0].bar(x, summ["mod_complementarity"], color=[pal[n] for n in niches],
                     yerr=err(summ, "mod_complementarity"), error_kw=dict(lw=0.6))
        ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels(niches)
        ax[1, 0].set_ylabel("module complementarity (collective / per-species)")
    ax[1, 0].set_title("C  Metabolic division of labour", loc="left", fontweight="bold")

    if len(tr):
        key = [t for t in cfg["community"]["key_traits"] if t in set(tr["trait"])]
        w = 0.8 / max(len(niches), 1)
        for i, n in enumerate(niches):
            sub = tr[tr["niche"] == n].set_index("trait").reindex(key)
            ax[1, 1].bar(np.arange(len(key)) + (i - (len(niches) - 1) / 2) * w,
                         sub["fraction"], w, color=pal[n], label=n)
        ax[1, 1].set_xticks(range(len(key)))
        ax[1, 1].set_xticklabels(key, rotation=30, ha="right", fontsize=5)
        ax[1, 1].set_ylabel("fraction of community")
        ax[1, 1].legend(frameon=False, fontsize=5)
    ax[1, 1].set_title("D  Community trait composition", loc="left", fontweight="bold")

    fig.suptitle("Human gut as an interacting community", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, args.out, cfg)

    # --- second figure: heatmaps ---
    fig2, ax2 = plt.subplots(1, 2, figsize=(mm(185), mm(90)))
    if len(aa):
        piv = aa.pivot_table(index="amino_acid", columns="niche", values="frac_auxotroph")[niches]
        im = ax2[0].imshow(piv.values, aspect="auto", cmap="magma_r", vmin=0, vmax=1)
        ax2[0].set_xticks(range(len(niches))); ax2[0].set_xticklabels(niches)
        ax2[0].set_yticks(range(len(piv))); ax2[0].set_yticklabels(piv.index, fontsize=5)
        fig2.colorbar(im, ax=ax2[0], fraction=0.04, pad=0.02, label="fraction auxotroph")
        ax2[0].set_title("Auxotrophy by amino acid", loc="left", fontweight="bold")
    if len(cf):
        piv = cf.pivot_table(index="metabolite", columns="niche", values="potential_norm")[niches]
        im = ax2[1].imshow(piv.values, aspect="auto", cmap="viridis")
        ax2[1].set_xticks(range(len(niches))); ax2[1].set_xticklabels(niches)
        ax2[1].set_yticks(range(len(piv))); ax2[1].set_yticklabels(piv.index, fontsize=5)
        fig2.colorbar(im, ax=ax2[1], fraction=0.04, pad=0.02, label="cross-feed potential")
        ax2[1].set_title("Cross-feeding by metabolite", loc="left", fontweight="bold")
    fig2.tight_layout()
    save(fig2, args.out + "_heatmaps", cfg)


if __name__ == "__main__":
    main()
