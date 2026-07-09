#!/usr/bin/env python3
"""
fig_ordination_varpart.py
-------------------------
Figure: multivariate functional structure and what explains it.

  A  PCoA of species functional repertoires (Jaccard), coloured by niche, with
     the PERMANOVA niche R2 and p annotated.
  B  Variation partitioning: the functional variance uniquely attributable to
     niche, to phylogeny, to genome size and to quality, plus shared fractions.
     The niche-unique bar is the headline number: functional differentiation
     that is NOT explained by phylogeny.
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
    ap.add_argument("--coords", required=True)
    ap.add_argument("--varexp", required=True)
    ap.add_argument("--permanova", required=True)
    ap.add_argument("--varpart", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)

    coords = pd.read_csv(args.coords, sep="\t")
    try:
        ve = pd.read_csv(args.varexp, sep="\t", header=None)
        ve = dict(zip(ve[0], ve[1]))
    except Exception:
        ve = {}

    fig, ax = plt.subplots(1, 2, figsize=(mm(180), mm(80)))

    # A PCoA
    for n, d in coords.groupby("niche"):
        ax[0].scatter(d["PCo1"], d["PCo2"], s=4, alpha=0.5, color=pal.get(n, "#777"),
                      label=n, rasterized=True)
    ax[0].set_xlabel(f"PCo1 ({ve.get('PCo1','')}%)")
    ax[0].set_ylabel(f"PCo2 ({ve.get('PCo2','')}%)")
    ax[0].set_title("A  Functional PCoA", loc="left", fontweight="bold")
    ax[0].legend(frameon=False, markerscale=2)
    try:
        perm = pd.read_csv(args.permanova, sep="\t", index_col=0)
        r2 = perm.loc["niche", "R2"] if "niche" in perm.index else np.nan
        pv = perm.loc["niche", "Pr(>F)"] if "niche" in perm.index else np.nan
        ax[0].text(0.03, 0.95, f"PERMANOVA niche R²={r2:.3f}, p={pv:.3g}",
                   transform=ax[0].transAxes, fontsize=5, va="top")
    except Exception:
        pass

    # B varpart
    try:
        vp = pd.read_csv(args.varpart, sep="\t")
        # vegan indfract: rows like [a]..[o]; keep individual fractions with labels
        vp = vp.dropna(subset=["Adj.R.square"]) if "Adj.R.square" in vp.columns else vp
        labels = vp.iloc[:, 0].astype(str).tolist()
        vals = pd.to_numeric(vp.get("Adj.R.square", vp.iloc[:, -1]), errors="coerce").fillna(0)
        vals = vals.clip(lower=0)
        ax[1].barh(range(len(vals)), vals.values, color="#4477AA")
        ax[1].set_yticks(range(len(labels)))
        ax[1].set_yticklabels(labels, fontsize=5)
        ax[1].set_xlabel("adj. R² (functional variance)")
        ax[1].set_title("B  Variation partitioning", loc="left", fontweight="bold")
    except Exception as e:
        ax[1].text(0.5, 0.5, f"varpart unavailable\n{e}", ha="center", fontsize=6)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
