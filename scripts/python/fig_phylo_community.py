#!/usr/bin/env python3
"""
fig_phylo_community.py
----------------------
Figure: phylogenetic structure of each niche.
  A  Faith's PD standardised effect size (SES) per niche.
  B  NRI (deep clustering) and NTI (tip clustering) per niche.
Positive SES/NRI/NTI = phylogenetically clustered (a restricted slice of the
tree); negative = overdispersed. Significant bars are outlined. This is the
phylogenetic statement of how "special" each niche's membership is.
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
    ap.add_argument("--phylo-community", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    d = pd.read_csv(args.phylo_community, sep="\t")
    niches = [n for n in cfg["inputs"]["niche_levels"] if n in set(d["niche"])]
    d = d.set_index("niche").reindex(niches)

    fig, ax = plt.subplots(1, 2, figsize=(mm(170), mm(70)))
    # A PD SES
    ax[0].bar(range(len(d)), d["PD_ses_z"], color=[pal[n] for n in d.index],
              edgecolor=["black" if p < 0.05 else "none" for p in d["PD_p"]], linewidth=1)
    ax[0].axhline(0, lw=0.5, color="#333333")
    ax[0].set_xticks(range(len(d))); ax[0].set_xticklabels(d.index)
    ax[0].set_ylabel("Faith's PD SES (z)")
    ax[0].set_title("A  Phylogenetic diversity vs null", loc="left", fontweight="bold")

    # B NRI/NTI
    x = np.arange(len(d)); w = 0.38
    ax[1].bar(x - w / 2, d["NRI"], w, label="NRI (deep)", color="#4477AA",
              edgecolor=["black" if p < 0.05 else "none" for p in d["MPD_p"]], linewidth=1)
    ax[1].bar(x + w / 2, d["NTI"], w, label="NTI (tips)", color="#CC6677",
              edgecolor=["black" if p < 0.05 else "none" for p in d["MNTD_p"]], linewidth=1)
    ax[1].axhline(0, lw=0.5, color="#333333")
    ax[1].set_xticks(x); ax[1].set_xticklabels(d.index)
    ax[1].set_ylabel("clustering index")
    ax[1].set_title("B  NRI / NTI (outlined = p<0.05)", loc="left", fontweight="bold")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
