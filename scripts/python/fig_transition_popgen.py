#!/usr/bin/env python3
"""
fig_transition_popgen.py
------------------------
Per-species population-genetic evidence for a recent acquisition.
  A  within-niche nucleotide diversity (pi) with bootstrap CI.
  B  within-niche Tajima's D with CI (line at 0; negative = expansion).
  C  folded site-frequency spectra per niche (excess of rare variants in a
     recently founded niche).
All on strain-dereplicated, equal-n-subsampled populations.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, niche_palette, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--popgen-dir", required=True)
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)

    div = pd.read_csv(f"{args.popgen_dir}/popgen_diversity.tsv", sep="\t")
    niches = [n for n in cfg["inputs"]["niche_levels"] if n in set(div["niche"])]
    div = div.set_index("niche").reindex(niches)

    fig, ax = plt.subplots(1, 3, figsize=(mm(180), mm(62)))
    x = np.arange(len(niches))

    def err(mean, lo, hi):  # clip tiny float-noise negatives to 0
        return np.vstack([np.clip(mean - lo, 0, None), np.clip(hi - mean, 0, None)])

    ax[0].bar(x, div["pi_mean"], color=[pal[n] for n in niches],
              yerr=err(div["pi_mean"], div["pi_lo"], div["pi_hi"]), error_kw=dict(lw=0.6))
    ax[0].set_xticks(x); ax[0].set_xticklabels(niches)
    ax[0].set_ylabel("nucleotide diversity pi")
    ax[0].set_title("A  Diversity", loc="left", fontweight="bold")

    ax[1].bar(x, div["tajimaD_mean"], color=[pal[n] for n in niches],
              yerr=err(div["tajimaD_mean"], div["tajimaD_lo"], div["tajimaD_hi"]),
              error_kw=dict(lw=0.6))
    ax[1].axhline(0, lw=0.5, color="#333333")
    ax[1].set_xticks(x); ax[1].set_xticklabels(niches)
    ax[1].set_ylabel("Tajima's D")
    ax[1].set_title("B  Tajima's D", loc="left", fontweight="bold")

    for n in niches:
        f = f"{args.popgen_dir}/sfs_{n}.tsv"
        if os.path.exists(f):
            s = pd.read_csv(f, sep="\t")
            tot = s["mean_sites"].sum()
            if tot > 0:
                ax[2].plot(s["minor_allele_count"], s["mean_sites"] / tot, "-o",
                           color=pal[n], ms=2, label=n)
    ax[2].set_xlabel("minor-allele count"); ax[2].set_ylabel("proportion of SNPs")
    ax[2].set_title("C  Folded SFS", loc="left", fontweight="bold")
    ax[2].legend(frameon=False)
    fig.suptitle(args.species_id, x=0.02, ha="left", fontsize=7)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
