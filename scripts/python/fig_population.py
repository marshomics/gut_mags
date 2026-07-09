#!/usr/bin/env python3
"""
fig_population.py
-----------------
Western vs non-Western: functional conservation despite species turnover.
  A  taxonomic vs functional Sorensen dissimilarity per layer (paired bars). The
     headline: taxonomic high, functional low = functions conserved across
     populations.
  B  carrier substitution: distribution of carrier-species Jaccard for functions
     shared by both populations (low = same function, different species).
  C  % shared species vs % shared functions.
Inert-safe (shows a note if the population analysis did not run). PNG + editable SVG.
"""
import argparse
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, niche_palette, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--turnover", required=True)
    ap.add_argument("--substitution-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    t = pd.read_csv(args.turnover, sep="\t")
    if "layer" not in t.columns:
        fig, ax = plt.subplots(figsize=(mm(120), mm(50)))
        ax.text(0.5, 0.5, "Western/non-Western analysis not run\n(add the western_nonwestern "
                "column and set population.enabled: true)", ha="center", fontsize=7)
        save(fig, args.out, cfg); return

    fig, ax = plt.subplots(1, 3, figsize=(mm(185), mm(62)))
    x = np.arange(len(t)); w = 0.38
    ax[0].bar(x - w / 2, t["taxonomic_sorensen"], w, color="#0072B2", label="taxonomic")
    ax[0].bar(x + w / 2, t["functional_sorensen"], w, color="#009E73", label="functional")
    ax[0].set_xticks(x); ax[0].set_xticklabels(t["layer"], rotation=20, ha="right")
    ax[0].set_ylabel("Sorensen dissimilarity"); ax[0].set_ylim(0, 1)
    ax[0].set_title("A  Taxonomic vs functional turnover", loc="left", fontweight="bold")
    ax[0].legend(frameon=False)

    # B carrier substitution (pool all layers' shared-function carrier jaccard)
    sub = []
    for f in glob.glob(f"{args.substitution_dir}/carrier_substitution_*.tsv"):
        try:
            sub.append(pd.read_csv(f, sep="\t"))
        except Exception:
            pass
    if sub:
        cj = pd.concat(sub, ignore_index=True)["carrier_jaccard"].dropna()
        ax[1].hist(cj, bins=30, color="#882255")
        ax[1].axvline(cj.median(), ls="--", lw=0.6, color="#333")
        ax[1].set_xlabel("carrier-species Jaccard (shared functions)")
        ax[1].set_ylabel("functions")
        ax[1].text(0.95, 0.92, f"median {cj.median():.2f}", transform=ax[1].transAxes,
                   ha="right", fontsize=6)
    else:
        ax[1].text(0.5, 0.5, "no shared-function carriers", ha="center")
    ax[1].set_title("B  Carrier substitution", loc="left", fontweight="bold")

    # C % shared
    ax[2].bar(x - w / 2, t["pct_shared_species"], w, color="#0072B2", label="species")
    ax[2].bar(x + w / 2, t["pct_shared_functions"], w, color="#009E73", label="functions")
    ax[2].set_xticks(x); ax[2].set_xticklabels(t["layer"], rotation=20, ha="right")
    ax[2].set_ylabel("% shared (W & non-W)")
    ax[2].set_title("C  Shared species vs functions", loc="left", fontweight="bold")
    ax[2].legend(frameon=False)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
