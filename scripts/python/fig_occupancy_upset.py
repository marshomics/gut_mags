#!/usr/bin/env python3
"""
fig_occupancy_upset.py
----------------------
Figure 3: niche specificity.

  A  UpSet of species occupancy across niches (named species; placeholder
     clusters reported separately because they are not cross-database
     comparable). Shows how niche-specific species are and how small the overlap
     is.
  B  Standardised Levins' niche breadth distribution; most species sit near 0
     (specialists), quantifying the "niche-specific with minor overlap" claim.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--upset", required=True)
    ap.add_argument("--breadth", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    niches = cfg["inputs"]["niche_levels"]

    up = pd.read_csv(args.upset, sep="\t")
    # build a boolean-indexed Series of named-species counts for UpSet
    idx_rows, counts = [], []
    for _, r in up.iterrows():
        patt = str(r["occupancy_pattern"]).split("+") if str(r["occupancy_pattern"]) else []
        idx_rows.append(tuple(n in patt for n in niches))
        counts.append(int(r.get("n_species_named", 0)))
    mi = pd.MultiIndex.from_tuples(idx_rows, names=niches)
    series = pd.Series(counts, index=mi)
    series = series[series > 0]

    fig = plt.figure(figsize=(mm(180), mm(80)))
    # UpSet draws onto the current figure via a dict of axes
    try:
        UpSet(series, sort_by="cardinality", show_counts=True,
              element_size=None).plot(fig=fig)
        fig.suptitle("A  Species niche occupancy (named species)",
                     x=0.02, ha="left", fontweight="bold")
    except Exception as e:
        ax = fig.add_subplot(111); ax.text(0.5, 0.5, f"UpSet failed: {e}", ha="center")

    fig.savefig  # no-op guard
    save(fig, args.out, cfg)

    # B breadth as a separate small figure (kept independent for layout freedom)
    br = pd.read_csv(args.breadth, sep="\t")
    fig2, ax2 = plt.subplots(figsize=(mm(85), mm(60)))
    ax2.hist(br["levins_B_std"].dropna(), bins=40, color="#555555")
    ax2.axvline(0.25, ls="--", lw=0.6, color="#cc3311")
    ax2.set_xlabel("standardised niche breadth (0=specialist, 1=even)")
    ax2.set_ylabel("species")
    ax2.set_title("B  Niche breadth", loc="left", fontweight="bold")
    save(fig2, args.out + "_breadth", cfg)


if __name__ == "__main__":
    main()
