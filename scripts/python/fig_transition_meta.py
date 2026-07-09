#!/usr/bin/env python3
"""
fig_transition_meta.py
----------------------
Cross-species result: are recent niche acquisitions consistent in direction?
  A  number of species with a supported recent acquisition into each niche.
  B  directional counts per niche pair (source->derived both ways) with the
     binomial test p for a consistent direction.
  C  distribution of transition depths across species (shallow = recent).
This panel is the headline: repeated, independent, same-direction acquisitions.
"""
import argparse
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, niche_palette, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--directionality", required=True)
    ap.add_argument("--all-calls", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)

    dirn = pd.read_csv(args.directionality, sep="\t") if __import__("os").path.getsize(args.directionality) else pd.DataFrame()
    calls = pd.read_csv(args.all_calls, sep="\t") if __import__("os").path.exists(args.all_calls) else pd.DataFrame()
    summary = json.load(open(args.summary))

    fig, ax = plt.subplots(1, 3, figsize=(mm(180), mm(62)))

    # A acquisitions per derived niche
    counts = dict(summary.get("recent_acquisition_niche_counts", []))
    niches = [n for n in cfg["inputs"]["niche_levels"] if n in counts] or list(counts)
    ax[0].bar(range(len(niches)), [counts.get(n, 0) for n in niches],
              color=[pal.get(n, "#777") for n in niches])
    ax[0].set_xticks(range(len(niches))); ax[0].set_xticklabels(niches)
    ax[0].set_ylabel("species (supported)")
    ax[0].set_title("A  Recent acquisition by niche", loc="left", fontweight="bold")

    # B directional counts per pair
    if not dirn.empty:
        labels, left, right, ps = [], [], [], []
        for _, r in dirn.iterrows():
            cols = [c for c in dirn.columns if "->" in c]
            if len(cols) >= 2:
                labels.append(r["pair"]); left.append(r[cols[0]]); right.append(r[cols[1]])
                ps.append(r.get("binomial_p", np.nan))
        y = np.arange(len(labels))
        ax[1].barh(y - 0.2, left, 0.4, color="#4477AA", label="dir 1")
        ax[1].barh(y + 0.2, right, 0.4, color="#CC6677", label="dir 2")
        ax[1].set_yticks(y); ax[1].set_yticklabels(labels, fontsize=5)
        ax[1].set_xlabel("species")
        for i, p in enumerate(ps):
            if pd.notna(p):
                ax[1].text(max(left[i], right[i]), i, f" p={p:.2g}", va="center", fontsize=5)
        ax[1].set_title("B  Direction per pair", loc="left", fontweight="bold")
        ax[1].legend(frameon=False, fontsize=5)
    else:
        ax[1].text(0.5, 0.5, "no directional calls", ha="center")

    # C transition depth distribution
    if not calls.empty and "transition_depth" in calls.columns:
        d = pd.to_numeric(calls["transition_depth"], errors="coerce").dropna()
        ax[2].hist(d, bins=20, color="#555555")
        ax[2].axvline(0.5, ls="--", lw=0.6, color="#cc3311")
        ax[2].set_xlabel("transition depth (0 recent .. 1 deep)")
        ax[2].set_ylabel("acquisitions")
        ax[2].set_title("C  Recency", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
