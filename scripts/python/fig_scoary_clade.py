#!/usr/bin/env python3
"""
fig_scoary_clade.py
-------------------
Clade-stratified Scoary2 summary: number of genera with niche-associated ortholog
families per contrast, and the total number of significant families. Shows
whether niche-associated gene families recur across independent genera.
"""
import argparse

import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    try:
        s = pd.read_csv(args.summary, sep="\t")
    except Exception:
        s = pd.DataFrame()

    fig, ax = plt.subplots(1, 2, figsize=(mm(160), mm(62)))
    if s.empty or "contrast" not in s.columns:
        ax[0].text(0.5, 0.5, "no clade results", ha="center")
        save(fig, args.out, cfg); return
    ax[0].bar(range(len(s)), s["n_clades_with_hits"], color="#117733")
    ax[0].set_xticks(range(len(s))); ax[0].set_xticklabels(s["contrast"], rotation=30, ha="right")
    ax[0].set_ylabel("genera with hits")
    ax[0].set_title("A  Clades with niche-associated families", loc="left", fontweight="bold")

    ax[1].bar(range(len(s)), s["total_sig_families"], color="#882255")
    ax[1].set_xticks(range(len(s))); ax[1].set_xticklabels(s["contrast"], rotation=30, ha="right")
    ax[1].set_ylabel("significant families (sum)")
    ax[1].set_title("B  Significant families", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
