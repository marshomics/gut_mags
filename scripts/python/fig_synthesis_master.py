#!/usr/bin/env python3
"""
fig_synthesis_master.py
-----------------------
The master figure tying the manuscript together.
  A  counts of human-specific species, genes and functions.
  B  adaptation-mode breakdown of the human signatures (acquisition / loss / ...).
  C  human-gut selective pressures that are enriched (fold enrichment).
  D  diversification rate (DR) by niche.
PNG + editable-text SVG.
"""
import argparse
import json

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--adaptation-summary", required=True)
    ap.add_argument("--ecological", required=True)
    ap.add_argument("--diversification", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    try:
        counts = json.load(open(args.counts))
    except Exception:
        counts = {}

    fig, ax = plt.subplots(2, 2, figsize=(mm(180), mm(150)))

    # A counts
    keys = [("n_human_specific_species", "species"),
            ("n_human_specific_genes", "genes/features"),
            ("n_human_specific_functions", "functions")]
    vals = [counts.get(k, 0) for k, _ in keys]
    ax[0, 0].bar(range(3), vals, color=["#0072B2", "#56B4E9", "#009E73"])
    ax[0, 0].set_xticks(range(3)); ax[0, 0].set_xticklabels([l for _, l in keys])
    ax[0, 0].set_ylabel("count")
    for i, v in enumerate(vals):
        ax[0, 0].text(i, v, str(v), ha="center", va="bottom", fontsize=6)
    ax[0, 0].set_title("A  Human-specific catalogue", loc="left", fontweight="bold")

    # B adaptation mode
    ad = load(args.adaptation_summary)
    if len(ad) and "mode" in ad.columns:
        mc = ad.groupby("mode")["n"].sum().sort_values(ascending=False)
        ax[0, 1].bar(range(len(mc)), mc.values, color="#CC79A7")
        ax[0, 1].set_xticks(range(len(mc))); ax[0, 1].set_xticklabels(mc.index, rotation=25, ha="right", fontsize=5)
        ax[0, 1].set_ylabel("signature features")
    else:
        ax[0, 1].text(0.5, 0.5, "adaptation modes n/a", ha="center")
    ax[0, 1].set_title("B  Adaptation mode", loc="left", fontweight="bold")

    # C ecological pressures
    eco = load(args.ecological)
    if len(eco) and "fold_enrichment" in eco.columns:
        e = eco.sort_values("fold_enrichment", ascending=True)
        # 'enriched' is absent when no pressure reached significance; e.get(...)
        # would return the scalar False and iterating it raises.
        enr = (e["enriched"].astype(bool) if "enriched" in e.columns
               else pd.Series(False, index=e.index))
        colors = ["#D55E00" if s else "#bbb" for s in enr]
        ax[1, 0].barh(range(len(e)), e["fold_enrichment"], color=colors)
        ax[1, 0].set_yticks(range(len(e))); ax[1, 0].set_yticklabels(e["pressure"], fontsize=5)
        ax[1, 0].set_xlabel("fold enrichment (human signatures)")
        ax[1, 0].axvline(1, ls="--", lw=0.5, color="#333")
    else:
        ax[1, 0].text(0.5, 0.5, "pressures n/a", ha="center")
    ax[1, 0].set_title("C  Human-gut selective pressures", loc="left", fontweight="bold")

    # D diversification
    dv = load(args.diversification)
    if len(dv) and "DR" in dv.columns:
        niches = [n for n in cfg["inputs"]["niche_levels"] if n in set(dv["niche"])]
        data = [dv.loc[dv["niche"] == n, "DR"].dropna().values for n in niches]
        parts = ax[1, 1].violinplot(data, showmedians=True)
        for i, b in enumerate(parts["bodies"]):
            b.set_facecolor(pal.get(niches[i], "#777")); b.set_alpha(0.6)
        ax[1, 1].set_xticks(range(1, len(niches) + 1)); ax[1, 1].set_xticklabels(niches)
        ax[1, 1].set_ylabel("DR (tip diversification rate)")
    else:
        ax[1, 1].text(0.5, 0.5, "diversification n/a", ha="center")
    ax[1, 1].set_title("D  Diversification by niche", loc="left", fontweight="bold")

    fig.suptitle("What makes the human gut human", x=0.02, ha="left", fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
