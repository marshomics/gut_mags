#!/usr/bin/env python3
"""
fig_redundancy.py
-----------------
Functional redundancy of the human gut.
  A  occupancy distribution: how many species carry each function (log), with the
     core / intermediate / rare split.
  B  accumulation: unique functions vs species sampled; the dashed line is y=x
     (species). Functions saturating below it is redundancy.
  C  Ricotta relative functional redundancy (1 - Q/D) per layer (and population).
  D  carrier spread: number of families carrying each function (robust vs single-
     clade).
Primary layer for A/B/D is the first redundancy layer (KO). PNG + editable SVG.
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def load(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--redundancy-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    L = cfg["redundancy"]["layers"][0]
    d = args.redundancy_dir

    occ = load(f"{d}/occupancy_{L}.tsv")
    acc = load(f"{d}/accumulation_{L}.tsv")
    summ = load(f"{d}/redundancy_summary.tsv")
    spread = load(f"{d}/spread_{L}.tsv")

    fig, ax = plt.subplots(2, 2, figsize=(mm(180), mm(150)))

    # A occupancy
    if len(occ):
        ax[0, 0].hist(occ["n_species"], bins=40, color="#0072B2")
        ax[0, 0].set_yscale("log")
        ax[0, 0].set_xlabel(f"species carrying a {L} function")
        ax[0, 0].set_ylabel("functions (log)")
        if "class" in occ.columns:
            cc = occ["class"].value_counts()
            txt = " | ".join(f"{k}:{v}" for k, v in cc.items())
            ax[0, 0].text(0.5, 0.95, txt, transform=ax[0, 0].transAxes, fontsize=5, va="top")
    ax[0, 0].set_title("A  Function occupancy", loc="left", fontweight="bold")

    # B accumulation
    if len(acc):
        ax[0, 1].plot(acc["k"], acc["functions_mean"], color="#009E73", lw=1.2, label="functions")
        ax[0, 1].fill_between(acc["k"], acc["lo"], acc["hi"], color="#009E73", alpha=0.2, lw=0)
        mx = acc["k"].max()
        ax[0, 1].plot([1, mx], [1, mx], ls="--", lw=0.6, color="#888", label="species (y=x)")
        ax[0, 1].set_xlabel("species sampled"); ax[0, 1].set_ylabel(f"unique {L} functions")
        ax[0, 1].legend(frameon=False)
    ax[0, 1].set_title("B  Taxonomic vs functional accumulation", loc="left", fontweight="bold")

    # C Ricotta relFR
    if len(summ) and "ricotta_relFR" in summ.columns:
        s = summ[summ["community"] == "human"] if "community" in summ.columns else summ
        ax[1, 0].bar(range(len(s)), s["ricotta_relFR"], color="#CC79A7")
        ax[1, 0].set_xticks(range(len(s))); ax[1, 0].set_xticklabels(s["layer"], rotation=20, ha="right")
        ax[1, 0].set_ylabel("relative FR (1 - Q/D)"); ax[1, 0].set_ylim(0, 1)
    ax[1, 0].set_title("C  Ricotta functional redundancy", loc="left", fontweight="bold")

    # D carrier spread
    if len(spread) and "n_families" in spread.columns:
        ax[1, 1].hist(spread["n_families"], bins=40, color="#D55E00")
        ax[1, 1].set_yscale("log")
        ax[1, 1].set_xlabel("families carrying a function"); ax[1, 1].set_ylabel("functions (log)")
    ax[1, 1].set_title("D  Carrier phylogenetic spread", loc="left", fontweight="bold")

    fig.suptitle("Functional redundancy of the human gut", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
