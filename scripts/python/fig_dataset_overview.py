#!/usr/bin/env python3
"""
fig_dataset_overview.py
-----------------------
Figure 1: the dataset and its built-in imbalances, shown openly. Four panels,
each one a confounder the analysis has to control:

  A  genomes-per-species rank curve (log-log)  -> strains-per-species skew
  B  genomes vs species per niche               -> species-per-niche imbalance
  C  genome quality (completeness vs contam.)   -> quality covariate
  D  animal-niche host composition              -> mouse dominance

Putting these on the first figure makes the rest of the pipeline legible: every
later control maps to a panel here.
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
    ap.add_argument("--samples", required=True)
    ap.add_argument("--strains-per-species", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]

    samples = pd.read_parquet(args.samples)
    spc = pd.read_csv(args.strains_per_species, sep="\t")

    fig, ax = plt.subplots(2, 2, figsize=(mm(180), mm(150)))

    # A: rank curve
    s = spc.sort_values("n_genomes", ascending=False).reset_index(drop=True)
    ax[0, 0].plot(np.arange(1, len(s) + 1), s["n_genomes"], lw=0.8, color="#333333")
    ax[0, 0].set_xscale("log"); ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlabel("species rank"); ax[0, 0].set_ylabel("genomes per species")
    ax[0, 0].set_title("A  Strains-per-species skew", loc="left", fontweight="bold")
    med = s["n_genomes"].median()
    ax[0, 0].axhline(med, ls="--", lw=0.6, color="#888888")
    ax[0, 0].text(1.2, med * 1.3, f"median {med:.0f}", fontsize=6, color="#555555")
    for i in range(min(3, len(s))):
        ax[0, 0].annotate(s.loc[i, "species"].replace("s__", ""),
                          (i + 1, s.loc[i, "n_genomes"]), fontsize=5,
                          xytext=(6, 0), textcoords="offset points", va="center")

    # B: genomes vs species per niche
    gper = samples["niche"].value_counts().reindex(niches)
    sper = samples.groupby("niche")["species"].nunique().reindex(niches)
    x = np.arange(len(niches)); w = 0.38
    ax[0, 1].bar(x - w / 2, gper.values, w, label="genomes",
                 color=[pal[n] for n in niches], alpha=0.95)
    ax[0, 1].bar(x + w / 2, sper.values, w, label="species",
                 color=[pal[n] for n in niches], alpha=0.5, hatch="//")
    ax[0, 1].set_yscale("log"); ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(niches)
    ax[0, 1].set_ylabel("count (log)")
    ax[0, 1].set_title("B  Species-per-niche imbalance", loc="left", fontweight="bold")
    ax[0, 1].legend(frameon=False, loc="upper right")

    # C: quality
    for n in niches:
        d = samples[samples["niche"] == n]
        ax[1, 0].scatter(d["completeness"], d["contamination"], s=1, alpha=0.05,
                         color=pal[n], rasterized=True, label=n)
    ax[1, 0].set_xlabel("completeness (%)"); ax[1, 0].set_ylabel("contamination (%)")
    ax[1, 0].set_title("C  Genome quality", loc="left", fontweight="bold")
    lg = ax[1, 0].legend(frameon=False, markerscale=6, loc="upper left")
    for h in lg.legend_handles:
        h.set_alpha(1)

    # D: host composition (animal)
    host_col = cfg["inputs"]["columns"]["host_common_name"]
    if host_col in samples.columns:
        hosts = (samples.loc[samples["niche"] == "animal", host_col]
                 .replace("", np.nan).dropna().value_counts().head(10))
        ax[1, 1].barh(hosts.index[::-1], hosts.values[::-1],
                      color=pal["animal"])
        ax[1, 1].set_xlabel("genomes"); ax[1, 1].set_xscale("log")
        ax[1, 1].set_title("D  Animal host dominance", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
