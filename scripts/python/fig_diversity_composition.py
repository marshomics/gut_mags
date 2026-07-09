#!/usr/bin/env python3
"""
fig_diversity_composition.py
----------------------------
Figure 2: niche taxonomic structure, corrected for sampling.

  A  rarefaction curves (species richness vs matched genome effort) with 95% CI
  B  Hill numbers q=0,1,2 per niche (richness / exp-Shannon / inv-Simpson) + CI
  C  phylum composition, species-weighted vs genome-weighted side by side, so
     the strain-sampling bias in the naive (genome-weighted) view is visible.

The contrast in C is the point: "free-living looks the most diverse" is partly a
sampling artefact, which is why rarefied/Hill comparisons (A,B) are used.
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
    ap.add_argument("--rarefaction", required=True)
    ap.add_argument("--hill", required=True)
    ap.add_argument("--composition-phylum", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]

    rare = pd.read_csv(args.rarefaction, sep="\t")
    hill = pd.read_csv(args.hill, sep="\t")
    comp = pd.read_csv(args.composition_phylum, sep="\t")

    fig = plt.figure(figsize=(mm(180), mm(75)))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.3], wspace=0.45)

    # A rarefaction
    axA = fig.add_subplot(gs[0, 0])
    for n in niches:
        d = rare[rare["niche"] == n]
        axA.plot(d["effort"], d["richness_mean"], color=pal[n], lw=1, label=n)
        axA.fill_between(d["effort"], d["lo"], d["hi"], color=pal[n], alpha=0.2, lw=0)
    axA.set_xlabel("genomes sampled"); axA.set_ylabel("species observed")
    axA.set_title("A  Rarefaction", loc="left", fontweight="bold")
    axA.legend(frameon=False)

    # B Hill
    axB = fig.add_subplot(gs[0, 1])
    qs = [0, 1, 2]; x = np.arange(len(qs)); w = 0.25
    for i, n in enumerate(niches):
        d = hill[hill["niche"] == n].set_index("q").reindex(qs)
        axB.bar(x + (i - 1) * w, d["hill_mean"], w, color=pal[n],
                yerr=[np.clip(d["hill_mean"] - d["lo"], 0, None),
                      np.clip(d["hi"] - d["hill_mean"], 0, None)],
                error_kw=dict(lw=0.5), label=n)
    axB.set_yscale("log"); axB.set_xticks(x)
    axB.set_xticklabels([f"q={q}" for q in qs])
    axB.set_ylabel("Hill number (log)")
    axB.set_title("B  Hill diversity", loc="left", fontweight="bold")

    # C composition (top phyla), species- vs genome-weighted
    axC = fig.add_subplot(gs[0, 2])
    top = (comp.groupby("taxon")["proportion"].max().sort_values(ascending=False)
           .head(8).index.tolist())
    cmap = plt.get_cmap("tab20")
    color_for = {t: cmap(i % 20) for i, t in enumerate(top)}
    bar_x, labels = [], []
    pos = 0
    for n in niches:
        for wgt, hatch in [("species_weighted", None), ("genome_weighted", "//")]:
            sub = comp[(comp["niche"] == n) & (comp["weighting"] == wgt)]
            sub = sub[sub["taxon"].isin(top)].set_index("taxon")["proportion"]
            bottom = 0
            for t in top:
                v = float(sub.get(t, 0))
                axC.bar(pos, v, 0.8, bottom=bottom, color=color_for[t], hatch=hatch,
                        edgecolor="white", lw=0.2)
                bottom += v
            bar_x.append(pos); labels.append(f"{n}\n{'sp' if wgt.startswith('species') else 'gn'}")
            pos += 1
        pos += 0.4
    axC.set_xticks(bar_x); axC.set_xticklabels(labels, fontsize=5)
    axC.set_ylabel("proportion"); axC.set_ylim(0, 1)
    axC.set_title("C  Phylum composition (sp vs gn weighted)", loc="left", fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_for[t]) for t in top]
    axC.legend(handles, [t.replace("p__", "") for t in top], frameon=False,
               fontsize=5, ncol=1, bbox_to_anchor=(1.01, 1), loc="upper left")
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
