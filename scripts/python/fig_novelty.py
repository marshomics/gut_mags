#!/usr/bin/env python3
"""
fig_novelty.py
--------------
Figure: undescribed diversity per niche.
  A  novel (placeholder) species fraction per niche: observed vs rarefied to a
     common sampling effort (with CI). The rarefied bars are the fair comparison.
  B  placeholder-taxa fraction across ranks per niche (where novelty sits in the
     hierarchy).
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
    ap.add_argument("--by-niche", required=True)
    ap.add_argument("--by-rank", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    pal = niche_palette(cfg)
    niches = cfg["inputs"]["niche_levels"]

    bn = pd.read_csv(args.by_niche, sep="\t").set_index("niche").reindex(niches)
    br = pd.read_csv(args.by_rank, sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(mm(170), mm(70)))
    x = np.arange(len(niches)); w = 0.38
    ax[0].bar(x - w / 2, bn["frac_novel"], w, color=[pal[n] for n in niches],
              alpha=0.5, label="observed")
    err = [np.clip(bn["frac_novel_rarefied_mean"] - bn["lo"], 0, None),
           np.clip(bn["hi"] - bn["frac_novel_rarefied_mean"], 0, None)]
    ax[0].bar(x + w / 2, bn["frac_novel_rarefied_mean"], w,
              color=[pal[n] for n in niches], yerr=err, error_kw=dict(lw=0.6),
              label="rarefied")
    ax[0].set_xticks(x); ax[0].set_xticklabels(niches)
    ax[0].set_ylabel("novel species fraction")
    ax[0].set_title("A  Undescribed species", loc="left", fontweight="bold")
    ax[0].legend(frameon=False)

    ranks = list(dict.fromkeys(br["rank"]))
    for n in niches:
        vals = [br[(br["rank"] == r) & (br["niche"] == n)]["frac_placeholder"].mean()
                for r in ranks]
        ax[1].plot(range(len(ranks)), vals, "-o", color=pal[n], ms=3, label=n)
    ax[1].set_xticks(range(len(ranks))); ax[1].set_xticklabels(ranks, rotation=30, ha="right")
    ax[1].set_ylabel("placeholder-taxa fraction")
    ax[1].set_title("B  Novelty across ranks", loc="left", fontweight="bold")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
