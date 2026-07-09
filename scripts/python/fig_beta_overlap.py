#!/usr/bin/env python3
"""
fig_beta_overlap.py
-------------------
Figure: between-niche structure.
  A  Sorensen dissimilarity per niche pair split into turnover vs nestedness
     (stacked). High turnover = niches hold different lineages, not subsets.
  B  Cross-niche species overlap: observed vs null, as the standardised effect
     size (SES) per niche pair / triple. Negative SES = less overlap than chance
     (niche specificity).
"""
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--beta-pairwise", required=True)
    ap.add_argument("--overlap-null", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)

    fig, ax = plt.subplots(1, 2, figsize=(mm(170), mm(72)))

    # A turnover vs nestedness
    b = pd.read_csv(args.beta_pairwise, sep="\t")
    b["pair"] = b["niche_a"] + "\n" + b["niche_b"]
    x = np.arange(len(b))
    ax[0].bar(x, b["turnover_sim"], label="turnover", color="#4477AA")
    ax[0].bar(x, b["nestedness_sne"], bottom=b["turnover_sim"], label="nestedness",
              color="#DDAA33")
    ax[0].set_xticks(x); ax[0].set_xticklabels(b["pair"], fontsize=5)
    ax[0].set_ylabel("Sorensen dissimilarity")
    ax[0].set_title("A  Turnover vs nestedness", loc="left", fontweight="bold")
    ax[0].legend(frameon=False)

    # B overlap SES
    o = pd.read_csv(args.overlap_null, sep="\t")
    colors = ["#CC3311" if s < 0 else "#117733" for s in o["SES"]]
    ax[1].barh(range(len(o))[::-1], o["SES"], color=colors)
    ax[1].axvline(0, lw=0.5, color="#333333")
    ax[1].set_yticks(range(len(o))[::-1])
    ax[1].set_yticklabels(o["set"], fontsize=5)
    ax[1].set_xlabel("overlap SES (obs vs null)")
    ax[1].set_title("B  Overlap vs chance", loc="left", fontweight="bold")
    for i, (_, r) in enumerate(o.iterrows()):
        if r["p_two_sided"] < 0.05:
            ax[1].text(r["SES"], len(o) - 1 - i, " *", va="center", fontsize=7)
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
