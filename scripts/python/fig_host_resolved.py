#!/usr/bin/env python3
"""
fig_host_resolved.py
--------------------
Figure: the animal niche resolved by host.
  A  per-host species richness rarefied to a common genome effort (top hosts),
     so mouse's dominance does not masquerade as higher diversity.
  B  host x phylum composition (top hosts), species-weighted.
The figure makes explicit how much of the 'animal' niche is mouse.
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
    ap.add_argument("--host-summary", required=True)
    ap.add_argument("--host-phylum", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)

    hs = pd.read_csv(args.host_summary, sep="\t")
    hp = pd.read_csv(args.host_phylum, sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(mm(185), mm(80)),
                           gridspec_kw={"width_ratios": [1, 1.2]})

    # A rarefied richness (hosts that qualified)
    q = hs.dropna(subset=["rarefied_richness_mean"]).sort_values(
        "rarefied_richness_mean", ascending=True)
    ax[0].barh(range(len(q)), q["rarefied_richness_mean"],
               xerr=[np.clip(q["rarefied_richness_mean"] - q["lo"], 0, None),
                     np.clip(q["hi"] - q["rarefied_richness_mean"], 0, None)],
               color="#D55E00", error_kw=dict(lw=0.5))
    ax[0].set_yticks(range(len(q))); ax[0].set_yticklabels(q["host"], fontsize=5)
    ax[0].set_xlabel("species richness (rarefied)")
    ax[0].set_title("A  Host richness (matched effort)", loc="left", fontweight="bold")

    # B host x phylum composition
    hosts = (hp.groupby("host")["n_species"].sum().sort_values(ascending=False)
             .head(8).index.tolist())
    top_phyla = (hp.groupby("gtdb_phylum")["proportion"].max().sort_values(ascending=False)
                 .head(8).index.tolist())
    cmap = plt.get_cmap("tab20")
    cfor = {p: cmap(i % 20) for i, p in enumerate(top_phyla)}
    for xi, host in enumerate(hosts):
        sub = hp[hp["host"] == host].set_index("gtdb_phylum")["proportion"]
        bottom = 0
        for p in top_phyla:
            v = float(sub.get(p, 0))
            ax[1].bar(xi, v, 0.8, bottom=bottom, color=cfor[p], edgecolor="white", lw=0.2)
            bottom += v
    ax[1].set_xticks(range(len(hosts))); ax[1].set_xticklabels(hosts, rotation=40, ha="right", fontsize=5)
    ax[1].set_ylabel("species proportion"); ax[1].set_ylim(0, 1)
    ax[1].set_title("B  Host phylum composition", loc="left", fontweight="bold")
    handles = [plt.Rectangle((0, 0), 1, 1, color=cfor[p]) for p in top_phyla]
    ax[1].legend(handles, [p.replace("p__", "") for p in top_phyla], frameon=False,
                 fontsize=5, bbox_to_anchor=(1.01, 1), loc="upper left")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
