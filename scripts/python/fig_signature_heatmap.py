#!/usr/bin/env python3
"""
fig_signature_heatmap.py
------------------------
Figure: the human functional signature at a glance. Rows are the top consensus
signature features (across all layers), columns are niches, colour is the
species-weighted prevalence (fraction of that niche's species carrying the
feature). Row side-colours mark the functional layer (KO/CAZyme/BGC/AMR/...).
Features are hierarchically clustered so co-occurring functions group together.

Species-weighting (not genome-weighting) keeps the strains-per-species bias out
of the prevalences shown here.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list

from hgn_utils import load_config
from plotting_theme import apply_theme, save, mm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--signatures-combined", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    niches = cfg["inputs"]["niche_levels"]
    maxf = cfg["figures"]["max_heatmap_features"]

    sig = pd.read_csv(args.signatures_combined, sep="\t")
    sig["consensus_signature"] = sig["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    sig = sig[sig["consensus_signature"]].copy()
    sig["consensus_log2or"] = pd.to_numeric(sig["consensus_log2or"], errors="coerce")
    sig = sig.reindex(sig["consensus_log2or"].abs().sort_values(ascending=False).index).head(maxf)
    feat_layer = dict(zip(sig["feature"], sig.get("layer", "na")))
    feats = sig["feature"].tolist()

    sp = pd.read_csv(args.species_table, sep="\t")[["species", "niche_primary"]]
    niche_n = sp["niche_primary"].value_counts()

    # per-niche species prevalence for each signature feature
    rows = {}
    for path in glob.glob(f"{args.profiles_dir}/prevalence_*.parquet"):
        p = pd.read_parquet(path)
        p = p[p["feature"].isin(feats) & (p["present"] == 1)]
        if p.empty:
            continue
        p = p.merge(sp, on="species", how="left")
        cnt = p.groupby(["feature", "niche_primary"])["species"].nunique().unstack(fill_value=0)
        for niche in niches:
            if niche not in cnt.columns:
                cnt[niche] = 0
        prev = cnt[niches].div([niche_n.get(n, 1) for n in niches], axis=1)
        for f, r in prev.iterrows():
            rows[f] = r
    mat = pd.DataFrame(rows).T.reindex(feats).fillna(0)[niches]
    if mat.empty:
        fig, ax = plt.subplots(); ax.text(0.5, 0.5, "no signature features", ha="center")
        save(fig, args.out, cfg); return

    order = leaves_list(linkage(mat.values, method="average")) if len(mat) > 2 else range(len(mat))
    mat = mat.iloc[list(order)]

    fig, ax = plt.subplots(figsize=(mm(95), mm(min(230, 4 + 2.6 * len(mat)))))
    im = ax.imshow(mat.values, aspect="auto", cmap="rocket_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(niches))); ax.set_xticklabels(niches)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels([f"{f}" for f in mat.index], fontsize=5)
    # row side colours by layer
    layers = sorted(set(feat_layer.get(f, "na") for f in mat.index))
    lcol = {l: plt.get_cmap("tab10")(i) for i, l in enumerate(layers)}
    for i, f in enumerate(mat.index):
        ax.add_patch(plt.Rectangle((-0.7, i - 0.5), 0.25, 1,
                                   color=lcol[feat_layer.get(f, "na")], clip_on=False))
    ax.set_xlim(-0.8, len(niches) - 0.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("species prevalence")
    handles = [plt.Rectangle((0, 0), 1, 1, color=lcol[l]) for l in layers]
    ax.legend(handles, layers, frameon=False, fontsize=5,
              bbox_to_anchor=(1.25, 1), loc="upper left", title="layer")
    ax.set_title("Human functional signature (consensus features)",
                 loc="left", fontweight="bold")
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
