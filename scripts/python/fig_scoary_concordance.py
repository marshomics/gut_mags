#!/usr/bin/env python3
"""
fig_scoary_concordance.py
-------------------------
Scoary2 as a fourth method, for one layer x contrast.
  A  phyloglm log2 OR vs Scoary2 log2 OR, points coloured by whether the feature
     is a 4-method consensus signature. Agreement on the diagonal means Scoary2's
     pairwise-comparison correction reaches the same calls as the phylogenetic
     regression.
  B  number of features called by each method and by the consensus, so the added
     value (and stringency) of requiring Scoary2 too is visible.
Reads the consensus table (which now carries the Scoary2 columns).
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
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    apply_theme(cfg)
    df = pd.read_csv(args.signatures, sep="\t")
    if "scoary_log2or" not in df.columns:
        fig, ax = plt.subplots(figsize=(mm(90), mm(60)))
        ax.text(0.5, 0.5, "Scoary2 not run for this contrast", ha="center")
        save(fig, args.out, cfg); return
    df["consensus_signature"] = df["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    for c in ["pg_log2or", "scoary_log2or"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    fig, ax = plt.subplots(1, 2, figsize=(mm(170), mm(75)))
    sub = df.dropna(subset=["pg_log2or", "scoary_log2or"])
    ax[0].axhline(0, lw=0.4, color="#ccc"); ax[0].axvline(0, lw=0.4, color="#ccc")
    ax[0].scatter(sub.loc[~sub["consensus_signature"], "pg_log2or"],
                  sub.loc[~sub["consensus_signature"], "scoary_log2or"],
                  s=3, color="#bbb", alpha=0.5, rasterized=True)
    ax[0].scatter(sub.loc[sub["consensus_signature"], "pg_log2or"],
                  sub.loc[sub["consensus_signature"], "scoary_log2or"],
                  s=8, color=cfg["figures"]["palette"]["human"])
    if len(sub) > 3:
        r = np.corrcoef(sub["pg_log2or"], sub["scoary_log2or"])[0, 1]
        lim = np.nanpercentile(np.abs(sub[["pg_log2or", "scoary_log2or"]].values), 99)
        ax[0].plot([-lim, lim], [-lim, lim], ":", lw=0.6, color="#333")
        ax[0].text(0.05, 0.92, f"r = {r:.2f}", transform=ax[0].transAxes, fontsize=6)
    ax[0].set_xlabel("phyloglm log2 OR"); ax[0].set_ylabel("Scoary2 log2 OR")
    ax[0].set_title(f"A  Concordance — {args.layer}/{args.contrast}", loc="left", fontweight="bold")

    flags = {"phyloglm": "pg_sig", "CMH": "cmh_sig", "resampling": "rs_sig",
             "Scoary2": "scoary_sig", "consensus": "consensus_signature"}
    counts = [int(df[c].astype(str).isin(["True", "TRUE", "1"]).sum()) if c in df.columns else 0
              for c in flags.values()]
    ax[1].bar(range(len(flags)), counts, color="#4477AA")
    ax[1].set_xticks(range(len(flags))); ax[1].set_xticklabels(list(flags), rotation=30, ha="right")
    ax[1].set_ylabel("features called")
    ax[1].set_title("B  Calls per method", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
