#!/usr/bin/env python3
"""
fig_volcano_concordance.py
--------------------------
Results figure for one layer x contrast.

  A  Volcano: effect size (consensus log2 OR, + = human-enriched) vs evidence
     (-log10 phyloglm q). Points coloured by tier; consensus signatures filled.
  B  Method concordance: phyloglm log2 OR vs CMH log2 OR. Because these are two
     independent phylogenetic controls (covariance model vs clade matching),
     points on the diagonal are robust calls. The correlation in B is the
     figure's defensive argument: the signal is not an artefact of one method.

Top consensus features are labelled in both panels.
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
    for c in ["pg_log2or", "cmh_log2or", "consensus_log2or", "pg_q"]:
        df[c] = pd.to_numeric(df.get(c), errors="coerce")
    df["consensus_signature"] = df["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])

    fig, ax = plt.subplots(1, 2, figsize=(mm(180), mm(80)))

    # A volcano
    q = df["pg_q"].clip(lower=1e-300)
    y = -np.log10(q)
    ns = ~df["consensus_signature"]
    ax[0].scatter(df.loc[ns, "consensus_log2or"], y[ns], s=3, color="#bbbbbb",
                  alpha=0.5, rasterized=True, label="not consensus")
    cs = df["consensus_signature"]
    enr = cs & (df["consensus_log2or"] > 0)
    dep = cs & (df["consensus_log2or"] < 0)
    ax[0].scatter(df.loc[enr, "consensus_log2or"], y[enr], s=7,
                  color=cfg["figures"]["palette"]["human"], label="human-enriched")
    ax[0].scatter(df.loc[dep, "consensus_log2or"], y[dep], s=7,
                  color="#999999", label="human-depleted")
    eff = cfg["stats"]["min_abs_log2_or"]
    ax[0].axvline(eff, ls="--", lw=0.5, color="#cc3311")
    ax[0].axvline(-eff, ls="--", lw=0.5, color="#cc3311")
    ax[0].axhline(-np.log10(cfg["stats"]["fdr_alpha"]), ls="--", lw=0.5, color="#cc3311")
    ax[0].set_xlabel("consensus log2 OR (human vs comparator)")
    ax[0].set_ylabel("-log10 phyloglm q")
    ax[0].set_title(f"A  Volcano — {args.layer} / {args.contrast}",
                    loc="left", fontweight="bold")
    ax[0].legend(frameon=False, fontsize=5)

    top = df[cs].reindex(df[cs]["consensus_log2or"].abs().sort_values(ascending=False).index).head(8)
    for _, r in top.iterrows():
        ax[0].annotate(str(r["feature"]), (r["consensus_log2or"], -np.log10(max(r["pg_q"], 1e-300))),
                       fontsize=5, xytext=(3, 2), textcoords="offset points")

    # B concordance
    sub = df.dropna(subset=["pg_log2or", "cmh_log2or"])
    ax[1].axhline(0, lw=0.4, color="#cccccc"); ax[1].axvline(0, lw=0.4, color="#cccccc")
    ax[1].scatter(sub.loc[~sub["consensus_signature"], "pg_log2or"],
                  sub.loc[~sub["consensus_signature"], "cmh_log2or"],
                  s=3, color="#bbbbbb", alpha=0.5, rasterized=True)
    ax[1].scatter(sub.loc[sub["consensus_signature"], "pg_log2or"],
                  sub.loc[sub["consensus_signature"], "cmh_log2or"],
                  s=7, color=cfg["figures"]["palette"]["human"])
    lim = np.nanpercentile(np.abs(sub[["pg_log2or", "cmh_log2or"]].values), 99)
    ax[1].plot([-lim, lim], [-lim, lim], ls=":", lw=0.6, color="#333333")
    if len(sub) > 3:
        r = np.corrcoef(sub["pg_log2or"], sub["cmh_log2or"])[0, 1]
        ax[1].text(0.05, 0.92, f"Pearson r = {r:.2f}", transform=ax[1].transAxes, fontsize=6)
    ax[1].set_xlabel("phyloglm log2 OR (phylogenetic covariance)")
    ax[1].set_ylabel("CMH log2 OR (clade-matched)")
    ax[1].set_title("B  Two phylogenetic controls agree", loc="left", fontweight="bold")
    fig.tight_layout()
    save(fig, args.out, cfg)


if __name__ == "__main__":
    main()
