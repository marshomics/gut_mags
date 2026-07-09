#!/usr/bin/env python3
"""
ecological_pressure_test.py
---------------------------
Tests the "why the human gut" hypotheses explicitly. Curated gene sets encoding
known human-gut selective pressures (host-glycan/mucin foraging, bile tolerance,
oxidative-stress handling, SCFA fermentation, dietary-fiber and starch
degradation) are tested for over-representation among the human-enriched
signature features, across the KO and CAZyme layers together.

Background is the tested features (from the combined signatures across layers)
that belong to any curated pressure set, so enrichment is judged against what
was actually examined, not the whole database. Hypergeometric test, BH-corrected
across pressures. The pressures that come out enriched are the candidate reasons
these species occupy the human gut; the curated map is editable.

Output: ecological_pressure_<contrast>.tsv (pressure, fold, p, q, overlap genes)
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from hgn_utils import load_config, get_logger

log = get_logger("eco-pressure")


def positive_label(contrast):
    return "host" if contrast == "host_vs_free" else contrast.split("_vs_")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--signatures", required=True, help="combined signatures_all for the contrast")
    ap.add_argument("--pressure-map", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    alpha = cfg["enrichment"]["fdr_alpha"]
    pos = positive_label(args.contrast)

    pm = pd.read_csv(args.pressure_map, sep="\t")
    pm["feature"] = pm["feature"].astype(str)
    pressure_feats = pm.groupby("category_id")["feature"].apply(set).to_dict()
    pressure_name = dict(zip(pm["category_id"], pm["category_name"]))
    all_pressure_features = set(pm["feature"])

    sig = pd.read_csv(args.signatures, sep="\t")
    sig["feature"] = sig["feature"].astype(str)
    sig["consensus_signature"] = sig["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    tested = set(sig["feature"])
    annotatable = tested & all_pressure_features
    N = len(annotatable)
    fg = set(sig.loc[sig["consensus_signature"] & (sig["direction"] == f"{pos}_enriched"),
                     "feature"]) & annotatable
    n = len(fg)

    rows = []
    for pid, feats in pressure_feats.items():
        K = len(feats & annotatable)
        if K == 0:
            continue
        k = len(fg & feats)
        p = float(hypergeom.sf(k - 1, N, K, n)) if (n and k) else 1.0
        fold = (k / n) / (K / N) if (n and K and N) else np.nan
        rows.append({"pressure": pid, "name": pressure_name.get(pid, pid),
                     "k_signature": k, "n_signature": n, "K_background": K,
                     "N_background": N, "fold_enrichment": fold, "p": p,
                     "overlap": ",".join(sorted(fg & feats))})
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = multipletests(out["p"], method="fdr_bh")[1]
        out["enriched"] = out["q"] < alpha
    out.sort_values("p").to_csv(args.out, sep="\t", index=False)
    log.info("%s ecological pressures: %d tested, %d enriched (q<%.2f); N_bg=%d, fg=%d",
             args.contrast, len(out), int(out["enriched"].sum()) if len(out) else 0,
             alpha, N, n)


if __name__ == "__main__":
    main()
