#!/usr/bin/env python3
"""
transition_verdict.py
---------------------
Integrate the per-species evidence into a recent-acquisition call for each
non-ancestral niche. No single statistic is trusted; a call is only as strong as
the number of independent lines that agree:

  gate   niche membership is phylogenetically structured (Slatkin-Maddison p<0.05)
  e1     transitions into the niche are shallow (near the tips)         -> recent
  e2     lower within-niche diversity than the source (founder effect)
  e3     Tajima's D negative in the niche (post-founder expansion)
  e4     the niche's variation is nested in the source (directionality)
  e5     accessory genes were gained in the niche (colonisation cargo)

The ancestral niche is the tree root state; derived niches are the others. Tier:
strong (gate + nestedness + >=3 of e1-e3,e5), moderate (gate + >=2), else weak.
The script is defensive to missing optional inputs.

Output: transition_verdict.tsv (one row per species x derived niche)
"""
import argparse
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("verdict")


def load(path):
    return pd.read_csv(path, sep="\t") if os.path.exists(path) and os.path.getsize(path) else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-id", required=True)
    ap.add_argument("--ancestral", required=True)     # ancestral_summary.tsv
    ap.add_argument("--popgen", required=True)        # popgen_diversity.tsv
    ap.add_argument("--divergence", required=True)
    ap.add_argument("--directionality", required=True)
    ap.add_argument("--accessory", required=True)     # accessory_summary.tsv
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    anc = load(args.ancestral)
    pg = load(args.popgen)
    dirn = load(args.directionality)
    acc = load(args.accessory)
    if anc.empty:
        pd.DataFrame([{"species_id": args.species_id, "status": "no_ancestral"}]).to_csv(
            args.out, sep="\t", index=False); return
    a = anc.iloc[0]
    root = a.get("root_state")
    sm_p = float(a.get("sm_p", np.nan))
    niches = str(a.get("niches", "")).split(",")
    gate = (sm_p < cfg["stats"]["fdr_alpha"])

    pi = dict(zip(pg["niche"], pg["pi_mean"])) if "pi_mean" in pg.columns else {}
    taj = dict(zip(pg["niche"], pg["tajimaD_mean"])) if "tajimaD_mean" in pg.columns else {}

    rows = []
    for d in [n for n in niches if n and n != root]:
        depth = a.get(f"mean_depth_into_{d}", np.nan)
        src = root
        # directionality derived call for the (src,d) pair
        nested_derived = None
        if not dirn.empty:
            for _, r in dirn.iterrows():
                if {r.get("X"), r.get("Y")} == {src, d}:
                    nested_derived = r.get("derived_call")
        # accessory gains into d
        acc_gain = np.nan
        if not acc.empty:
            col = f"enriched_in_{d}"
            for _, r in acc.iterrows():
                if col in acc.columns and {r["pair"].split("|")[0], r["pair"].split("|")[1]} == {src, d}:
                    acc_gain = r.get(col, np.nan)
        e1 = (depth < 0.5) if pd.notna(depth) else False
        e2 = (pi.get(d, np.inf) < pi.get(src, -np.inf)) if pi else False
        e3 = (taj.get(d, np.inf) < 0) if taj else False
        e4 = (nested_derived == d) if nested_derived else False
        e5 = (acc_gain and acc_gain > 0) if pd.notna(acc_gain) else False
        n_support = int(e1) + int(e2) + int(e3) + int(e5)
        if gate and e4 and n_support >= 3:
            tier = "strong"
        elif gate and (n_support + int(e4)) >= 2:
            tier = "moderate"
        else:
            tier = "weak"
        rows.append({
            "species_id": args.species_id, "ancestral_niche": root,
            "derived_niche": d, "structured_gate_smp": round(sm_p, 4),
            "transition_depth": round(depth, 3) if pd.notna(depth) else np.nan,
            "pi_derived": pi.get(d, np.nan), "pi_source": pi.get(src, np.nan),
            "tajimaD_derived": taj.get(d, np.nan),
            "nested_derived_call": nested_derived,
            "accessory_genes_gained": acc_gain,
            "e1_shallow": e1, "e2_lower_div": e2, "e3_neg_tajD": e3,
            "e4_nested": e4, "e5_gene_gain": e5,
            "n_support": n_support + int(e4), "tier": tier})
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    log.info("[%s] root=%s -> %s", args.species_id, root,
             [(r["derived_niche"], r["tier"]) for r in rows])


if __name__ == "__main__":
    main()
