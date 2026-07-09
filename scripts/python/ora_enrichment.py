#!/usr/bin/env python3
"""
ora_enrichment.py
-----------------
Hypergeometric over-representation analysis of the niche-signature features
within each functional category system. Two defensible choices are enforced:

  * Background = the features actually TESTED for this layer x contrast that are
    annotatable in the system (have >=1 category). Using the tested, annotatable
    universe (not the whole database) is what keeps the enrichment honest.
  * Direction is handled separately: 'up' = features enriched in the contrast's
    positive niche, 'down' = depleted, plus 'all'. Mixing directions would let
    opposing signals cancel.

For each category: k = signature features in it, n = signature features (in
background), K = background features in it, N = background size; p from the
hypergeometric survival function, fold = (k/n)/(K/N), BH-corrected within each
system x direction. Overlapping (leading) features are listed.

BH is applied over EVERY category that passes the size filter, including those
with no overlap (p = 1), because those categories were hypotheses too. Dropping
them before correction, as is common, shrinks the denominator and inflates
significance; only the k >= 1 rows are written out.

Output: enrichment_ora_<layer>_<contrast>.tsv
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from hgn_utils import load_config, get_logger

log = get_logger("ora")


def positive_label(contrast):
    return "host" if contrast == "host_vs_free" else contrast.split("_vs_")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--genesets", required=True)
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    en = cfg["enrichment"]
    alpha = en["fdr_alpha"]
    pos = positive_label(args.contrast)

    sig = pd.read_csv(args.signatures, sep="\t")
    sig["feature"] = sig["feature"].astype(str)
    sig["consensus_signature"] = sig["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    tested = set(sig["feature"])
    up = set(sig.loc[sig["consensus_signature"] & (sig["direction"] == f"{pos}_enriched"), "feature"])
    down = set(sig.loc[sig["consensus_signature"] & (sig["direction"] == f"{pos}_depleted"), "feature"])
    allsig = up | down
    fg_sets = {"up": up, "down": down, "all": allsig}

    gs = pd.read_csv(args.genesets, sep="\t")
    gs["feature"] = gs["feature"].astype(str)
    name_of = dict(zip(gs["category_id"], gs["category_name"]))

    rows = []
    for system, sysdf in gs.groupby("system"):
        members = sysdf.groupby("category_id")["feature"].apply(set).to_dict()
        annotatable = set().union(*members.values()) & tested
        N = len(annotatable)
        if N < en["min_set_size"]:
            continue
        for direction, fg in fg_sets.items():
            fg_ann = fg & annotatable
            n = len(fg_ann)
            if n == 0:
                continue
            recs = []
            for cid, feats in members.items():
                K = len(feats & annotatable)
                if K < en["min_set_size"] or K > en["max_set_size"]:
                    continue
                k = len(fg_ann & feats)
                # k == 0 categories are still hypotheses: keep them for BH (p=1),
                # drop them only from the written output.
                p = float(hypergeom.sf(k - 1, N, K, n)) if k else 1.0
                fold = (k / n) / (K / N) if K and n else np.nan
                recs.append({"system": system, "category_id": cid,
                             "category_name": name_of.get(cid, cid),
                             "direction": direction, "k": k, "n": n, "K": K, "N": N,
                             "fold_enrichment": fold, "p": p,
                             "overlap": ",".join(sorted(fg_ann & feats)[:25])})
            if recs:
                rd = pd.DataFrame(recs)
                rd["q"] = multipletests(rd["p"], method="fdr_bh")[1]
                rd["n_categories_tested"] = len(rd)
                rows.append(rd[rd["k"] > 0])
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["system", "category_id", "category_name", "direction", "k", "n",
                 "K", "N", "fold_enrichment", "p", "q", "overlap", "n_categories_tested"])
    out.sort_values(["system", "direction", "q"]).to_csv(args.out, sep="\t", index=False)
    log.info("%s/%s ORA: %d category tests, %d at q<%.2f",
             args.layer, args.contrast, len(out), int((out["q"] < alpha).sum()) if len(out) else 0, alpha)


if __name__ == "__main__":
    main()
