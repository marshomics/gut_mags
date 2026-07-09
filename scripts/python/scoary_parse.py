#!/usr/bin/env python3
"""
scoary_parse.py
---------------
Parse a Scoary2 run for one layer x contrast into the uniform signature schema
used by the consensus. Scoary2 writes traits/<contrast>/results.tsv with, per
gene: odds_ratio, fisher_p, fisher_q, empirical_p, and the pairwise-comparison
columns contrasting / supporting / opposing / best / worst (Roder et al. 2024).

A Scoary2 "hit" requires all of: corrected Fisher q below threshold, the
label-switching permutation p below threshold (so the signal is not merely
lineage-specific), enough supporting pairwise comparisons (the structure-aware
evidence), and the 'best' pairwise p below threshold. Direction comes from the
odds ratio (>1 => enriched in the contrast's positive side).

A criterion is applied only when Scoary2 actually reported it. empirical_p is
written only when Scoary2 was run with permutations, and the pairwise columns
only when a tree was supplied; treating a column that was never produced as a
failed criterion would silently make every gene non-significant and so mute the
fourth consensus method entirely. Which criteria were applied is written to
scoary_<layer>_<contrast>.criteria.txt and logged, and if the two
structure-aware criteria (empirical_p, supporting pairs) are both missing the
run is flagged, because a plain Fisher test is not a pan-GWAS.

Output: scoary_<layer>_<contrast>.tsv with columns the consensus reads:
  feature, scoary_log2or, scoary_fisher_q, scoary_empirical_p,
  scoary_supporting, scoary_best_p, scoary_sig, scoary_dir
"""
import argparse
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("scoary-parse")
LN2 = np.log(2.0)


def col(df, *cands):
    low = {c.lower(): c for c in df.columns}
    for c in cands:
        if c.lower() in low:
            return low[c.lower()]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--results", required=True, help="traits/<contrast>/results.tsv")
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    sc = cfg["scoary"]
    if not os.path.exists(args.results) or os.path.getsize(args.results) == 0:
        pd.DataFrame(columns=["feature", "scoary_log2or", "scoary_fisher_q",
                              "scoary_empirical_p", "scoary_supporting",
                              "scoary_best_p", "scoary_sig", "scoary_dir"]).to_csv(
            args.out, sep="\t", index=False)
        log.warning("No Scoary2 results at %s; wrote empty.", args.results)
        return

    d = pd.read_csv(args.results, sep="\t")
    gene = col(d, "Gene", "gene")
    orr = col(d, "odds_ratio", "Odds_ratio")
    fq = col(d, "fisher_q", "Benjamini_H_p", "qval")
    ep = col(d, "empirical_p", "Empirical_p")
    sup = col(d, "supporting", "Max_supporting_pairs")
    best = col(d, "best", "Best_pairwise_comp_p")

    out = pd.DataFrame({"feature": d[gene]})
    o = pd.to_numeric(d[orr], errors="coerce") if orr else np.nan
    out["scoary_log2or"] = np.log2(o.clip(lower=1e-6)).clip(-10, 10) if orr else np.nan
    out["scoary_fisher_q"] = pd.to_numeric(d[fq], errors="coerce") if fq else np.nan
    out["scoary_empirical_p"] = pd.to_numeric(d[ep], errors="coerce") if ep else np.nan
    out["scoary_supporting"] = pd.to_numeric(d[sup], errors="coerce") if sup else np.nan
    out["scoary_best_p"] = pd.to_numeric(d[best], errors="coerce") if best else np.nan

    # Apply a criterion only if Scoary2 reported the column it needs. A column
    # Scoary2 never wrote is a criterion that was never evaluated, not one that
    # was failed; NaN WITHIN a reported column is still a failure, as Scoary2
    # emits a value for every gene it tested.
    if fq is None:
        log.error("Scoary2 results at %s have no Fisher q column; nothing can be "
                  "called. Writing all-NA.", args.results)
    crit, applied = [], []
    if fq is not None:
        crit.append(out["scoary_fisher_q"] < sc["fisher_q_max"]); applied.append("fisher_q")
    if ep is not None:
        crit.append(out["scoary_empirical_p"] < sc["empirical_p_max"]); applied.append("empirical_p")
    if sup is not None:
        crit.append(out["scoary_supporting"] >= sc["min_supporting_pairs"]); applied.append("supporting_pairs")
    if best is not None:
        crit.append(out["scoary_best_p"] < sc["best_pairwise_p_max"]); applied.append("best_pairwise_p")

    if crit:
        sig = crit[0].fillna(False)
        for c in crit[1:]:
            sig &= c.fillna(False)
        out["scoary_sig"] = sig
    else:
        out["scoary_sig"] = False
    out["scoary_dir"] = np.sign(out["scoary_log2or"]).replace(0, np.nan)
    out.sort_values("scoary_fisher_q").to_csv(args.out, sep="\t", index=False)

    with open(args.out.replace(".tsv", "") + ".criteria.txt", "w") as fh:
        fh.write("\n".join(applied) + "\n")
    if not {"empirical_p", "supporting_pairs"} & set(applied):
        log.warning("Scoary2 at %s reported neither empirical_p nor supporting pairs: "
                    "the population-structure correction was NOT applied and these "
                    "calls are ordinary Fisher tests. Re-run Scoary2 with "
                    "--n-permut and a tree, or drop 'scoary' from "
                    "stats.consensus.require_methods.", args.results)
    log.info("%s/%s: %d genes, %d Scoary2 hits; criteria applied: %s",
             args.layer, args.contrast, len(out), int(out["scoary_sig"].sum()),
             ", ".join(applied) or "none")


if __name__ == "__main__":
    main()
