#!/usr/bin/env python3
"""
accessory_differentiation.py
----------------------------
Which genes were gained or lost with the niche. From the Panaroo
gene_presence_absence matrix restricted to the species' dereplicated strains,
each accessory gene is tested for differential presence between the niche
populations (Fisher's exact test, log2 odds ratio, BH FDR). Genes enriched in
the derived niche are candidate acquisitions that accompanied colonisation; this
connects the strain-level transition story to gene content and feeds the later
gene/function stage.

Dereplication (one genome per within-niche clonal cluster) keeps a clonal
expansion from inflating a gene's apparent niche association.

Two confounders are controlled explicitly.

  Strain imbalance. Fisher's exact test on the full data is valid, but its POWER
  depends on sample size, so with unequal strain counts more genes reach
  significance in the larger population and the raw "genes gained in X vs Y"
  counts are not comparable. Since that comparison is precisely the gene-gain
  evidence line in transition_verdict.py, the counts are recomputed on BALANCED
  draws: both populations subsampled to the smaller, BH within each draw, and a
  gene counted when it is significant in at least `min_support` of the draws.
  With equal row margins the exact two-sided p is exactly twice the smaller
  hypergeometric tail (the null distribution of the cell count is symmetric), so
  the balanced test is computed in closed form and vectorised over genes.

  Gene-family multiplicity. BH runs over every accessory gene tested, including
  those with no prevalence difference; nothing is filtered before correction.

Outputs:
  accessory_differentiation.tsv  per gene: full-data Fisher stats + balanced
      support fractions (support_X, support_Y) and the balanced call
  accessory_summary.tsv          per pair: enriched_in_<niche> from the BALANCED
      analysis (used downstream), plus the unbalanced full-data counts
"""
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("accessory")


def balanced_two_sided_p(a, c, n):
    """Exact two-sided p for a 2x2 table with EQUAL row margins n.

    Conditioning on the column margin m = a + c, the count a is hypergeometric
    with population 2n, m successes and n draws. That distribution is symmetric
    about m/2 when both rows have size n, so doubling the smaller tail is exact
    rather than an approximation. Vectorised over genes.
    """
    m = a + c
    upper = hypergeom.sf(a - 1, 2 * n, m, n)     # P(A >= a)
    lower = hypergeom.cdf(a, 2 * n, m, n)        # P(A <= a)
    return np.minimum(1.0, 2.0 * np.minimum(upper, lower))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--rtab", required=True, help="Panaroo gene_presence_absence.Rtab")
    ap.add_argument("--niche-map", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fdr = cfg["stats"]["fdr_alpha"]
    acfg = cfg["transition"].get("accessory", {})
    B = int(acfg.get("balanced_iterations", 200))
    min_support = float(acfg.get("min_support", 0.5))

    pa = pd.read_csv(args.rtab, sep="\t", index_col=0)
    pa.columns = [c.split("/")[-1].replace(".gff", "").replace(".fna", "") for c in pa.columns]
    nm = pd.read_csv(args.niche_map, sep="\t")
    nm = nm[nm["role"] == "focal"]
    pops = {n: [g for g in d["genome"] if g in pa.columns]
            for n, d in nm.groupby("niche")}
    pops = {n: v for n, v in pops.items() if len(v) >= 3}
    niches = sorted(pops)

    rows, summ = [], []
    for X, Y in itertools.combinations(niches, 2):
        gx, gy = pops[X], pops[Y]
        nX, nY = len(gx), len(gy)
        sub = pa[gx + gy]
        # drop core (present in all) and absent-everywhere genes
        psum = sub.sum(axis=1)
        acc = sub[(psum > 0) & (psum < (nX + nY))]
        if acc.empty:
            continue

        # ---- full data: exact test, unequal n (reported, not counted) --------
        AX = acc[gx].to_numpy(int)
        AY = acc[gy].to_numpy(int)
        a_full, c_full = AX.sum(axis=1), AY.sum(axis=1)
        p_full = np.array([fisher_exact([[a, nX - a], [c, nY - c]])[1]
                           for a, c in zip(a_full, c_full)])
        log2or = np.log2(((a_full + 0.5) * (nY - c_full + 0.5)) /
                         ((nX - a_full + 0.5) * (c_full + 0.5)))
        rd = pd.DataFrame({"pair": f"{X}|{Y}", "gene": acc.index,
                           "present_X": a_full, "present_Y": c_full,
                           "nX": nX, "nY": nY, "log2_or": log2or, "p": p_full})
        rd["q"] = multipletests(rd["p"], method="fdr_bh")[1]
        rd["enriched_in"] = np.where(rd["log2_or"] > 0, X, Y)

        # ---- balanced draws: equal n, so the counts are comparable -----------
        n_min = min(nX, nY)
        sup_x = np.zeros(len(acc)); sup_y = np.zeros(len(acc))
        for b in range(B):
            rng = np.random.default_rng(derive_seed(cfg["seed"], "accessory", X, Y, b))
            a = AX[:, rng.choice(nX, n_min, replace=False)].sum(axis=1)
            c = AY[:, rng.choice(nY, n_min, replace=False)].sum(axis=1)
            p = balanced_two_sided_p(a, c, n_min)
            q = multipletests(p, method="fdr_bh")[1]
            sig = q < fdr
            sup_x += sig & (a > c)
            sup_y += sig & (c > a)
        rd["support_X"] = sup_x / B
        rd["support_Y"] = sup_y / B
        rd["balanced_call"] = np.where(rd["support_X"] >= min_support, X,
                               np.where(rd["support_Y"] >= min_support, Y, "ns"))
        rows.append(rd)

        summ.append({
            "pair": f"{X}|{Y}", "n_accessory_tested": len(rd),
            "n_strains_X": nX, "n_strains_Y": nY, "n_balanced_subsample": n_min,
            "n_balanced_iterations": B, "min_support": min_support,
            # balanced counts: what transition_verdict.py consumes
            f"enriched_in_{X}": int((rd["balanced_call"] == X).sum()),
            f"enriched_in_{Y}": int((rd["balanced_call"] == Y).sum()),
            "n_diff_q05": int((rd["balanced_call"] != "ns").sum()),
            # full-data counts, retained for reference; power differs with n
            f"unbalanced_enriched_in_{X}": int(((rd["q"] < fdr) & (rd["enriched_in"] == X)).sum()),
            f"unbalanced_enriched_in_{Y}": int(((rd["q"] < fdr) & (rd["enriched_in"] == Y)).sum()),
            "unbalanced_n_diff_q05": int((rd["q"] < fdr).sum())})

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["pair", "gene", "present_X", "present_Y", "nX", "nY", "log2_or",
                 "p", "q", "enriched_in", "support_X", "support_Y", "balanced_call"])
    if len(out):
        out = out.sort_values(["pair", "q"])
    out.to_csv(args.out, sep="\t", index=False)
    pd.DataFrame(summ).to_csv(args.summary, sep="\t", index=False)
    log.info("Accessory differentiation: %d pairs, %d gene tests, %d balanced calls",
             len(summ), len(out), int((out["balanced_call"] != "ns").sum()) if len(out) else 0)


if __name__ == "__main__":
    main()
