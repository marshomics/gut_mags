#!/usr/bin/env python3
"""
permutation_null.py
-------------------
Empirical calibration of the signature count. Group labels are permuted across
species and the fast enrichment screen (prevalence-difference with a Fisher
exact test) is recomputed, giving a null distribution for the number of
"significant" features. Comparing the observed count to this null yields an
empirical false-discovery estimate that does not rely on the parametric
assumptions of any single method.

Permutation is PHYLOGENY-AWARE: labels are shuffled WITHIN phylum (configurable
rank) so the null preserves the phylogenetic clustering of niche membership.
Shuffling labels globally would destroy that structure and make the null too
easy to beat, overstating significance; within-phylum shuffling is the
conservative choice.

Output JSON:
  observed_n_sig, null_mean, null_p95, empirical_fdr, n_permutations
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("permnull")


def screen(M, group, alpha, eff_log2):
    """Return indices of features 'significant' by Fisher + effect size."""
    n1 = (group == 1).sum(); n0 = (group == 0).sum()
    p1 = M[group == 1].mean(axis=0); p0 = M[group == 0].mean(axis=0)
    a = p1 * n1 + 0.5; b = (1 - p1) * n1 + 0.5
    c = p0 * n0 + 0.5; d = (1 - p0) * n0 + 0.5
    log2or = np.log2((a * d) / (b * c))
    sig = np.zeros(M.shape[1], bool)
    cand = np.where(np.abs(log2or) >= eff_log2)[0]
    for j in cand:
        a1 = int(round(p1[j] * n1)); a0 = int(round(p0[j] * n0))
        tbl = [[a1, n1 - a1], [a0, n0 - a0]]
        _, pv = fisher_exact(tbl)
        if pv < alpha:
            sig[j] = True
    return sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rank", default="family")
    args = ap.parse_args()

    cfg = load_config(args.config)
    if not cfg["stats"]["permutation_null"]["enabled"]:
        json.dump({"enabled": False}, open(args.out, "w")); return
    B = cfg["stats"]["permutation_null"]["iterations"]
    alpha = cfg["stats"]["fdr_alpha"]
    eff = cfg["stats"]["min_abs_log2_or"]

    meta = pd.read_csv(args.analysis, sep="\t")
    rank = args.rank if args.rank in meta.columns else "family"
    meta[rank] = meta[rank].fillna("unknown").replace("", "unknown")
    pres = pd.read_parquet(args.presence)
    wide = (pres.pivot_table(index="species", columns="feature",
                             values="present", fill_value=0).astype(np.int8))
    meta = meta[meta["species"].isin(wide.index)].copy()
    wide = wide.loc[meta["species"]]
    M = wide.to_numpy()
    group = meta["group"].to_numpy().astype(int)
    strata = meta[rank].to_numpy()

    obs = int(screen(M, group, alpha, eff).sum())

    # Labels are shuffled only WITHIN a clade, so the clade-niche association is
    # preserved and the null is conservative. A stratum holding a single group
    # cannot contribute: its labels are invariant under permutation. If almost no
    # species are permutable the null is degenerate and its empirical FDR is
    # uninformative, so both counts are reported and a warning is raised.
    informative = [s for s in np.unique(strata) if len(set(group[strata == s])) > 1]
    n_permutable = int(sum((strata == s).sum() for s in informative))

    null_counts = []
    for b in range(B):
        rng = np.random.default_rng(derive_seed(cfg["seed"], "perm", args.contrast, b))
        perm = group.copy()
        for s in np.unique(strata):          # shuffle within phylum/family
            idx = np.where(strata == s)[0]
            perm[idx] = rng.permutation(group[idx])
        null_counts.append(int(screen(M, perm, alpha, eff).sum()))
    null_counts = np.array(null_counts)

    # empirical FDR is a proportion, so cap it at 1; and report the one-sided
    # permutation p for "more signal than the null", with the add-one correction.
    emp_fdr = float(min(1.0, null_counts.mean() / obs)) if obs > 0 else np.nan
    p_emp = float((1 + (null_counts >= obs).sum()) / (B + 1))
    out = {
        "contrast": args.contrast, "n_permutations": B,
        "observed_n_sig": obs,
        "null_mean": float(null_counts.mean()),
        "null_p95": float(np.percentile(null_counts, 95)),
        "empirical_fdr": emp_fdr,
        "p_empirical": p_emp,
        "permutation_rank": rank,
        "n_species": int(len(group)),
        "n_informative_strata": len(informative),
        "n_permutable_species": n_permutable,
        "null_degenerate": bool(n_permutable < 0.1 * len(group)),
    }
    json.dump(out, open(args.out, "w"), indent=2)
    log.info("Permutation null: observed=%d null_mean=%.1f emp_FDR=%.3f p=%.4g "
             "(%d/%d species permutable across %d informative strata)",
             obs, out["null_mean"], emp_fdr if emp_fdr == emp_fdr else float("nan"),
             p_emp, n_permutable, len(group), len(informative))
    if out["null_degenerate"]:
        log.warning("Fewer than 10%% of species are permutable at rank '%s': niche is "
                    "nearly nested within clade, so the within-clade null cannot move "
                    "labels and its empirical FDR is not interpretable. Use a coarser "
                    "rank (--rank phylum), or read this null as a lower bound.", rank)


if __name__ == "__main__":
    main()
