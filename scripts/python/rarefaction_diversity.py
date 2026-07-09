#!/usr/bin/env python3
"""
rarefaction_diversity.py
------------------------
Compare taxonomic richness across niches WITHOUT being fooled by the fact that
the niches are sampled to very different depths (human 6k species vs free 46k).

Two complementary outputs:

  1. Sample-size-based rarefaction of species richness: repeatedly subsample an
     equal number of GENOMES from each niche and count distinct species, so
     richness is compared at matched sampling effort. Bootstrap CIs included.

  2. Hill numbers (q = 0, 1, 2) on the species-abundance (genomes-per-species)
     distribution within each niche, with coverage estimated via the Chao
     sample-coverage estimator. Hill q=0 is richness, q=1 is exp(Shannon),
     q=2 is inverse Simpson; together they describe how evenly genomes are
     spread across species in each niche.

Why this matters: a raw "free-living has more species" statement is an artefact
of sampling. Rarefied/Hill comparison is the defensible version.

Outputs:
  rarefaction_curves.tsv   niche, effort(genomes), richness_mean, lo, hi
  hill_numbers.tsv         niche, q, hill_mean, lo, hi, coverage
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import (load_config, get_logger, set_global_seed, derive_seed)

log = get_logger("diversity")


def hill(counts, q):
    counts = np.asarray(counts, float)
    p = counts / counts.sum()
    p = p[p > 0]
    if q == 1:
        return float(np.exp(-(p * np.log(p)).sum()))
    return float((np.power(p, q).sum()) ** (1.0 / (1.0 - q)))


def chao1(counts):
    """Bias-corrected Chao1 asymptotic richness + lognormal 95% CI."""
    counts = np.asarray(counts, float)
    S = int((counts > 0).sum())
    f1 = int((counts == 1).sum()); f2 = int((counts == 2).sum())
    est = S + (f1 * (f1 - 1)) / (2 * (f2 + 1))
    # variance (Chao 1987) and lognormal CI
    if f1 > 0:
        var = (f1 * (f1 - 1) / (2 * (f2 + 1))
               + f1 * (2 * f1 - 1) ** 2 / (4 * (f2 + 1) ** 2)
               + f1 ** 2 * f2 * (f1 - 1) ** 2 / (4 * (f2 + 1) ** 4))
    else:
        var = 0.0
    d = est - S
    if d > 0 and var > 0:
        C = np.exp(1.96 * np.sqrt(np.log(1 + var / d ** 2)))
        lo, hi = S + d / C, S + d * C
    else:
        lo = hi = est
    return est, lo, hi


def ace(counts, rare_cutoff=10):
    """Abundance-based Coverage Estimator of richness."""
    counts = np.asarray(counts, float)
    S_abund = int((counts > rare_cutoff).sum())
    rare = counts[(counts > 0) & (counts <= rare_cutoff)]
    S_rare = int(len(rare))
    N_rare = rare.sum()
    f1 = int((counts == 1).sum())
    if N_rare == 0:
        return float(S_abund + S_rare)
    C_ace = 1 - f1 / N_rare if N_rare > 0 else 1
    if C_ace == 0:
        return float(S_abund + S_rare)
    g = 0.0
    denom = sum(i * (i - 1) * (counts == i).sum() for i in range(1, rare_cutoff + 1))
    gamma = max((S_rare / C_ace) * denom / (N_rare * (N_rare - 1)) - 1, 0) if N_rare > 1 else 0
    return float(S_abund + S_rare / C_ace + (f1 / C_ace) * gamma)


def chao_coverage(counts):
    """Good-Turing / Chao sample coverage estimate C_hat."""
    counts = np.asarray(counts, float)
    n = counts.sum()
    f1 = int((counts == 1).sum())
    f2 = int((counts == 2).sum())
    if n == 0:
        return np.nan
    if f2 > 0:
        return 1.0 - (f1 / n) * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    if f1 > 0:
        return 1.0 - (f1 / n) * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    return 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-steps", type=int, default=25)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    niches = cfg["inputs"]["niche_levels"]

    df = pd.read_parquet(args.samples)[["genome", "species", "niche"]]

    # --- rarefaction (sample-size based, matched effort) -----------------------
    min_genomes = df.groupby("niche").size().min()
    efforts = np.unique(np.linspace(1, min_genomes, args.n_steps).astype(int))
    rows = []
    for n in niches:
        sp = df.loc[df["niche"] == n, "species"].to_numpy()
        for e in efforts:
            vals = []
            for b in range(args.n_boot):
                rng = np.random.default_rng(derive_seed(cfg["seed"], "rare", n, e, b))
                pick = rng.choice(sp, size=e, replace=False)
                vals.append(len(set(pick)))
            vals = np.array(vals)
            rows.append({"niche": n, "effort": int(e),
                         "richness_mean": vals.mean(),
                         "lo": np.percentile(vals, 2.5),
                         "hi": np.percentile(vals, 97.5)})
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}_rarefaction_curves.tsv",
                              sep="\t", index=False)
    log.info("Rarefaction done; matched effort up to %d genomes", min_genomes)

    # --- Hill numbers with bootstrap -------------------------------------------
    hrows = []
    for n in niches:
        counts = df.loc[df["niche"] == n].groupby("species").size().to_numpy()
        cov = chao_coverage(counts)
        for q in (0, 1, 2):
            point = hill(counts, q)
            boots = []
            for b in range(args.n_boot):
                rng = np.random.default_rng(derive_seed(cfg["seed"], "hill", n, q, b))
                # multinomial resample of the genome->species assignment
                resampled = rng.multinomial(counts.sum(), counts / counts.sum())
                boots.append(hill(resampled[resampled > 0], q))
            boots = np.array(boots)
            hrows.append({"niche": n, "q": q, "hill_mean": point,
                          "lo": np.percentile(boots, 2.5),
                          "hi": np.percentile(boots, 97.5),
                          "coverage": cov})
    pd.DataFrame(hrows).to_csv(f"{args.out_prefix}_hill_numbers.tsv",
                               sep="\t", index=False)
    log.info("Hill numbers done.")

    # --- asymptotic richness estimators (Chao1, ACE) ---------------------------
    erows = []
    for n in niches:
        counts = df.loc[df["niche"] == n].groupby("species").size().to_numpy()
        c, lo, hi = chao1(counts)
        erows.append({"niche": n, "S_observed": int((counts > 0).sum()),
                      "chao1": round(c, 1), "chao1_lo": round(lo, 1),
                      "chao1_hi": round(hi, 1), "ace": round(ace(counts), 1),
                      "coverage": round(chao_coverage(counts), 4)})
    pd.DataFrame(erows).to_csv(f"{args.out_prefix}_richness_estimators.tsv",
                               sep="\t", index=False)
    log.info("Richness estimators done: %s",
             {r["niche"]: r["chao1"] for r in erows})


if __name__ == "__main__":
    main()
