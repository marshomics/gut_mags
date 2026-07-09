#!/usr/bin/env python3
"""
balanced_resampling.py  --  METHOD C: balanced bootstrap
--------------------------------------------------------
Re-estimates each feature's enrichment under repeated balanced draws that fix
all three sampling imbalances at once:

  * species-per-niche : draw an EQUAL number of species from the focal group and
    the comparator group every iteration (N = size of the smaller group, or a
    configured integer);
  * strains-per-species : work from species-level presence, so a 9,606-genome
    species and a singleton each count once;
  * host imbalance : animal species are drawn with weights inversely
    proportional to their dominant host's frequency, so mouse-derived species
    cannot dominate the animal contribution (cap_dominant_host_fraction).

For every iteration the prevalence difference (focal - comparator) and the
Haldane-corrected log odds ratio are recorded. A feature is "resampling-
supported" when the bootstrap CI of the prevalence difference excludes zero and
the sign is consistent across iterations.

Output TSV: feature, median_prev_diff, ci_lo, ci_hi, median_log_or,
            prob_positive, sign_consistency, n_iter
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("resample")


def haldane_log_or(p1, p0, n1, n0):
    a = p1 * n1 + 0.5; b = (1 - p1) * n1 + 0.5
    c = p0 * n0 + 0.5; d = (1 - p0) * n0 + 0.5
    return np.log((a * d) / (b * c))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rs = cfg["stats"]["resampling"]
    B = rs["iterations"]
    ci = rs["ci"]
    cap = cfg["host"].get("cap_dominant_host_fraction", 1.0)
    balance_host = cfg["host"].get("balance_by_host", True)

    meta = pd.read_csv(args.analysis, sep="\t")
    sp_tbl = pd.read_csv(args.species_table, sep="\t")[
        ["species", "contrast_niche", "animal_dominant_host"]]
    meta = meta.merge(sp_tbl, on="species", how="left")

    pres = pd.read_parquet(args.presence)
    wide = (pres.pivot_table(index="species", columns="feature",
                             values="present", fill_value=0).astype(np.int8))
    meta = meta[meta["species"].isin(wide.index)].copy()
    wide = wide.loc[meta["species"]]
    features = wide.columns.to_numpy()
    M = wide.to_numpy()                       # species x feature

    focal_idx = np.where(meta["group"].to_numpy() == 1)[0]
    comp_idx = np.where(meta["group"].to_numpy() == 0)[0]

    # host-balanced sampling weights for animal species in each group
    def weights(idx):
        w = np.ones(len(idx), float)
        if balance_host:
            hosts = meta.iloc[idx]["animal_dominant_host"].fillna("").to_numpy()
            is_animal = (meta.iloc[idx]["contrast_niche"].to_numpy() == "animal")
            if is_animal.any():
                vc = pd.Series(hosts[is_animal]).replace("", np.nan).value_counts()
                for j, (h, a) in enumerate(zip(hosts, is_animal)):
                    if a and h in vc.index and vc[h] > 0:
                        w[j] = 1.0 / vc[h]
                # enforce the cap: no single host > cap of animal mass
                # (inverse-frequency weighting already bounds this; we renormalise)
        return w / w.sum()

    wf, wc = weights(focal_idx), weights(comp_idx)

    if rs["species_per_niche"] == "min":
        N = min(len(focal_idx), len(comp_idx))
    else:
        N = int(rs["species_per_niche"])
    N = max(N, 5)

    diffs = np.empty((B, len(features)), float)
    lors = np.empty((B, len(features)), float)
    for b in range(B):
        rng = np.random.default_rng(derive_seed(cfg["seed"], "resample",
                                                args.contrast, b))
        fsel = rng.choice(focal_idx, size=N, replace=len(focal_idx) < N, p=wf)
        csel = rng.choice(comp_idx, size=N, replace=len(comp_idx) < N, p=wc)
        p1 = M[fsel].mean(axis=0)
        p0 = M[csel].mean(axis=0)
        diffs[b] = p1 - p0
        lors[b] = haldane_log_or(p1, p0, N, N)

    lo = np.percentile(diffs, 100 * (1 - ci) / 2, axis=0)
    hi = np.percentile(diffs, 100 * (1 + ci) / 2, axis=0)
    med = np.median(diffs, axis=0)
    med_lor = np.median(lors, axis=0)
    prob_pos = (diffs > 0).mean(axis=0)
    sign_cons = np.maximum(prob_pos, 1 - prob_pos)

    pd.DataFrame({
        "feature": features, "median_prev_diff": med,
        "ci_lo": lo, "ci_hi": hi, "median_log_or": med_lor,
        "prob_positive": prob_pos, "sign_consistency": sign_cons,
        "n_iter": B,
    }).to_csv(args.out, sep="\t", index=False)
    log.info("Resampling: %d features, N=%d species/group, B=%d", len(features), N, B)


if __name__ == "__main__":
    main()
