#!/usr/bin/env python3
"""
functional_redundancy.py
------------------------
How functionally redundant is the human gut: are functions carried by many
species (robust) or few (fragile)? Species is the unit throughout, so over-
sequenced species do not inflate redundancy, and comparisons are rarefied.

Per functional layer, for the human gut (and optionally within Western /
non-Western):
  * occupancy        per function, how many species carry it; the fraction that
                     are core (>= core_fraction of species) vs rare (<= 1 species).
  * accumulation     rarefied curve of unique functions as species are added; a
                     curve that saturates while species keep accumulating is the
                     signature of redundancy. (The species curve is y=k.)
  * Ricotta FR       community functional redundancy = Gini-Simpson diversity D
                     minus Rao's quadratic entropy Q on functional distances;
                     FR=D-Q, relative FR = 1-Q/D. Computed on a species subsample
                     with bootstraps (pairwise distances are O(S^2)).
  * carrier spread   per function, the number of distinct families carrying it
                     (a phylogenetic-breadth proxy: robust vs single-clade).

Outputs (per layer + a combined summary): occupancy_<layer>.tsv,
accumulation_<layer>.tsv, spread_<layer>.tsv, redundancy_summary.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("redundancy")


def occupancy(present_long, species_set):
    p = present_long[present_long["species"].isin(species_set)]
    occ = p.groupby("feature")["species"].nunique()
    return occ


def accumulation(present_long, species_list, steps, nboot, seed, tag):
    sl = np.array(sorted(species_list))
    S = len(sl)
    if S < 3:
        return pd.DataFrame(columns=["k", "functions_mean", "lo", "hi"])
    ks = np.unique(np.linspace(1, S, steps).astype(int))
    feat_by_sp = present_long.groupby("species")["feature"].apply(set).to_dict()
    rows = []
    for k in ks:
        vals = []
        for b in range(nboot):
            rng = np.random.default_rng(derive_seed(seed, "accum", tag, k, b))
            pick = rng.choice(sl, size=k, replace=False)
            u = set()
            for s in pick:
                u |= feat_by_sp.get(s, set())
            vals.append(len(u))
        vals = np.array(vals)
        rows.append({"k": int(k), "functions_mean": vals.mean(),
                     "lo": np.percentile(vals, 2.5), "hi": np.percentile(vals, 97.5)})
    return pd.DataFrame(rows)


def ricotta_fr(present_long, species_list, cfg_r, seed, tag):
    """Gini-Simpson D, Rao Q, FR=D-Q on a species subsample (bootstrapped)."""
    maxs = cfg_r.get("max_species", 1000)
    nboot = 5
    sl = sorted(species_list)
    feat_by_sp = present_long.groupby("species")["feature"].apply(set).to_dict()
    Ds, Qs, FRs = [], [], []
    for b in range(nboot):
        rng = np.random.default_rng(derive_seed(seed, "ricotta", tag, b))
        sub = list(rng.choice(sl, size=min(maxs, len(sl)), replace=False)) if len(sl) > maxs else sl
        feats = sorted(set().union(*[feat_by_sp.get(s, set()) for s in sub]) or {"_"})
        fidx = {f: i for i, f in enumerate(feats)}
        M = np.zeros((len(sub), len(feats)), dtype=bool)
        for i, s in enumerate(sub):
            for f in feat_by_sp.get(s, set()):
                M[i, fidx[f]] = True
        if M.shape[0] < 3:
            continue
        d = squareform(pdist(M, metric=cfg_r.get("distance", "jaccard")))
        p = np.full(M.shape[0], 1.0 / M.shape[0])     # equal species weights
        D = 1.0 - np.sum(p ** 2)
        Q = float(p @ d @ p)
        Ds.append(D); Qs.append(Q); FRs.append(D - Q)
    if not Ds:
        return dict(ricotta_D=np.nan, ricotta_Q=np.nan, ricotta_FR=np.nan, ricotta_relFR=np.nan)
    D, Q, FR = np.mean(Ds), np.mean(Qs), np.mean(FRs)
    return dict(ricotta_D=round(D, 4), ricotta_Q=round(Q, 4), ricotta_FR=round(FR, 4),
                ricotta_relFR=round(1 - Q / D, 4) if D > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rc = cfg["redundancy"]
    seed = cfg["seed"]
    rank = rc["phylo_redundancy_rank"]
    os.makedirs(args.out_dir, exist_ok=True)

    sp = pd.read_csv(args.species_table, sep="\t")
    fam_of = dict(zip(sp["species"], sp.get(rank, pd.Series(index=sp.index, dtype=str))))
    gen_of = dict(zip(sp["species"], sp.get("gtdb_genus", pd.Series(index=sp.index, dtype=str))))

    # communities: human gut, and optionally within Western / non-Western
    comms = {"human": set(sp.loc[sp["n_human"] > 0, "species"])}
    if rc.get("by_population") and "n_western" in sp.columns:
        w = set(sp.loc[sp["n_western"] > 0, "species"])
        nw = set(sp.loc[sp["n_non_western"] > 0, "species"])
        if w:
            comms["western"] = w
        if nw:
            comms["non_western"] = nw

    summary = []
    for layer in rc["layers"]:
        path = f"{args.profiles_dir}/prevalence_{layer}.parquet"
        if not os.path.exists(path):
            continue
        prev = pd.read_parquet(path)
        prev = prev[prev["present"] == 1][["species", "feature"]]

        for comm, sset in comms.items():
            pl = prev[prev["species"].isin(sset)]
            if pl.empty:
                continue
            # One consistent denominator everywhere: the species of this community
            # that carry at least one feature in this layer (i.e. are annotated for
            # it). Species with no annotation in a layer are excluded from that
            # layer's occupancy, accumulation and Ricotta FR rather than being
            # counted as "lacking every function", which would inflate apparent
            # complementarity. n_species_pool vs n_species is reported.
            annot = sorted(set(pl["species"]))
            S = len(annot)
            if S < 3:
                continue
            occ = occupancy(pl, set(annot))
            tag = f"{layer}_{comm}"
            # occupancy table + classes (only write for the main human community)
            if comm == "human":
                occdf = occ.rename("n_species").reset_index()
                occdf["prop_species"] = occdf["n_species"] / S
                occdf["class"] = np.where(occdf["prop_species"] >= rc["core_fraction"], "core",
                                  np.where(occdf["n_species"] <= rc["rare_max_species"], "rare",
                                           "intermediate"))
                occdf.sort_values("n_species", ascending=False).to_csv(
                    f"{args.out_dir}/occupancy_{layer}.tsv", sep="\t", index=False)
                accumulation(pl, annot, rc["accumulation"]["steps"],
                             rc["accumulation"]["bootstrap"], seed, tag).to_csv(
                    f"{args.out_dir}/accumulation_{layer}.tsv", sep="\t", index=False)
                # carrier spread (families / genera per function)
                car = pl.copy()
                car["family"] = car["species"].map(fam_of)
                car["genus"] = car["species"].map(gen_of)
                spread = car.groupby("feature").agg(
                    n_species=("species", "nunique"),
                    n_families=("family", "nunique"),
                    n_genera=("genus", "nunique")).reset_index()
                spread.sort_values("n_families", ascending=False).to_csv(
                    f"{args.out_dir}/spread_{layer}.tsv", sep="\t", index=False)

            ric = ricotta_fr(pl, annot, rc["ricotta"], seed, tag)
            prop = occ / S
            summary.append({"layer": layer, "community": comm,
                            "n_species_pool": len(sset), "n_species": S,
                            "n_features": int(occ.shape[0]),
                            "median_occupancy": float(occ.median()),
                            "pct_core": float((prop >= rc["core_fraction"]).mean() * 100),
                            "pct_rare": float((occ <= rc["rare_max_species"]).mean() * 100),
                            **ric})
            log.info("%s/%s: %d species, %d features, relFR=%s",
                     layer, comm, S, occ.shape[0], ric["ricotta_relFR"])

    pd.DataFrame(summary).to_csv(f"{args.out_dir}/redundancy_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
