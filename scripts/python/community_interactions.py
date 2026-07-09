#!/usr/bin/env python3
"""
community_interactions.py
-------------------------
Community-level (niche-pool) metabolic-interaction analyses on the species-
collapsed data. The species-per-niche confounder is controlled by rarefying every
niche's species pool to the same number of species (bootstrapped); the strains-
per-species confounder was already removed at ingest (species-level prevalence).

Per niche (species occupying the niche), rarefied to equal species number:
  * auxotrophy dependency      mean auxotrophies per species; fraction auxotrophic;
                               per amino acid the prototroph:auxotroph ratio (the
                               community's capacity to provision auxotrophs).
  * carbon cross-feeding       for each exchange metabolite, potential
                               producer -> consumer links; number of active
                               exchanges and total cross-feeding potential. Acetate
                               uses a real producer module (M00579); other organic
                               acids use a fermenter-trait producer proxy (map).
  * module complementarity     collective module coverage vs mean per species
                               (division of labour) of the pool.
  * trait composition          fraction of the pool that is anaerobe / spore-
                               forming / bile-susceptible / motile / ...

Outputs: community_summary.tsv, auxotrophy_by_aa.tsv,
crossfeeding_by_metabolite.tsv, trait_composition.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("community")


def present_matrix(path):
    """species x feature 0/1 matrix from a prevalence parquet."""
    d = pd.read_parquet(path)
    m = (d[d["present"] == 1].pivot_table(index="species", columns="feature",
                                          values="present", fill_value=0)
         .astype(np.int8))
    return m


def ci(a):
    a = np.asarray(a, float)
    return float(np.nanmean(a)), float(np.nanpercentile(a, 2.5)), float(np.nanpercentile(a, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--metabolite-map", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cc = cfg["community"]
    seed = cfg["seed"]
    niches = cfg["inputs"]["niche_levels"]
    B = cc["rarefaction"]["bootstrap"]
    os.makedirs(args.out_dir, exist_ok=True)

    sp = pd.read_csv(args.species_table, sep="\t")
    pools = {n: set(sp.loc[sp[f"n_{n}"] > 0, "species"]) for n in niches
             if f"n_{n}" in sp.columns}
    pools = {n: v for n, v in pools.items() if len(v) >= 10}
    if not pools:
        log.warning("No niche pools; abort."); return
    k = min(len(v) for v in pools.values())

    aux = present_matrix(f"{args.profiles_dir}/prevalence_auxotrophy.parquet")   # present=auxotroph
    carb = present_matrix(f"{args.profiles_dir}/prevalence_carbon.parquet")
    mods = present_matrix(f"{args.profiles_dir}/prevalence_modulep.parquet")
    trait = present_matrix(f"{args.profiles_dir}/prevalence_trait.parquet")

    mmap = pd.read_csv(args.metabolite_map, sep="\t", dtype=str, keep_default_na=False)

    def sub(mat, species):
        """Rows for exactly these species, zero-filled.

        A species with no 'present' feature is absent from the pivoted matrix.
        Dropping it would compute means and fractions over only the species that
        have >=1 feature, inflating (e.g.) mean auxotrophies and trait fractions.
        Reindexing keeps the denominator equal to the rarefied pool size.
        """
        return mat.reindex(index=list(species), fill_value=0)

    summ, aa_rows, cf_rows, tr_rows = [], [], [], []
    for niche, pool in pools.items():
        pool = sorted(pool)
        # accumulate bootstrap draws
        acc = {m: [] for m in ["mean_aux", "frac_aux_sp", "carbon_breadth",
                               "n_active_exchange", "crossfeed_potential",
                               "mod_collective", "mod_mean_per_sp", "mod_complementarity"]}
        aa_frac = {a: [] for a in aux.columns}
        aa_ratio = {a: [] for a in aux.columns}
        cf_prod = {m: [] for m in mmap["metabolite"]}
        cf_cons = {m: [] for m in mmap["metabolite"]}
        cf_pot = {m: [] for m in mmap["metabolite"]}
        tr_frac = {t: [] for t in trait.columns}

        for b in range(B):
            rng = np.random.default_rng(derive_seed(seed, "community", niche, b))
            pick = list(rng.choice(pool, size=k, replace=False))
            A = sub(aux, pick); C = sub(carb, pick); M = sub(mods, pick); T = sub(trait, pick)

            # auxotrophy
            if len(A):
                per = A.sum(axis=1)
                acc["mean_aux"].append(per.mean())
                acc["frac_aux_sp"].append((per > 0).mean())
                for a in aux.columns:
                    n_aux = int(A[a].sum()); n_pro = len(A) - n_aux
                    aa_frac[a].append(n_aux / len(A))
                    aa_ratio[a].append(n_pro / n_aux if n_aux > 0 else np.nan)
            # carbon breadth
            if len(C):
                acc["carbon_breadth"].append(C.sum(axis=1).mean())
            # cross-feeding
            n_active, tot_pot = 0, 0.0
            for _, r in mmap.iterrows():
                met = r["metabolite"]; cons_col = r["consumer_carbon"]
                consumers = int(C[cons_col].sum()) if (len(C) and cons_col in C.columns) else 0
                if r["producer_module"] and len(M) and r["producer_module"] in M.columns:
                    producers = int(M[r["producer_module"]].sum())
                elif r["producer_trait"] and len(T) and r["producer_trait"] in T.columns:
                    producers = int(T[r["producer_trait"]].sum())
                else:
                    producers = 0
                pot = producers * consumers / (k * k)
                cf_prod[met].append(producers); cf_cons[met].append(consumers); cf_pot[met].append(pot)
                if producers > 0 and consumers > 0:
                    n_active += 1; tot_pot += pot
            acc["n_active_exchange"].append(n_active)
            acc["crossfeed_potential"].append(tot_pot)
            # module complementarity
            if len(M):
                collective = int((M.sum(axis=0) > 0).sum())
                mean_sp = M.sum(axis=1).mean()
                acc["mod_collective"].append(collective)
                acc["mod_mean_per_sp"].append(mean_sp)
                acc["mod_complementarity"].append(collective / mean_sp if mean_sp > 0 else np.nan)
            # traits
            if len(T):
                for t in trait.columns:
                    tr_frac[t].append(T[t].mean())

        row = {"niche": niche, "n_species_pool": len(pool), "rarefied_to": k}
        for m, vals in acc.items():
            mean, lo, hi = ci(vals) if vals else (np.nan, np.nan, np.nan)
            row[m] = round(mean, 4); row[f"{m}_lo"] = round(lo, 4); row[f"{m}_hi"] = round(hi, 4)
        summ.append(row)
        for a in aux.columns:
            aa_rows.append({"niche": niche, "amino_acid": a,
                            "frac_auxotroph": round(np.nanmean(aa_frac[a]), 4),
                            "prototroph_auxotroph_ratio": round(np.nanmean(aa_ratio[a]), 3)})
        for met in mmap["metabolite"]:
            cf_rows.append({"niche": niche, "metabolite": met,
                            "n_producers": round(np.nanmean(cf_prod[met]), 1),
                            "n_consumers": round(np.nanmean(cf_cons[met]), 1),
                            "potential_norm": round(np.nanmean(cf_pot[met]), 5)})
        for t in trait.columns:
            m, lo, hi = ci(tr_frac[t]) if tr_frac[t] else (np.nan, np.nan, np.nan)
            tr_rows.append({"niche": niche, "trait": t, "fraction": round(m, 4),
                            "lo": round(lo, 4), "hi": round(hi, 4)})

    pd.DataFrame(summ).to_csv(f"{args.out_dir}/community_summary.tsv", sep="\t", index=False)
    pd.DataFrame(aa_rows).to_csv(f"{args.out_dir}/auxotrophy_by_aa.tsv", sep="\t", index=False)
    pd.DataFrame(cf_rows).to_csv(f"{args.out_dir}/crossfeeding_by_metabolite.tsv", sep="\t", index=False)
    pd.DataFrame(tr_rows).to_csv(f"{args.out_dir}/trait_composition.tsv", sep="\t", index=False)
    log.info("Community interactions done: niches=%s rarefied to %d species",
             list(pools), k)


if __name__ == "__main__":
    main()
