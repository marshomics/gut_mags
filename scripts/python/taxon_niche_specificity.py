#!/usr/bin/env python3
"""
taxon_niche_specificity.py
--------------------------
How niche-specific is each taxon, at every rank? For each taxon (phylum down to
genus) this computes its distribution of species across niches, its niche
breadth, and whether it is more niche-restricted than expected by chance.

Two views are always produced so the strain-sampling bias is visible:
  * species-weighted: distinct species of the taxon occupying each niche (primary);
  * genome-weighted: genomes of the taxon in each niche (the biased view).

Specificity test: a taxon with k species is compared against a null in which
species' niche labels are shuffled (preserving niche sizes). The taxon's
observed standardised Levins breadth is compared to the null distribution; a
small p means the taxon's species cluster into one niche more tightly than
random. p-values are FDR-corrected within each rank. This converts "family X
looks human-specific" into a test, which is what a reviewer will want.

Outputs (one per rank, plus a summary):
  taxon_specificity_<rank>.tsv
  specialist_summary.tsv
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, set_global_seed
from statsmodels.stats.multitest import multipletests

log = get_logger("taxonspec")


def breadth_std(counts_2d, K):
    """Standardised Levins breadth per row of a (taxa x niche) count matrix."""
    tot = counts_2d.sum(axis=1, keepdims=True)
    p = np.divide(counts_2d, tot, out=np.zeros_like(counts_2d, float), where=tot > 0)
    B = 1.0 / np.clip((p ** 2).sum(axis=1), 1e-12, None)
    return (B - 1.0) / (K - 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    tcfg = cfg["taxonomy"]
    niches = cfg["inputs"]["niche_levels"]
    K = len(niches)
    N = tcfg["null_iterations"]
    min_sp = tcfg["enrichment"]["min_taxon_species"]
    cutoff = tcfg["specialist_breadth_cutoff"]

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    samples = pd.read_parquet(args.samples)
    sptab = pd.read_csv(args.species_table, sep="\t")
    # one row per species: taxonomy (all ranks) + primary niche
    rankcols = tcfg["ranks"]
    sp_tax = samples.groupby("species")[rankcols].first()
    sp_prim = sptab.set_index("species")["niche_primary"]
    sp = sp_tax.join(sp_prim).dropna(subset=["niche_primary"])
    sp["niche_primary"] = pd.Categorical(sp["niche_primary"], categories=niches)

    # species occupancy (distinct species x niche) for species-weighted counts
    occ_pairs = samples[["species", "niche"]].drop_duplicates()

    summary = []
    for rank in rankcols:
        if rank == "gtdb_species":
            continue  # specificity is degenerate at species level (use occupancy elsewhere)
        # species-weighted: distinct species of taxon occupying each niche
        occ_tax = occ_pairs.merge(sp_tax[[rank]], left_on="species", right_index=True)
        sw = (occ_tax.groupby([rank, "niche"])["species"].nunique()
              .unstack(fill_value=0).reindex(columns=niches, fill_value=0))
        # genome-weighted
        gw = (samples.groupby([rank, "niche"]).size()
              .unstack(fill_value=0).reindex(columns=niches, fill_value=0))
        n_species = sp.groupby(rank).size().reindex(sw.index).fillna(0).astype(int)

        # observed breadth from primary-niche species counts (one niche per species)
        prim_ct = (sp.groupby([rank, "niche_primary"], observed=False).size()
                   .unstack(fill_value=0).reindex(columns=niches, fill_value=0))
        prim_ct = prim_ct.reindex(sw.index, fill_value=0)
        obs_bstd = breadth_std(prim_ct.to_numpy(float), K)

        # permutation null on species' primary niche labels
        tax_codes = sp[rank].astype("category")
        codes = tax_codes.cat.codes.to_numpy()
        n_tax = len(tax_codes.cat.categories)
        niche_codes = sp["niche_primary"].cat.codes.to_numpy()
        # map taxon category order to sw.index order
        cat_order = list(tax_codes.cat.categories)
        order_idx = [cat_order.index(t) for t in sw.index]

        ge = np.zeros(n_tax)  # count of null breadth <= observed
        # observed per taxon-category (aligned to cat_order)
        obs_by_cat = pd.Series(obs_bstd, index=sw.index).reindex(cat_order).to_numpy()
        valid = (n_species.reindex(cat_order).fillna(0).to_numpy() >= min_sp)
        rng = np.random.default_rng(cfg["seed"])
        for _ in range(N):
            perm = rng.permutation(niche_codes)
            M = np.zeros((n_tax, K))
            np.add.at(M, (codes, perm), 1)
            bstd = breadth_std(M, K)
            ge += (bstd <= obs_by_cat + 1e-12)
        p = (1.0 + ge) / (N + 1.0)
        p = pd.Series(p, index=cat_order).reindex(sw.index).to_numpy()
        p[~pd.Series(valid, index=cat_order).reindex(sw.index).to_numpy()] = np.nan
        q = np.full(len(p), np.nan)
        m = ~np.isnan(p)
        if m.sum():
            q[m] = multipletests(p[m], method="fdr_bh")[1]

        dominant = prim_ct.idxmax(axis=1)
        specificity = (prim_ct.max(axis=1) / prim_ct.sum(axis=1).clip(lower=1))
        out = pd.DataFrame({
            "taxon": sw.index, "rank": rank.replace("gtdb_", ""),
            "n_species": n_species.values,
            **{f"species_{n}": sw[n].values for n in niches},
            **{f"genomes_{n}": gw.reindex(sw.index, fill_value=0)[n].values for n in niches},
            "levins_B_std": obs_bstd,
            "specificity_max_prop": specificity.values,
            "dominant_niche": dominant.values,
            "specialist": obs_bstd < cutoff,
            "p_specificity": p, "q_specificity": q,
        })
        out.sort_values(["specialist", "n_species"], ascending=[False, False]) \
           .to_csv(f"{args.out_dir}/taxon_specificity_{rank.replace('gtdb_','')}.tsv",
                   sep="\t", index=False)

        # observed vs null specialist fraction at this rank
        valid_mask = out["n_species"] >= min_sp
        obs_spec_frac = out.loc[valid_mask, "specialist"].mean() if valid_mask.any() else np.nan
        sig = (out["q_specificity"] < 0.05)
        for n in niches:
            cnt = int(((out["dominant_niche"] == n) & sig).sum())
            summary.append({"rank": rank.replace("gtdb_", ""), "niche": n,
                            "n_significant_specific_taxa": cnt})
        log.info("%s: %d taxa (>=%d sp), specialist frac=%.2f, %d FDR-sig specific",
                 rank, int(valid_mask.sum()), min_sp,
                 obs_spec_frac if obs_spec_frac == obs_spec_frac else float("nan"),
                 int(sig.sum()))

    pd.DataFrame(summary).to_csv(f"{args.out_dir}/specialist_summary.tsv",
                                 sep="\t", index=False)


if __name__ == "__main__":
    main()
