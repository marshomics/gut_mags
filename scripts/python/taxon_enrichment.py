#!/usr/bin/env python3
"""
taxon_enrichment.py
-------------------
Which taxa are over- or under-represented in each niche, at every rank, on a
species basis (each species counted once). For niche n and taxon t the 2x2 table
is (species in t vs not) x (species with primary niche n vs the rest); Fisher's
exact test gives an odds ratio and p, FDR-corrected across all taxon x niche
tests within a rank. A rank-level G-test of independence with standardised
Pearson residuals is also reported, so both the global pattern and the specific
taxa are covered.

Species-weighting removes the strain bias; using each species' primary niche
keeps every species in exactly one cell so the contingency tables are clean.

Outputs:
  enrichment_<rank>.tsv      taxon, niche, n_taxon_in_niche, log2_or, p, q, direction
  gtest_<rank>.tsv           rank-level G statistic, df, p, and residual matrix
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from statsmodels.stats.multitest import multipletests

from hgn_utils import load_config, get_logger


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    log = get_logger("enrich")
    tcfg = cfg["taxonomy"]
    niches = cfg["inputs"]["niche_levels"]
    min_sp = tcfg["enrichment"]["min_taxon_species"]
    use_g = tcfg["enrichment"]["test"] == "gtest"

    import os
    os.makedirs(args.out_dir, exist_ok=True)

    samples = pd.read_parquet(args.samples)
    sptab = pd.read_csv(args.species_table, sep="\t").set_index("species")
    sp_tax = samples.groupby("species")[tcfg["ranks"]].first()
    sp = sp_tax.join(sptab["niche_primary"]).dropna(subset=["niche_primary"])
    sp = sp[sp["niche_primary"].isin(niches)]
    total = len(sp)

    for rank in tcfg["ranks"]:
        if rank == "gtdb_species":
            continue
        ct = (sp.groupby([rank, "niche_primary"]).size()
              .unstack(fill_value=0).reindex(columns=niches, fill_value=0))
        ct = ct[ct.sum(axis=1) >= min_sp]
        if ct.empty:
            continue
        niche_tot = sp["niche_primary"].value_counts().reindex(niches).fillna(0)

        rows = []
        for taxon, r in ct.iterrows():
            t_tot = r.sum()
            for n in niches:
                a = int(r[n])                       # taxon & niche
                b = int(t_tot - a)                  # taxon & not niche
                c = int(niche_tot[n] - a)           # not taxon & niche
                d = int(total - t_tot - c)          # not taxon & not niche
                if use_g:
                    tbl = np.array([[a, b], [c, d]]) + 0.5
                    _, p, _, _ = chi2_contingency(tbl, lambda_="log-likelihood")
                    orr = (tbl[0, 0] * tbl[1, 1]) / (tbl[0, 1] * tbl[1, 0])
                else:
                    orr, p = fisher_exact([[a, b], [c, d]])
                    orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
                rows.append({"taxon": taxon, "rank": rank.replace("gtdb_", ""),
                             "niche": n, "n_taxon_in_niche": a,
                             "n_taxon_total": int(t_tot),
                             "log2_or": float(np.log2(orr)), "p": p})
        df = pd.DataFrame(rows)
        df["q"] = multipletests(df["p"], method="fdr_bh")[1]
        df["direction"] = np.where(df["log2_or"] > 0, "enriched", "depleted")
        df.sort_values(["niche", "q", "log2_or"], ascending=[True, True, False]) \
          .to_csv(f"{args.out_dir}/enrichment_{rank.replace('gtdb_','')}.tsv",
                  sep="\t", index=False)

        # rank-level G-test of independence + standardised residuals
        g, p, dof, exp = chi2_contingency(ct.to_numpy() + 0.5, lambda_="log-likelihood")
        resid = (ct.to_numpy() - exp) / np.sqrt(exp)
        res_df = pd.DataFrame(resid, index=ct.index, columns=[f"resid_{n}" for n in niches])
        res_df.insert(0, "G_statistic", g); res_df.insert(1, "df", dof); res_df.insert(2, "p", p)
        res_df.reset_index().to_csv(
            f"{args.out_dir}/gtest_{rank.replace('gtdb_','')}.tsv", sep="\t", index=False)
        log.info("%s: %d taxa tested, %d enriched calls (q<0.05)",
                 rank, len(ct), int((df["q"] < 0.05).sum()))


if __name__ == "__main__":
    main()
