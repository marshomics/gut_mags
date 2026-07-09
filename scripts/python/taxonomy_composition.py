#!/usr/bin/env python3
"""
taxonomy_composition.py
-----------------------
Taxonomic composition of each niche, computed TWO ways on purpose:

  * genome-weighted : proportion of genomes in each taxon (the naive view, which
                      is dominated by over-sequenced human commensals).
  * species-weighted: proportion of *species* in each taxon (each species counts
                      once). This is the view used for interpretation because it
                      removes the strains-per-species sampling bias.

Reporting both, side by side, makes the confounder visible rather than hidden,
and lets a reviewer see exactly how much the bias matters.

Outputs (long-format, ready for plotting):
  composition_<rank>.tsv   columns: niche, taxon, weighting, proportion, n
"""
import argparse

import pandas as pd

from hgn_utils import load_config, get_logger


def comp(df, rank, group, weighting):
    """Proportion of `group` within each niche at taxonomic `rank`."""
    g = df.groupby(["niche", rank]).size().rename("n").reset_index()
    tot = g.groupby("niche")["n"].transform("sum")
    g["proportion"] = g["n"] / tot
    g["weighting"] = weighting
    g = g.rename(columns={rank: "taxon"})
    return g[["niche", "taxon", "weighting", "proportion", "n"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    log = get_logger("composition")
    ranks = cfg["inputs"]["columns"]["gtdb_ranks"]

    genomes = pd.read_parquet(args.samples)
    species = pd.read_csv(args.species_table, sep="\t")

    # species-weighted frame: one row per (species, niche occupied)
    niches = cfg["inputs"]["niche_levels"]
    occ_long = []
    for n in niches:
        sub = species[species[f"n_{n}"] > 0].copy()
        sub["niche"] = n
        occ_long.append(sub)
    occ_long = pd.concat(occ_long, ignore_index=True)

    import os
    os.makedirs(args.out_dir, exist_ok=True)
    for rank in ranks[:-1]:   # all ranks except species itself
        gw = comp(genomes, rank, "genome", "genome_weighted")
        sw = comp(occ_long, rank, rank, "species_weighted")
        out = pd.concat([gw, sw], ignore_index=True)
        short = rank.replace("gtdb_", "")
        out.to_csv(f"{args.out_dir}/composition_{short}.tsv", sep="\t", index=False)
        log.info("Wrote composition for rank %s", short)


if __name__ == "__main__":
    main()
