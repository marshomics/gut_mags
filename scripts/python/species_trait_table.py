#!/usr/bin/env python3
"""
species_trait_table.py
----------------------
Per-species continuous traits for the phylogenetic comparative methods:
functional richness (number of distinct features) and load (summed mean copy
number) for each layer, plus genome-level traits. One row per species, niche
labelled, so PGLS / signal tests treat each species once.

Output: species_traits.tsv
  species, domain, genus, family, niche_primary, specialist_niche,
  completeness_mean, gc_mean, genome_size_mean, cds_number_mean,
  <layer>_richness, <layer>_load   for layer in ko,pfam,cog,cazyme,bgc,amr
"""
import argparse
import glob
import os

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("traits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    samples = pd.read_parquet(args.samples)
    sp = pd.read_csv(args.species_table, sep="\t")

    gstats = (samples.groupby("species")
              .agg(completeness_mean=("completeness", "mean"),
                   gc_mean=("gc", "mean"),
                   genome_size_mean=("genome_size", "mean"),
                   cds_number_mean=("cds_number", "mean"))
              .reset_index())

    traits = sp[["species", "domain", "gtdb_genus", "gtdb_family",
                 "niche_primary", "specialist_niche"]].rename(
        columns={"gtdb_genus": "genus", "gtdb_family": "family"})
    traits = traits.merge(gstats, on="species", how="left")

    for path in sorted(glob.glob(f"{args.profiles_dir}/prevalence_*.parquet")):
        layer = os.path.basename(path).replace("prevalence_", "").replace(".parquet", "")
        p = pd.read_parquet(path)
        rich = (p[p["present"] == 1].groupby("species").size()
                .rename(f"{layer}_richness"))
        load = (p.groupby("species")["mean_copies"].sum()
                .rename(f"{layer}_load"))
        traits = traits.merge(rich, on="species", how="left") \
                       .merge(load, on="species", how="left")
        log.info("Added traits for layer %s", layer)

    for c in traits.columns:
        if c.endswith("_richness") or c.endswith("_load"):
            traits[c] = traits[c].fillna(0)
    traits.to_csv(args.out, sep="\t", index=False)
    log.info("Wrote %d species trait rows", len(traits))


if __name__ == "__main__":
    main()
