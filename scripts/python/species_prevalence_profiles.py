#!/usr/bin/env python3
"""
species_prevalence_profiles.py
------------------------------
Collapse per-genome feature tables to SPECIES-level profiles. This is the step
that removes the strains-per-species confounder from the functional data: a
species with 9,606 genomes and a singleton species each contribute exactly one
profile.

For each species s and feature f:
  prevalence(s,f) = (# conspecific genomes carrying f) / (# conspecific genomes)
  present(s,f)    = prevalence(s,f) >= species.presence_threshold
  mean_copies(s,f)= mean copy number of f across conspecific genomes (incl. 0s)

Source priority is applied here:
  * KO       : KOfam primary, eggNOG cross-check -> per-genome KO concordance.
  * CAZyme   : dbCAN primary, eggNOG cross-check -> per-genome CAZy concordance.
Concordance (Jaccard of feature sets per genome) is written out so the choice
of primary source is auditable rather than asserted.

Aggregation is streamed (map-reduce over Arrow record batches) so the billion-
row KO table never has to be held in memory at once.

Outputs (parquet, long format):
  prevalence_<layer>.parquet   species, feature, prevalence, present, mean_copies, n_genomes
  ko_source_concordance.tsv    genome-level Jaccard(KOfam, eggNOG)
  cazyme_source_concordance.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from hgn_utils import load_config, get_logger

log = get_logger("profiles")


def stream_species_counts(parquet_path, genome2species, layer_filter=None,
                          batch_size=2_000_000):
    """Map-reduce: return DataFrame species,feature -> (n_present, sum_copies)."""
    dataset = ds.dataset(parquet_path, format="parquet")
    partials = []

    def consolidate(parts):
        big = pd.concat(parts, ignore_index=True)
        return [big.groupby(["species", "feature"], observed=True, as_index=False)
                   .agg(n_present=("n_present", "sum"),
                        sum_copies=("sum_copies", "sum"))]

    for batch in dataset.to_batches(batch_size=batch_size):
        b = batch.to_pandas()
        if layer_filter is not None:
            b = b[b["layer"] == layer_filter].copy()
            if b.empty:
                continue
        else:
            b = b.copy()
        b["species"] = b["genome"].map(genome2species)
        b = b.dropna(subset=["species"])
        b["n_present"] = 1                      # one row == one genome carrying feature
        b["sum_copies"] = b["count"].astype("int64")
        part = (b.groupby(["species", "feature"], observed=True, as_index=False)
                  .agg(n_present=("n_present", "sum"),
                       sum_copies=("sum_copies", "sum")))
        partials.append(part)
        if len(partials) >= 50:
            partials = consolidate(partials)

    if not partials:
        return pd.DataFrame(columns=["species", "feature", "n_present", "sum_copies"])
    return consolidate(partials)[0]


def build_layer(counts, species_ngenomes, threshold):
    counts = counts.copy()
    counts["n_genomes"] = counts["species"].map(species_ngenomes)
    counts["prevalence"] = counts["n_present"] / counts["n_genomes"]
    counts["mean_copies"] = counts["sum_copies"] / counts["n_genomes"]
    counts["present"] = (counts["prevalence"] >= threshold).astype(int)
    return counts[["species", "feature", "prevalence", "present",
                   "mean_copies", "n_genomes"]]


def genome_jaccard(parquet_a, parquet_b, layer, genomes_subset=None):
    """Per-genome Jaccard of feature sets between two sources for one layer."""
    def feats(path):
        d = ds.dataset(path, format="parquet").to_table(
            filter=ds.field("layer") == layer,
            columns=["genome", "feature"]).to_pandas()
        return d.groupby("genome")["feature"].apply(set)
    A, B = feats(parquet_a), feats(parquet_b)
    common = set(A.index) & set(B.index)
    if genomes_subset is not None:
        common &= set(genomes_subset)
    rows = []
    for g in common:
        a, b = A[g], B[g]
        u = len(a | b)
        rows.append({"genome": g, "jaccard": (len(a & b) / u) if u else np.nan,
                     "n_primary": len(a), "n_crosscheck": len(b)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--kofam", required=True)
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--dbcan", required=True)
    ap.add_argument("--antismash", required=True)
    ap.add_argument("--amrfinder", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    thr = cfg["species"]["presence_threshold"]
    os.makedirs(args.out_dir, exist_ok=True)

    samples = pd.read_parquet(args.samples)[["genome", "species"]]
    g2s = dict(zip(samples["genome"], samples["species"]))
    species_ngenomes = samples.groupby("species")["genome"].nunique().to_dict()

    # layer -> (parquet path, layer tag within file)
    plan = {
        "ko":     (args.kofam, "ko"),       # KOfam primary
        "pfam":   (args.eggnog, "pfam"),
        "cog":    (args.eggnog, "cog"),
        "ec":     (args.eggnog, "ec"),
        "cazyme": (args.dbcan, "cazyme"),   # dbCAN primary
        "bgc":    (args.antismash, "bgc"),
        "amr":    (args.amrfinder, "amr"),
    }
    for out_layer, (path, tag) in plan.items():
        log.info("Building species prevalence for layer '%s' from %s", out_layer, path)
        counts = stream_species_counts(path, g2s, layer_filter=tag)
        prof = build_layer(counts, species_ngenomes, thr)
        prof.to_parquet(f"{args.out_dir}/prevalence_{out_layer}.parquet", index=False)
        log.info("  layer %s: %d species-feature rows, %d features",
                 out_layer, len(prof), prof["feature"].nunique())

    # --- source concordance (auditability) ------------------------------------
    log.info("Computing KO source concordance (KOfam vs eggNOG)")
    genome_jaccard(args.kofam, args.eggnog, "ko").to_csv(
        f"{args.out_dir}/ko_source_concordance.tsv", sep="\t", index=False)
    log.info("Computing CAZyme source concordance (dbCAN vs eggNOG)")
    genome_jaccard(args.dbcan, args.eggnog, "cazyme").to_csv(
        f"{args.out_dir}/cazyme_source_concordance.tsv", sep="\t", index=False)
    log.info("Done.")


if __name__ == "__main__":
    main()
