#!/usr/bin/env python3
"""
differential_prep.py
--------------------
Assemble the analysis-ready inputs shared by all three differential-function
methods (phyloglm, CMH, resampling), so they test exactly the same features on
exactly the same species with exactly the same covariates. Consistency across
methods is only meaningful if their inputs are identical; this script guarantees
that.

Produces, for one functional layer and one contrast:
  analysis_species_<contrast>.tsv
      species, group (1=focal/human, 0=comparator), genus, family, domain,
      completeness_mean, log10_genome_size_mean, gc_mean, n_genomes
  presence_<layer>_<contrast>.parquet
      species, feature, present   (restricted to tested species & features)
  tested_features_<layer>_<contrast>.txt
      feature ids passing the per-niche prevalence filter

Confounder handling baked in here:
  * species-level covariates (quality, size, GC) computed as the mean over a
    species' genomes, ready to enter every model;
  * prevalence filter applied PER GROUP (max across groups) so rare-feature
    removal does not favour the larger group;
  * only species with a defined contrast membership are kept (specialists for
    the primary contrast, per config).
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("diffprep")


def positive_label(contrast):
    """The niche/label that is group 1 (the 'enriched in' side) for a contrast."""
    if contrast == "host_vs_free":
        return "host"
    return contrast.split("_vs_")[0]


def build_groups(species_tbl, contrast, all_niches):
    """Map species -> group (1 = positive side, 0 = comparator, NaN = excluded).

    Supported contrasts:
      <niche>_vs_rest     niche vs all other niches
      <niche>_vs_<niche>  the two named niches only
      host_vs_free        host-associated (human+animal) vs free-living
    Niche membership uses the species' contrast_niche (specialist by default).
    """
    s = species_tbl.copy()
    niche = s["contrast_niche"].astype(object)
    if contrast == "western_vs_nonwestern":
        # within-human stratification by lifestyle (population column)
        pop = s.get("population", pd.Series("", index=s.index)).astype(object)
        grp = np.where(pop == "western", 1.0,
              np.where(pop == "non_western", 0.0, np.nan))
        s["group"] = grp
        return s.dropna(subset=["group"])
    if contrast == "host_vs_free":
        grp = np.where(niche.isin(["human", "animal"]), 1.0,
              np.where(niche == "free", 0.0, np.nan))
    elif contrast.endswith("_vs_rest"):
        x = contrast[:-len("_vs_rest")]
        others = [n for n in all_niches if n != x]
        grp = np.where(niche == x, 1.0,
              np.where(niche.isin(others), 0.0, np.nan))
    elif "_vs_" in contrast:
        x, y = contrast.split("_vs_")
        grp = np.where(niche == x, 1.0, np.where(niche == y, 0.0, np.nan))
    else:
        raise ValueError(f"unknown contrast {contrast}")
    s["group"] = grp
    return s.dropna(subset=["group"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--prevalence", required=True, help="prevalence_<layer>.parquet")
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    niches = cfg["inputs"]["niche_levels"]
    min_prev = cfg["functional_annotation"]["min_prevalence_any_niche"]

    samples = pd.read_parquet(args.samples)
    sp = pd.read_csv(args.species_table, sep="\t")

    # species-level covariates (means over conspecific genomes)
    cov = (samples.groupby("species")
                  .agg(completeness_mean=("completeness", "mean"),
                       log10_genome_size_mean=("log10_genome_size", "mean"),
                       gc_mean=("gc", "mean"),
                       n_genomes=("genome", "nunique"))
                  .reset_index())
    sp = sp.merge(cov, on="species", how="left")

    grouped = build_groups(sp, args.contrast, niches)
    keep_cols = ["species", "group", "gtdb_genus", "gtdb_family", "domain",
                 "completeness_mean", "log10_genome_size_mean", "gc_mean", "n_genomes"]
    keep_cols = [c for c in keep_cols if c in grouped.columns]
    analysis = grouped[keep_cols].rename(columns={"gtdb_genus": "genus",
                                                  "gtdb_family": "family"})
    analysis.to_csv(f"{args.out_prefix}_analysis_species.tsv", sep="\t", index=False)
    log.info("Contrast %s: %d focal vs %d comparator species",
             args.contrast, int((analysis['group'] == 1).sum()),
             int((analysis['group'] == 0).sum()))

    # presence matrix restricted to analysis species
    prev = pd.read_parquet(args.prevalence)[["species", "feature", "present"]]
    prev = prev.merge(analysis[["species", "group"]], on="species", how="inner")

    # per-group prevalence of each feature; keep if max across groups >= min_prev
    grp_sizes = analysis.groupby("group")["species"].nunique()
    fp = (prev[prev["present"] == 1]
          .groupby(["feature", "group"])["species"].nunique()
          .unstack(fill_value=0))
    for g in (0.0, 1.0):
        if g in fp.columns and g in grp_sizes.index:
            fp[g] = fp[g] / grp_sizes[g]
    fp["max_prev"] = fp.max(axis=1)
    tested = fp.index[fp["max_prev"] >= min_prev].tolist()

    out = prev[prev["feature"].isin(tested)][["species", "feature", "present"]]
    out.to_parquet(f"{args.out_prefix}_presence.parquet", index=False)
    with open(f"{args.out_prefix}_tested_features.txt", "w") as fh:
        fh.write("\n".join(map(str, tested)) + "\n")
    log.info("Layer %s: %d features pass prevalence>=%.2f in >=1 group",
             args.layer, len(tested), min_prev)


if __name__ == "__main__":
    main()
