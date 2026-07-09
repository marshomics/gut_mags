#!/usr/bin/env python3
"""
community_ingest.py
-------------------
Turn one per-genome metabolic/trait prediction table (input/) into a species-level
prevalence profile in the pipeline's standard schema, so it plugs straight into
the differential + consensus stack and the community analyses. This is where the
strains-per-species confounder is removed for these data: each genome is a strain,
and a capability is called present in a species when at least `presence_threshold`
of that species' strains have it.

Types:
  auxotrophy  amino-acid biosynthesis table (1 = can synthesise). INVERTED so the
              feature = auxotrophy (present = the species cannot synthesise it,
              i.e. depends on host/diet/community).
  carbon      carbon-source utilisation (1 = can use).
  trait       predicted phenotypes (anaerobe, spore formation, bile, ...).
  modulep     predicted KEGG module presence (wide).

Only genomes present in the QC-passed sample sheet are used (so species/niche are
defined consistently with the rest of the pipeline). The genome-id column is the
first column whatever its name (genome / file / unnamed).

Output: prevalence_<type>.parquet (species, feature, prevalence, present,
mean_copies, n_genomes)  written into 03_profiles.
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("comm-ingest")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--type", required=True,
                    choices=["auxotrophy", "carbon", "trait", "modulep"])
    ap.add_argument("--input", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cc = cfg["community"]
    thr = cc["presence_threshold"]

    g2s = (pd.read_parquet(args.samples)[["genome", "species"]]
           .assign(genome=lambda d: d["genome"].astype(str)))
    gmap = dict(zip(g2s["genome"], g2s["species"]))

    df = pd.read_csv(args.input, sep="\t", dtype=str, keep_default_na=False)
    gid = df.columns[0]
    df = df.rename(columns={gid: "genome"})
    df["genome"] = df["genome"].astype(str)
    df = df[df["genome"].isin(gmap)]
    df["species"] = df["genome"].map(gmap)

    feature_cols = [c for c in df.columns if c not in ("genome", "species")]
    if args.type == "auxotrophy":
        feature_cols = [c for c in feature_cols if c in cc["amino_acids"]]
    # to numeric 0/1
    X = df[feature_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    X["species"] = df["species"].values

    # per-species capability prevalence (fraction of strains with capability=1)
    grp = X.groupby("species")
    cap = grp.mean(numeric_only=True)
    ngen = grp.size().rename("n_genomes")

    long = cap.reset_index().melt(id_vars="species", var_name="feature",
                                  value_name="cap_prev")
    long = long.merge(ngen, on="species", how="left")
    if args.type == "auxotrophy":
        # feature = auxotrophy: present when the species mostly CANNOT synthesise
        long["prevalence"] = (1.0 - long["cap_prev"]).round(4)
        long["present"] = (long["cap_prev"] < thr).astype(int)
    else:
        long["prevalence"] = long["cap_prev"].round(4)
        long["present"] = (long["cap_prev"] >= thr).astype(int)
    long["mean_copies"] = long["prevalence"]
    out = long[["species", "feature", "prevalence", "present", "mean_copies", "n_genomes"]]
    out.to_parquet(args.out, index=False)
    log.info("[%s] %d genomes -> %d species x %d features; %d present calls",
             args.type, len(df), out["species"].nunique(), len(feature_cols),
             int(out["present"].sum()))


if __name__ == "__main__":
    main()
