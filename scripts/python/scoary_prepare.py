#!/usr/bin/env python3
"""
scoary_prepare.py
-----------------
Build Scoary2 inputs for one functional layer and one niche contrast from the
species-level material the differential stage already produced. Running Scoary2
on the species-level functional presence matrix (one row per species) plus the
GTDB species tree means the strains-per-species bias is already removed and the
pairwise-comparisons correction has a real tree to work with.

Writes three files with matching isolate (species) names:
  genes.tsv   rows = features, cols = species, 0/1  (Scoary2 gene-count format)
  traits.tsv  rows = species, one column = the contrast (0/1)  (binary format)
  tree.nwk    GTDB species tree pruned to these species, tips relabelled to ids

Only species placed on the GTDB scaffold are kept (the tree is required for the
pairwise test); coverage is reported.
"""
import argparse
import os
import re

import pandas as pd
import dendropy

from hgn_utils import load_config, get_logger

log = get_logger("scoary-prep")


def sid(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--tip-map", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    os.makedirs(args.out_dir, exist_ok=True)

    analysis = pd.read_csv(args.analysis, sep="\t")
    tipmap = pd.read_csv(args.tip_map, sep="\t")          # tip_label(accession), species
    sp2acc = dict(zip(tipmap["species"], tipmap["tip_label"]))

    # keep analysis species that are on the tree
    analysis = analysis[analysis["species"].isin(sp2acc)].copy()
    analysis["sid"] = analysis["species"].map(sid)
    keep_species = analysis["species"].tolist()
    if len(analysis) < 6:
        log.warning("Too few placed species (%d) for %s/%s; writing stubs.",
                    len(analysis), args.layer, args.contrast)

    # genes table: feature x species (0/1)
    pres = pd.read_parquet(args.presence)
    pres = pres[pres["species"].isin(keep_species)]
    wide = (pres.pivot_table(index="feature", columns="species", values="present",
                             fill_value=0).astype(int))
    wide = wide.reindex(columns=keep_species, fill_value=0)
    wide.columns = [sid(c) for c in wide.columns]
    wide.index.name = "Gene"
    wide.to_csv(f"{args.out_dir}/genes.tsv", sep="\t")

    # traits table: species x contrast (0/1)
    traits = analysis[["sid", "group"]].copy()
    traits["group"] = traits["group"].astype(int)
    traits = traits.rename(columns={"sid": "Isolate", "group": args.contrast}).set_index("Isolate")
    traits.to_csv(f"{args.out_dir}/traits.tsv", sep="\t")

    # tree pruned to these species, tips relabelled accession -> sid
    tree = dendropy.Tree.get(path=args.tree, schema="newick",
                             preserve_underscores=True)
    acc_keep = {sp2acc[s] for s in keep_species}
    taxa_keep = [t for t in tree.taxon_namespace if t.label in acc_keep]
    tree.retain_taxa(taxa_keep)
    acc2sid = {sp2acc[s]: sid(s) for s in keep_species}
    for leaf in tree.leaf_node_iter():
        if leaf.taxon and leaf.taxon.label in acc2sid:
            leaf.taxon.label = acc2sid[leaf.taxon.label]
    tree.write(path=f"{args.out_dir}/tree.nwk", schema="newick",
               suppress_rooting=True, unquoted_underscores=True)

    log.info("%s/%s: %d species, %d features -> Scoary2 inputs",
             args.layer, args.contrast, len(keep_species), wide.shape[0])


if __name__ == "__main__":
    main()
