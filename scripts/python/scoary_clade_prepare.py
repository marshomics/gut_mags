#!/usr/bin/env python3
"""
scoary_clade_prepare.py
-----------------------
Build Scoary2 inputs for one clade (genus) from its Panaroo pangenome, for one
contrast. Genes are real ortholog families (Panaroo gene_presence_absence.Rtab),
samples are the clade's species representatives, the trait is niche, and the tree
is the GTDB species subtree relabelled to the same species ids.

Outputs: genes.tsv, traits.tsv, tree.nwk
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import dendropy

from hgn_utils import load_config, get_logger

log = get_logger("scoary-clade-prep")


def sid(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def group_for(contrast, niche, niches):
    if contrast == "host_vs_free":
        if niche in ("human", "animal"):
            return 1
        if niche == "free":
            return 0
        return np.nan
    if contrast.endswith("_vs_rest"):
        x = contrast[:-len("_vs_rest")]
        return 1 if niche == x else (0 if niche in niches else np.nan)
    if "_vs_" in contrast:
        x, y = contrast.split("_vs_")
        return 1 if niche == x else (0 if niche == y else np.nan)
    return np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-genomes", required=True)
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--tip-map", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    niches = cfg["inputs"]["niche_levels"]
    os.makedirs(args.out_dir, exist_ok=True)

    sg = pd.read_csv(args.species_genomes, sep="\t")
    sg["group"] = sg["niche"].map(lambda n: group_for(args.contrast, n, niches))
    sg = sg.dropna(subset=["group"])
    sg["sid"] = sg["species"].map(sid)
    g2sid = dict(zip(sg["genome"].astype(str), sg["sid"]))

    # genes from Panaroo Rtab (genefamily x genome) -> genefamily x sid
    pa = pd.read_csv(args.rtab, sep="\t", index_col=0)
    pa.columns = [c.split("/")[-1].replace(".gff", "").replace(".fna", "") for c in pa.columns]
    cols = [c for c in pa.columns if c in g2sid]
    pa = pa[cols]
    pa.columns = [g2sid[c] for c in cols]
    pa = pa.loc[(pa.sum(axis=1) > 0) & (pa.sum(axis=1) < pa.shape[1])]   # accessory only
    pa.index.name = "Gene"
    pa.astype(int).to_csv(f"{args.out_dir}/genes.tsv", sep="\t")

    traits = sg[["sid", "group"]].drop_duplicates("sid")
    traits["group"] = traits["group"].astype(int)
    traits.rename(columns={"sid": "Isolate", "group": args.contrast}).set_index(
        "Isolate").to_csv(f"{args.out_dir}/traits.tsv", sep="\t")

    # tree: prune GTDB scaffold to these species, relabel accession->sid
    tipmap = pd.read_csv(args.tip_map, sep="\t")
    sp2acc = dict(zip(tipmap["species"], tipmap["tip_label"]))
    keep = sg[sg["species"].isin(sp2acc)]
    acc2sid = {sp2acc[s]: sid(s) for s in keep["species"]}
    if len(acc2sid) >= 4:
        tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True)
        tree.retain_taxa([t for t in tree.taxon_namespace if t.label in acc2sid])
        for leaf in tree.leaf_node_iter():
            if leaf.taxon and leaf.taxon.label in acc2sid:
                leaf.taxon.label = acc2sid[leaf.taxon.label]
        tree.write(path=f"{args.out_dir}/tree.nwk", schema="newick",
                   suppress_rooting=True, unquoted_underscores=True)
    else:
        open(f"{args.out_dir}/tree.nwk", "w").write("")   # Scoary2 will build one
    log.info("clade %s: %d species, %d accessory families", args.contrast,
             len(traits), pa.shape[0])


if __name__ == "__main__":
    main()
