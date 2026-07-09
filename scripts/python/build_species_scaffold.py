#!/usr/bin/env python3
"""
build_species_scaffold.py
-------------------------
Map each species to the GTDB reference phylogeny and emit a pruned species-level
tree for the comparative analyses. Using the GTDB reference tree as the scaffold
(rather than building a tree de novo) makes the phylogeny fully reproducible and
independent of this study's genomes.

Logic:
  * GTDB reference trees (bac120, ar53) have one tip per species: the species
    representative genome accession.
  * The GTDB taxonomy files map every accession to its species; filtering to
    accessions that are tips gives a species -> representative-accession map.
  * Each dataset species (gtdb_species string) is matched to that map.
  * Bacterial and archaeal trees are pruned separately (config separate_domains)
    because GTDB provides no single cross-domain tree; downstream phylogenetic
    tests run per domain, with bacteria as the primary analysis (>98% of genomes).

Species that are GTDB placeholder/de-novo clusters (no named GTDB tip) are
"unplaced". Per config phylogeny.unplaced_species they are either dropped from
phylogenetic tests (default, honest about coverage) or grafted with GTDB-Tk in a
separate rule. Coverage is always reported.

Outputs:
  species_tree.nwk            bacterial species tree, tips = accession
  species_tree_archaea.nwk    archaeal species tree
  tip_map.tsv                 tip_label(accession), species, domain
  scaffold_coverage.json      placed vs unplaced counts, per domain & niche
"""
import argparse
import json

import pandas as pd
from ete3 import Tree

from hgn_utils import load_config, get_logger

log = get_logger("scaffold")


def species_to_rep(taxonomy_path, tip_set):
    """species -> representative accession, for accessions that are tree tips."""
    s2a = {}
    with open(taxonomy_path) as fh:
        for line in fh:
            acc, tax = line.rstrip("\n").split("\t")[:2]
            if acc not in tip_set:
                continue
            sp = [t for t in tax.split(";") if t.startswith("s__")]
            if sp:
                s2a.setdefault(sp[0], acc)
    return s2a


def prune_relabel(tree_path, keep_accessions):
    t = Tree(tree_path, format=1, quoted_node_names=True)
    tips = set(t.get_leaf_names())
    keep = [a for a in keep_accessions if a in tips]
    if not keep:
        return None, []
    t.prune(keep, preserve_branch_length=True)
    return t, keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ph = cfg["phylogeny"]
    import os
    os.makedirs(args.out_dir, exist_ok=True)

    sp = pd.read_csv(args.species_table, sep="\t")
    bac = sp[sp["domain"] == "d__Bacteria"]["species"].unique()
    arc = sp[sp["domain"] == "d__Archaea"]["species"].unique()

    tip_rows, coverage = [], {}
    for dom, tree_path, tax_path, sp_list, out_name in [
        ("Bacteria", ph["bac120_tree"], ph["bac120_taxonomy"], bac, "species_tree.nwk"),
        ("Archaea", ph["ar53_tree"], ph["ar53_taxonomy"], arc, "species_tree_archaea.nwk"),
    ]:
        log.info("Loading %s reference tree: %s", dom, tree_path)
        t = Tree(tree_path, format=1, quoted_node_names=True)
        tip_set = set(t.get_leaf_names())
        s2a = species_to_rep(tax_path, tip_set)
        keep_acc = [s2a[s] for s in sp_list if s in s2a]
        placed = {s: s2a[s] for s in sp_list if s in s2a}
        log.info("%s: %d/%d species placed on scaffold", dom, len(placed), len(sp_list))

        if keep_acc:
            t.prune(keep_acc, preserve_branch_length=True)
            t.write(outfile=f"{args.out_dir}/{out_name}", format=1)
        for s, acc in placed.items():
            tip_rows.append({"tip_label": acc, "species": s, "domain": dom})
        coverage[dom] = {"n_species": int(len(sp_list)),
                         "n_placed": int(len(placed)),
                         "n_unplaced": int(len(sp_list) - len(placed)),
                         "fraction_placed": round(len(placed) / max(len(sp_list), 1), 4)}

    pd.DataFrame(tip_rows).to_csv(f"{args.out_dir}/tip_map.tsv", sep="\t", index=False)

    placed_species = {r["species"] for r in tip_rows}
    cov_by_niche = {}
    for niche in cfg["inputs"]["niche_levels"]:
        ss = set(sp.loc[sp["niche_primary"] == niche, "species"])
        cov_by_niche[niche] = {
            "n_species": len(ss),
            "n_placed": len(ss & placed_species),
            "fraction_placed": round(len(ss & placed_species) / max(len(ss), 1), 4)}
    coverage["by_primary_niche"] = cov_by_niche
    coverage["unplaced_policy"] = ph["unplaced_species"]
    json.dump(coverage, open(f"{args.out_dir}/scaffold_coverage.json", "w"), indent=2)
    log.info("Scaffold coverage by niche: %s",
             {k: v["fraction_placed"] for k, v in cov_by_niche.items()})


if __name__ == "__main__":
    main()
