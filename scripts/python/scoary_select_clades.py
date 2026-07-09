#!/usr/bin/env python3
"""
scoary_select_clades.py
-----------------------
Pick clades (genera by default) for clade-stratified Scoary2, where ortholog
families are comparable so the pan-GWAS tests real gene families rather than
functional categories. A genus qualifies when it has enough species in each of
at least two niches. One representative genome per species is used (species-level
unit, tractable Panaroo).

Outputs a checkpoint directory:
  manifest.tsv                       qualifying clades + per-niche species counts
  clades/<clade_id>/species_genomes.tsv   species, genome(rep), niche
"""
import argparse
import os
import re

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("scoary-clades")


def cid(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--representatives", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cc = cfg["scoary"]["clade"]
    niches = cfg["inputs"]["niche_levels"]
    rank = cc["rank"]
    assign = cfg["scoary"]["niche_assignment"]
    os.makedirs(f"{args.out_dir}/clades", exist_ok=True)

    sp = pd.read_csv(args.species_table, sep="\t")
    reps = pd.read_csv(args.representatives, sep="\t")
    sp = sp.merge(reps[["species", reps.columns[0] if "genome" not in reps.columns else "genome"]]
                  .rename(columns={reps.columns[0]: "genome"}) if "genome" not in reps.columns
                  else reps[["species", "genome"]], on="species", how="left")
    sp["niche"] = sp["specialist_niche"] if assign == "specialist" else sp["niche_primary"]
    sp = sp[sp["niche"].isin(niches) & sp["genome"].notna()]

    rows = []
    for clade, d in sp.groupby(rank):
        per = {n: (d["niche"] == n).sum() for n in niches}
        ok = [n for n in niches if per[n] >= cc["min_species_per_niche"]]
        if len(ok) < cc["min_niches"]:
            continue
        clid = cid(clade)
        sub = d[d["niche"].isin(ok)][["species", "genome", "niche"]]
        os.makedirs(f"{args.out_dir}/clades/{clid}", exist_ok=True)
        sub.to_csv(f"{args.out_dir}/clades/{clid}/species_genomes.tsv", sep="\t", index=False)
        rows.append({"clade_id": clid, "clade": clade,
                     **{f"n_{n}": per[n] for n in niches},
                     "qual_niches": ",".join(ok), "n_species": int(len(sub))})

    man = pd.DataFrame(rows).sort_values("n_species", ascending=False)
    if cc["top_clades"] and len(man) > cc["top_clades"]:
        man = man.head(cc["top_clades"])
        keep = set(man["clade_id"])
        import shutil
        for d in os.listdir(f"{args.out_dir}/clades"):
            if d not in keep:
                shutil.rmtree(f"{args.out_dir}/clades/{d}", ignore_errors=True)
    man.to_csv(f"{args.out_dir}/manifest.tsv", sep="\t", index=False)
    log.info("Clade-stratified Scoary2: %d %s clades qualify", len(man), rank)


if __name__ == "__main__":
    main()
