#!/usr/bin/env python3
"""
transition_select.py
--------------------
Pick the species worth a within-species niche-transition analysis and assemble
their inputs. A species qualifies when it has at least `min_strains_per_niche`
near-complete genomes (completeness >= completeness_min, contamination <=
contamination_max) in at least `min_niches_occupied` niches. Requiring adequate
sampling in EACH niche is the first defensible control: diversity and the SFS
are sample-size sensitive, so a niche with a handful of genomes cannot anchor a
"recent acquisition" claim.

For each qualifying species an outgroup is chosen for rooting (no root => no
directionality): a congeneric species present in the dataset with a high-quality
genome (the GTDB scaffold places congenerics as the nearest sisters), falling
back to the same family, then to midpoint rooting (flagged) if nothing closer
exists.

Outputs (a checkpoint directory):
  manifest.tsv                            qualifying species + counts + outgroup
  species/<species_id>/candidate_genomes.tsv   genome, niche, role(focal/outgroup)
"""
import argparse
import os
import re

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("trans-select")


def species_id(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def best_outgroup(genus, family, focal, hq_by_species, rep_quality):
    """Return (og_species, og_genome, level) using nearest taxonomy with an HQ genome."""
    for level, group in (("genus", genus), ("family", family)):
        cands = [s for s in group if s != focal and s in hq_by_species]
        if cands:
            # pick the species whose best HQ genome has the highest quality score
            best_s = max(cands, key=lambda s: rep_quality.get(s, (None, -1))[1])
            og_genome = rep_quality[best_s][0]
            return best_s, og_genome, level
    return None, None, "midpoint"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    tc = cfg["transition"]
    niches = cfg["inputs"]["niche_levels"]
    os.makedirs(f"{args.out_dir}/species", exist_ok=True)

    s = pd.read_parquet(args.samples)
    hq = s[(s["completeness"] >= tc["completeness_min"]) &
           (s["contamination"] <= tc["contamination_max"])].copy()

    # per-species best HQ genome (for outgroup quality) + genus/family map
    rep_quality, genus_of, family_of = {}, {}, {}
    for sp, d in hq.groupby("species"):
        top = d.loc[d["quality_score"].idxmax()]
        rep_quality[sp] = (top["genome"], float(top["quality_score"]))
        genus_of[sp] = d["gtdb_genus"].iloc[0]
        family_of[sp] = d["gtdb_family"].iloc[0]
    hq_by_species = set(rep_quality)
    by_genus = hq.groupby("gtdb_genus")["species"].apply(lambda x: sorted(set(x))).to_dict()
    by_family = hq.groupby("gtdb_family")["species"].apply(lambda x: sorted(set(x))).to_dict()

    # per-species HQ genomes per niche
    rows = []
    for sp, d in hq.groupby("species"):
        per = {n: d[d["niche"] == n] for n in niches}
        occ = {n: len(per[n]) for n in niches}
        qual_niches = [n for n in niches if occ[n] >= tc["min_strains_per_niche"]]
        if len(qual_niches) < tc["min_niches_occupied"]:
            continue
        sid = species_id(sp)
        og_sp, og_gen, level = best_outgroup(
            by_genus.get(genus_of[sp], []), by_family.get(family_of[sp], []),
            sp, hq_by_species, rep_quality)
        # write candidate genome list (focal niches that qualify + outgroup)
        recs = []
        for n in qual_niches:
            for g in per[n]["genome"]:
                recs.append({"genome": g, "niche": n, "role": "focal"})
        if og_gen is not None:
            recs.append({"genome": og_gen, "niche": "outgroup", "role": "outgroup"})
        os.makedirs(f"{args.out_dir}/species/{sid}", exist_ok=True)
        pd.DataFrame(recs).to_csv(
            f"{args.out_dir}/species/{sid}/candidate_genomes.tsv", sep="\t", index=False)
        rows.append({"species_id": sid, "species": sp,
                     "genus": genus_of[sp], "family": family_of[sp],
                     **{f"n_{n}": occ[n] for n in niches},
                     "qual_niches": ",".join(qual_niches),
                     "n_total": int(sum(occ[n] for n in qual_niches)),
                     "outgroup_species": og_sp or "", "outgroup_genome": og_gen or "",
                     "outgroup_level": level})

    man = pd.DataFrame(rows).sort_values("n_total", ascending=False)
    if tc["max_species"] and len(man) > tc["max_species"]:
        man = man.head(tc["max_species"])
        keep = set(man["species_id"])
        for d in os.listdir(f"{args.out_dir}/species"):
            if d not in keep:
                import shutil
                shutil.rmtree(f"{args.out_dir}/species/{d}", ignore_errors=True)
    man.to_csv(f"{args.out_dir}/manifest.tsv", sep="\t", index=False)
    log.info("Qualifying species: %d (outgroup level: %s)",
             len(man), man["outgroup_level"].value_counts().to_dict() if len(man) else {})
    if len(man):
        log.info("Top: %s", man.head(5)[["species", "qual_niches", "n_total"]].to_dict("records"))


if __name__ == "__main__":
    main()
