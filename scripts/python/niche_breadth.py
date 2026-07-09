#!/usr/bin/env python3
"""
niche_breadth.py
----------------
Quantify how niche-specific the species are, and prepare the overlap data for
the UpSet figure. Operates on species units (one species = one observation), so
the strains-per-species bias does not inflate any niche.

Metrics:
  * Levins' niche breadth B = 1 / sum(p_i^2), where p_i is the fraction of a
    species' genomes found in niche i. B near 1 = specialist, near 3 = even
    across all three niches. Standardised B_A = (B-1)/(N-1) in [0,1].
  * Specialist / generalist counts per niche.
  * Occupancy set sizes for UpSet (named vs placeholder species reported
    separately, because placeholder clusters are not cross-database comparable).

Outputs:
  species_breadth.tsv      species, B, B_std, occupancy_pattern, is_placeholder
  upset_sets.tsv           occupancy_pattern, n_species_named, n_species_placeholder
  breadth_summary.json
"""
import argparse
import json

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    log = get_logger("breadth")
    niches = cfg["inputs"]["niche_levels"]
    N = len(niches)

    sp = pd.read_csv(args.species_table, sep="\t")
    ncols = [f"n_{c}" for c in niches]
    counts = sp[ncols].to_numpy(float)
    tot = counts.sum(axis=1, keepdims=True)
    p = np.divide(counts, tot, out=np.zeros_like(counts), where=tot > 0)

    B = 1.0 / np.clip((p ** 2).sum(axis=1), 1e-12, None)
    sp["levins_B"] = B
    sp["levins_B_std"] = (B - 1.0) / (N - 1.0)

    occ = (counts >= cfg["occupancy"]["min_genomes_per_niche"]).astype(int)
    sp["occupancy_pattern"] = ["+".join([niches[i] for i in range(N) if row[i]])
                               for row in occ]

    sp[["species", "levins_B", "levins_B_std", "occupancy_pattern",
        "species_is_placeholder"]].to_csv(
        f"{args.out_prefix}_species_breadth.tsv", sep="\t", index=False)

    upset = (sp.groupby(["occupancy_pattern", "species_is_placeholder"])
               .size().rename("n").reset_index())
    upset = upset.pivot(index="occupancy_pattern", columns="species_is_placeholder",
                        values="n").fillna(0).astype(int)
    upset.columns = ["n_species_named" if c is False else "n_species_placeholder"
                     for c in upset.columns]
    for need in ["n_species_named", "n_species_placeholder"]:
        if need not in upset.columns:
            upset[need] = 0
    upset.reset_index().to_csv(f"{args.out_prefix}_upset_sets.tsv",
                               sep="\t", index=False)

    summary = {
        "n_species": int(len(sp)),
        "specialists": int((sp["levins_B_std"] < 0.25).sum()),
        "generalists_even": int((sp["levins_B_std"] > 0.75).sum()),
        "median_breadth_std": float(sp["levins_B_std"].median()),
        "occupancy_patterns": upset.reset_index().set_index("occupancy_pattern")
                                   .sum(axis=1).to_dict(),
    }
    with open(f"{args.out_prefix}_breadth_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("Niche breadth summarised: %s", summary["occupancy_patterns"])


if __name__ == "__main__":
    main()
