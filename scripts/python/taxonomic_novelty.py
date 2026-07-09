#!/usr/bin/env python3
"""
taxonomic_novelty.py
--------------------
How much undescribed diversity sits in each niche? GTDB placeholder labels
(species like 's__Genus [species N]', and bracketed/!-suffixed higher ranks)
mark lineages without a validly published name, i.e. novelty. This quantifies
novelty per niche and per rank, and -- critically -- at matched sampling effort,
because the niche with the most genomes will otherwise look the most novel for
purely sampling reasons.

For each niche:
  * observed fraction of novel (placeholder) species;
  * rarefied novel fraction: subsample each niche to the smallest niche's species
    count, recompute, bootstrap for a CI. This is the comparison that survives
    review.
Per rank: count and fraction of placeholder taxa per niche.

Outputs: novelty_by_niche.tsv, novelty_by_rank.tsv
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import (load_config, get_logger, set_global_seed, derive_seed,
                       is_placeholder_species)

log = get_logger("novelty")


def is_placeholder_label(x: str) -> bool:
    if x is None or x == "":
        return True
    x = str(x)
    # GTDB placeholder cues: brackets, or alphanumeric-only epithet codes
    return ("[" in x) or ("]" in x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-boot", type=int, default=500)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    niches = cfg["inputs"]["niche_levels"]
    ranks = cfg["taxonomy"]["ranks"]

    import os
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    samples = pd.read_parquet(args.samples)
    sptab = pd.read_csv(args.species_table, sep="\t").set_index("species")
    sp = samples.groupby("species")[ranks].first()
    sp["novel"] = sp.index.map(is_placeholder_species)
    sp = sp.join(sptab["niche_primary"]).dropna(subset=["niche_primary"])
    sp = sp[sp["niche_primary"].isin(niches)]

    # --- per niche, observed + rarefied novel species fraction ---
    min_n = sp["niche_primary"].value_counts().min()
    rows = []
    for n in niches:
        d = sp[sp["niche_primary"] == n]
        frac = d["novel"].mean()
        boots = []
        for b in range(args.n_boot):
            rng = np.random.default_rng(derive_seed(cfg["seed"], "novelty", n, b))
            idx = rng.choice(len(d), size=int(min_n), replace=False)
            boots.append(d["novel"].to_numpy()[idx].mean())
        boots = np.array(boots)
        rows.append({"niche": n, "n_species": len(d),
                     "n_novel_species": int(d["novel"].sum()),
                     "frac_novel": round(float(frac), 4),
                     "frac_novel_rarefied_mean": round(float(boots.mean()), 4),
                     "lo": round(float(np.percentile(boots, 2.5)), 4),
                     "hi": round(float(np.percentile(boots, 97.5)), 4),
                     "rarefied_to": int(min_n)})
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}_by_niche.tsv", sep="\t", index=False)

    # --- per rank, placeholder taxa per niche ---
    rrows = []
    for rank in ranks:
        # taxa present in each niche (distinct labels among species of that niche)
        for n in niches:
            d = sp[sp["niche_primary"] == n]
            taxa = d[rank].dropna().unique()
            n_taxa = len(taxa)
            if rank == "gtdb_species":
                n_ph = int(pd.Series(taxa).map(is_placeholder_species).sum())
            else:
                n_ph = int(pd.Series(taxa).map(is_placeholder_label).sum())
            rrows.append({"rank": rank.replace("gtdb_", ""), "niche": n,
                          "n_taxa": n_taxa, "n_placeholder": n_ph,
                          "frac_placeholder": round(n_ph / n_taxa, 4) if n_taxa else 0.0})
    pd.DataFrame(rrows).to_csv(f"{args.out_prefix}_by_rank.tsv", sep="\t", index=False)
    log.info("Novelty: rarefied novel fraction per niche = %s",
             {r["niche"]: r["frac_novel_rarefied_mean"] for r in rows})


if __name__ == "__main__":
    main()
