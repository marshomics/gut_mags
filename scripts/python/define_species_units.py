#!/usr/bin/env python3
"""
define_species_units.py
-----------------------
Turn the clean per-genome sample sheet into the SPECIES-level analysis units
that every comparative test uses. This is where the strains-per-species and
species-per-niche confounders are first neutralised: each species becomes one
row, carrying its niche occupancy and the metadata needed to weight it fairly.

Outputs:
  species_table.tsv        one row per species: occupancy, niche assignment,
                           n_genomes per niche, host breadth, domain, taxonomy
  occupancy_matrix.tsv     species x niche binary occupancy
  strains_per_species.tsv  full distribution (for the confounder figure)
  species_units_report.json

Definitions (all from config):
  occupies(niche)   : species has >= occupancy.min_genomes_per_niche genomes there
  specialist        : occupies exactly one niche
  generalist        : occupies >= 2 niches
  niche_primary     : niche holding the majority of the species' genomes
  contrast_niche    : the niche label used in the primary differential test,
                      governed by occupancy.primary_contrast
"""
import argparse
import json

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, set_global_seed, provenance_stamp

log = get_logger("species")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True, help="samples.parquet from ingest")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    niches = cfg["inputs"]["niche_levels"]
    occ_cfg = cfg["occupancy"]
    min_g = occ_cfg["min_genomes_per_niche"]

    df = pd.read_parquet(args.samples)
    log.info("%d genomes across %d species", len(df), df["species"].nunique())

    # genomes per species per niche
    counts = (df.groupby(["species", "niche"]).size()
                .unstack(fill_value=0)
                .reindex(columns=niches, fill_value=0))
    counts.columns = [f"n_{c}" for c in niches]

    occ = pd.DataFrame(index=counts.index)
    for c in niches:
        occ[c] = (counts[f"n_{c}"] >= min_g).astype(int)
    occ_sum = occ[niches].sum(axis=1)

    sp = counts.copy()
    sp["n_total"] = sp.sum(axis=1)
    sp["n_niches_occupied"] = occ_sum.values
    sp["specialist"] = (occ_sum == 1).values
    sp["generalist"] = (occ_sum >= 2).values
    # the single niche a specialist belongs to (NaN for generalists)
    sp["specialist_niche"] = np.where(
        sp["specialist"],
        occ[niches].idxmax(axis=1).where(occ_sum == 1),
        np.nan)
    # majority niche (works for everyone)
    sp["niche_primary"] = counts[[f"n_{c}" for c in niches]].idxmax(axis=1).str.replace("n_", "", regex=False)

    # contrast niche per config
    if occ_cfg["primary_contrast"] == "specialists_only":
        sp["contrast_niche"] = sp["specialist_niche"]
    else:
        sp["contrast_niche"] = sp["niche_primary"]

    # taxonomy + placeholder flag (take first genome's taxonomy; species-stable)
    tax_cols = cfg["inputs"]["columns"]["gtdb_ranks"]
    tax = df.groupby("species")[tax_cols + ["domain", "species_is_placeholder"]].first()
    sp = sp.join(tax)

    # host breadth (animal): number of distinct hosts contributing genomes
    host_col = cfg["inputs"]["columns"]["host_common_name"]
    if host_col in df.columns:
        host_breadth = (df[df["niche"] == "animal"]
                        .groupby("species")[host_col]
                        .apply(lambda s: s[s != ""].nunique()))
        sp["animal_host_breadth"] = host_breadth.reindex(sp.index).fillna(0).astype(int)
        # dominant host per species (for host-balanced weighting later)
        dom_host = (df[df["niche"] == "animal"]
                    .groupby("species")[host_col]
                    .agg(lambda s: s[s != ""].mode().iat[0] if (s != "").any() else ""))
        sp["animal_dominant_host"] = dom_host.reindex(sp.index).fillna("")
    else:
        sp["animal_host_breadth"] = 0
        sp["animal_dominant_host"] = ""

    # population (Western vs non-Western) occupancy among human genomes
    if "population" in df.columns:
        hp = df[(df["niche"] == "human") & (df["population"].isin(["western", "non_western"]))]
        pc = (hp.groupby(["species", "population"]).size().unstack(fill_value=0)
              if len(hp) else pd.DataFrame())
        for p in ["western", "non_western"]:
            sp[f"n_{p}"] = (pc[p].reindex(sp.index).fillna(0).astype(int)
                           if (len(pc) and p in pc.columns) else 0)
        def _popassign(w, nw):
            if w > 0 and nw > 0:
                return "shared"
            if w > 0:
                return "western"
            if nw > 0:
                return "non_western"
            return ""
        sp["population"] = [_popassign(w, nw) for w, nw in zip(sp["n_western"], sp["n_non_western"])]
    else:
        sp["n_western"] = 0; sp["n_non_western"] = 0; sp["population"] = ""

    sp = sp.reset_index().rename(columns={"index": "species"})

    # --- write -----------------------------------------------------------------
    sp.to_csv(f"{args.out_prefix}_species_table.tsv", sep="\t", index=False)
    occ.reset_index().to_csv(f"{args.out_prefix}_occupancy_matrix.tsv",
                             sep="\t", index=False)

    spc = df.groupby("species").size().rename("n_genomes").reset_index()
    spc = spc.merge(sp[["species", "contrast_niche", "niche_primary",
                        "specialist", "generalist"]], on="species", how="left")
    spc.sort_values("n_genomes", ascending=False).to_csv(
        f"{args.out_prefix}_strains_per_species.tsv", sep="\t", index=False)

    report = {
        "provenance": provenance_stamp(cfg, {"script": "define_species_units"}),
        "n_species": int(len(sp)),
        "n_specialists": int(sp["specialist"].sum()),
        "n_generalists": int(sp["generalist"].sum()),
        "specialists_per_niche": sp.loc[sp["specialist"], "specialist_niche"]
                                   .value_counts().to_dict(),
        "primary_niche_distribution": sp["niche_primary"].value_counts().to_dict(),
        "strains_per_species": {
            "min": int(spc["n_genomes"].min()),
            "median": float(spc["n_genomes"].median()),
            "mean": float(spc["n_genomes"].mean()),
            "max": int(spc["n_genomes"].max()),
            "singletons": int((spc["n_genomes"] == 1).sum()),
        },
        "placeholder_species": int(sp["species_is_placeholder"].sum()),
        "contrast_mode": occ_cfg["primary_contrast"],
    }
    with open(f"{args.out_prefix}_species_units_report.json", "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    log.info("Species: %d (%d specialists, %d generalists)",
             report["n_species"], report["n_specialists"], report["n_generalists"])
    log.info("Specialists per niche: %s", report["specialists_per_niche"])


if __name__ == "__main__":
    main()
