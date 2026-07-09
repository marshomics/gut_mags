#!/usr/bin/env python3
"""
synthesis_catalogues.py
-----------------------
Integrate the stages into the three "human-specific" catalogues the manuscript
needs, each with multi-line evidence.

  human_specific_species.tsv   niche specialists for the focal niche, annotated
      with strain count, novelty and whether their genus/family is a significant
      human indicator taxon (IndVal). Tiered by how many lines agree.
  human_specific_genes.tsv     consensus signature features enriched in the focal
      niche, across layers, with the number of supporting methods and the
      adaptation mode (acquisition / loss / ...).
  human_specific_functions.tsv categories enriched by both ORA and GSEA, plus the
      curated ecological pressures that are human-enriched.

Output also: catalogue_counts.json (headline numbers for the master figure).
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("synth-cat")


def load(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--indicator-dir", required=True)
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--adaptation", required=True)
    ap.add_argument("--enrichment-top", required=True)
    ap.add_argument("--ecological", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    focal = cfg["inputs"]["focal_niche"]
    contrast = cfg["synthesis"]["focal_contrast"]
    os.makedirs(args.out_dir, exist_ok=True)

    # --- species ---------------------------------------------------------------
    sp = pd.read_csv(args.species_table, sep="\t")
    hs = sp[sp["specialist_niche"] == focal].copy()
    # human-indicator genera/families (IndVal) for tiering
    ind_taxa = set()
    for f in glob.glob(f"{args.indicator_dir}/indicator_*.tsv"):
        d = load(f)
        if len(d) and "niche_combination" in d.columns:
            ind_taxa |= set(d.loc[(d["niche_combination"] == focal) &
                                  (d.get("q.value", 1) < 0.05), "taxon"])
    hs["genus_is_human_indicator"] = hs["gtdb_genus"].isin(ind_taxa)
    hs["family_is_human_indicator"] = hs["gtdb_family"].isin(ind_taxa)
    hs["tier"] = ["high" if (g or fam) else "moderate"
                  for g, fam in zip(hs["genus_is_human_indicator"], hs["family_is_human_indicator"])]
    cols = ["species", "gtdb_genus", "gtdb_family", "n_" + focal,
            "species_is_placeholder", "genus_is_human_indicator",
            "family_is_human_indicator", "tier"]
    cols = [c for c in cols if c in hs.columns]
    hs[cols].sort_values("tier").to_csv(f"{args.out_dir}/human_specific_species.tsv",
                                        sep="\t", index=False)

    # --- genes / features ------------------------------------------------------
    sig = pd.read_csv(args.signatures, sep="\t")
    sig["consensus_signature"] = sig["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    genes = sig[sig["consensus_signature"] & (sig["direction"] == f"{focal}_enriched")].copy()
    adapt = load(args.adaptation)
    if len(adapt):
        genes = genes.merge(adapt[["feature", "layer", "mode", "human_prev", "comparator_prev"]],
                            on=["feature", "layer"], how="left")
    keep = [c for c in ["feature", "layer", "consensus_log2or", "n_methods_support",
                        "mode", "human_prev", "comparator_prev"] if c in genes.columns]
    genes[keep].sort_values(["layer", "consensus_log2or"], ascending=[True, False]).to_csv(
        f"{args.out_dir}/human_specific_genes.tsv", sep="\t", index=False)

    # --- functions -------------------------------------------------------------
    et = load(args.enrichment_top)
    et = et[et["contrast"] == contrast] if "contrast" in et.columns else et
    eco = load(args.ecological)
    func_rows = []
    if len(et):
        for _, r in et.iterrows():
            func_rows.append({"source": "enrichment", "system": r.get("system"),
                              "category": r.get("category_name"), "direction": r.get("direction"),
                              "score": r.get("NES"), "q": r.get("ora_q")})
    if len(eco) and "enriched" in eco.columns:
        for _, r in eco[eco["enriched"].astype(bool)].iterrows():
            func_rows.append({"source": "ecological_pressure", "system": "pressure",
                              "category": r.get("name"), "direction": "up",
                              "score": r.get("fold_enrichment"), "q": r.get("q")})
    pd.DataFrame(func_rows).to_csv(f"{args.out_dir}/human_specific_functions.tsv",
                                   sep="\t", index=False)

    counts = {"focal_niche": focal,
              "n_human_specific_species": int(len(hs)),
              "n_high_tier_species": int((hs["tier"] == "high").sum()),
              "n_human_specific_genes": int(len(genes)),
              "n_human_specific_functions": int(len(func_rows)),
              "genes_per_layer": genes["layer"].value_counts().to_dict() if len(genes) else {}}
    json.dump(counts, open(f"{args.out_dir}/catalogue_counts.json", "w"), indent=2, default=str)
    log.info("Catalogues: %d species, %d genes, %d functions",
             counts["n_human_specific_species"], counts["n_human_specific_genes"],
             counts["n_human_specific_functions"])


if __name__ == "__main__":
    main()
