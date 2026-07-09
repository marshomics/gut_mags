#!/usr/bin/env python3
"""
population_turnover.py
----------------------
Western vs non-Western: is the functional repertoire conserved despite species
turnover? Puts taxonomic and functional turnover side by side and measures, for
each shared function, whether it is carried by the SAME or DIFFERENT species in
the two populations.

  taxonomic turnover   Sorensen dissimilarity of species composition, split into
                       turnover (replacement) and nestedness (Baselga); % shared.
  functional turnover  the same on the function sets per layer; % shared functions.
  carrier substitution for each function present in both populations, the Jaccard
                       overlap of its carrier species across populations. Low
                       overlap = same function, different species (functional
                       redundancy expressed as taxonomic substitution).

The headline contrast is taxonomic dissimilarity (expected high) versus
functional dissimilarity (low if functions are conserved). Both are also computed
rarefied to equal species number per population, since Western guts are usually
over-sampled.

Inert (writes a note) unless population.enabled and both populations are present.

Outputs: population_turnover.tsv, carrier_substitution_<layer>.tsv
"""
import argparse
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("pop-turnover")


def baselga_sorensen(A, B):
    a = len(A & B); b = len(A - B); c = len(B - A)
    denom = 2 * a + b + c
    if denom == 0:
        return dict(sorensen=np.nan, turnover=np.nan, nestedness=np.nan, pct_shared=np.nan)
    sor = (b + c) / denom
    sim = (min(b, c) / (a + min(b, c))) if (a + min(b, c)) else 0.0
    return dict(sorensen=round(sor, 4), turnover=round(sim, 4),
                nestedness=round(sor - sim, 4),
                pct_shared=round(100 * a / len(A | B), 2) if (A | B) else np.nan)


def rarefied_sorensen(A, B, k, nboot, seed, tag):
    A, B = list(A), list(B)
    vals = []
    for i in range(nboot):
        rng = np.random.default_rng(derive_seed(seed, "rare_sor", tag, i))
        sa = set(rng.choice(A, k, replace=False)); sb = set(rng.choice(B, k, replace=False))
        vals.append(baselga_sorensen(sa, sb)["sorensen"])
    return float(np.nanmean(vals)) if vals else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    pop = cfg.get("population", {})
    seed = cfg["seed"]
    os.makedirs(args.out_dir, exist_ok=True)
    note = f"{args.out_dir}/population_turnover.tsv"

    sp = pd.read_csv(args.species_table, sep="\t")
    if not pop.get("enabled") or "n_western" not in sp.columns:
        pd.DataFrame([{"status": "population analysis disabled or no western_nonwestern column"}]).to_csv(
            note, sep="\t", index=False)
        log.info("Population turnover inert (not enabled / column absent)."); return
    W = set(sp.loc[sp["n_western"] > 0, "species"])
    NW = set(sp.loc[sp["n_non_western"] > 0, "species"])
    if min(len(W), len(NW)) < pop.get("min_species_per_population", 20):
        pd.DataFrame([{"status": f"too few species per population (W={len(W)}, NW={len(NW)})"}]).to_csv(
            note, sep="\t", index=False)
        log.warning("Too few species per population."); return

    k = min(len(W), len(NW))
    tax = baselga_sorensen(W, NW)
    tax_rare = rarefied_sorensen(W, NW, k, 200, seed, "tax") if pop.get("rarefy") else np.nan

    rows = []
    for layer in cfg["redundancy"]["layers"]:
        path = f"{args.profiles_dir}/prevalence_{layer}.parquet"
        if not os.path.exists(path):
            continue
        prev = pd.read_parquet(path)
        prev = prev[prev["present"] == 1][["species", "feature"]]
        wf = set(prev.loc[prev["species"].isin(W), "feature"])
        nwf = set(prev.loc[prev["species"].isin(NW), "feature"])
        fun = baselga_sorensen(wf, nwf)

        # carrier substitution for shared functions
        wc = prev[prev["species"].isin(W)].groupby("feature")["species"].apply(set)
        nwc = prev[prev["species"].isin(NW)].groupby("feature")["species"].apply(set)
        shared = set(wc.index) & set(nwc.index)
        sub_rows = []
        for f in shared:
            u = wc[f] | nwc[f]
            sub_rows.append({"feature": f, "carrier_jaccard": len(wc[f] & nwc[f]) / len(u) if u else np.nan,
                             "n_carriers_western": len(wc[f]), "n_carriers_nonwestern": len(nwc[f])})
        subdf = pd.DataFrame(sub_rows)
        if len(subdf):
            subdf.sort_values("carrier_jaccard").to_csv(
                f"{args.out_dir}/carrier_substitution_{layer}.tsv", sep="\t", index=False)
        med_carrier_j = float(subdf["carrier_jaccard"].median()) if len(subdf) else np.nan

        rows.append({"layer": layer,
                     "taxonomic_sorensen": tax["sorensen"], "taxonomic_turnover": tax["turnover"],
                     "taxonomic_sorensen_rarefied": round(tax_rare, 4) if tax_rare == tax_rare else np.nan,
                     "pct_shared_species": tax["pct_shared"],
                     "functional_sorensen": fun["sorensen"], "functional_turnover": fun["turnover"],
                     "pct_shared_functions": fun["pct_shared"],
                     "median_carrier_jaccard_shared_fn": round(med_carrier_j, 4) if med_carrier_j == med_carrier_j else np.nan,
                     "interpretation": ("functions conserved despite species turnover"
                                        if (tax["sorensen"] or 0) - (fun["sorensen"] or 0) > 0.2
                                        else "functional and taxonomic turnover comparable")})
    pd.DataFrame(rows).to_csv(note, sep="\t", index=False)
    log.info("Population turnover: taxonomic Sorensen=%.3f; per-layer functional vs taxonomic written",
             tax["sorensen"] if tax["sorensen"] == tax["sorensen"] else float("nan"))


if __name__ == "__main__":
    main()
