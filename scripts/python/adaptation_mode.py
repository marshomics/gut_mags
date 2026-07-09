#!/usr/bin/env python3
"""
adaptation_mode.py
------------------
Decompose how each human signature is achieved, the mechanistic answer to "what
drives niche adaptation". For every consensus human-signature feature it compares
the species-level prevalence in human vs comparator niches and assigns a mode:

  acquisition          enriched in human, rare in comparator (gene gain / HGT)
  differential_retention  enriched in human, also moderately present elsewhere
  loss_in_human        depleted in human, common in comparator (gene loss)
  depletion            depleted in human, also uncommon elsewhere

Sequence-level change is added from the HyPhy selection scan (a separate axis:
the gene is shared but evolves differently), and genome architecture (size, GC,
coding density) from PGLS. The output is the relative contribution of
acquisition, loss, selection and architecture, which together are the answer.

Outputs:
  adaptation_mode_features.tsv   per feature: prevalences, mode
  adaptation_mode_summary.tsv    counts per layer x mode
  adaptation_drivers.json        headline breakdown + architecture + selection
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("adapt-mode")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--signatures", required=True, help="combined signatures_all for the focal contrast")
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--profiles-dir", required=True)
    ap.add_argument("--pgls", default=None)
    ap.add_argument("--selection", default=None, help="selection aggregate tsv (optional)")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    am = cfg["synthesis"]["adaptation_mode"]
    focal = cfg["inputs"]["focal_niche"]
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)

    sp = pd.read_csv(args.species_table, sep="\t")
    sp["niche"] = sp["specialist_niche"].fillna(sp["niche_primary"])
    human = set(sp.loc[sp["niche"] == focal, "species"])
    comparator = set(sp.loc[sp["niche"].isin([n for n in cfg["inputs"]["niche_levels"]
                                              if n != focal]), "species"])

    sig = pd.read_csv(args.signatures, sep="\t")
    sig["consensus_signature"] = sig["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])
    sig = sig[sig["consensus_signature"]].copy()

    rows = []
    for layer, sub in sig.groupby("layer"):
        path = f"{args.profiles_dir}/prevalence_{layer}.parquet"
        if not os.path.exists(path):
            continue
        prev = pd.read_parquet(path)
        prev = prev[prev["feature"].isin(set(sub["feature"])) & (prev["present"] == 1)]
        hp = prev[prev["species"].isin(human)].groupby("feature")["species"].nunique()
        cp = prev[prev["species"].isin(comparator)].groupby("feature")["species"].nunique()
        nH, nC = max(len(human), 1), max(len(comparator), 1)
        for _, r in sub.iterrows():
            f = r["feature"]
            human_prev = hp.get(f, 0) / nH
            comp_prev = cp.get(f, 0) / nC
            direction = str(r.get("direction", ""))
            if direction.endswith("enriched"):
                if comp_prev <= am["comparator_prev_gain_max"] and human_prev >= am["human_prev_min"]:
                    mode = "acquisition"
                else:
                    mode = "differential_retention"
            else:
                if comp_prev >= am["comparator_prev_loss_min"]:
                    mode = "loss_in_human"
                else:
                    mode = "depletion"
            rows.append({"feature": f, "layer": layer, "direction": direction,
                         "human_prev": round(human_prev, 3),
                         "comparator_prev": round(comp_prev, 3),
                         "consensus_log2or": r.get("consensus_log2or"), "mode": mode})
    feats = pd.DataFrame(rows)
    feats.to_csv(f"{args.out_prefix}_features.tsv", sep="\t", index=False)
    summ = (feats.groupby(["layer", "mode"]).size().rename("n").reset_index()
            if len(feats) else pd.DataFrame(columns=["layer", "mode", "n"]))
    summ.to_csv(f"{args.out_prefix}_summary.tsv", sep="\t", index=False)

    drivers = {"n_signatures": int(len(feats)),
               "mode_counts": feats["mode"].value_counts().to_dict() if len(feats) else {}}
    # genome architecture from PGLS (human term)
    if args.pgls and os.path.exists(args.pgls):
        pg = pd.read_csv(args.pgls, sep="\t")
        arch = pg[pg["term"].astype(str).str.contains("niche", case=False)]
        drivers["genome_architecture_pgls"] = arch[["trait", "term", "estimate", "p"]].to_dict("records")
    # sequence selection count
    if args.selection and os.path.exists(args.selection):
        try:
            sel = pd.read_csv(args.selection, sep="\t")
            drivers["n_families_positive_selection"] = int(
                sel.get("busted_sig", pd.Series(dtype=bool)).astype(str).isin(["True", "TRUE", "1"]).sum())
            drivers["n_families_intensified_selection"] = int(
                (pd.to_numeric(sel.get("relax_K"), errors="coerce") > 1).sum())
        except Exception:
            pass
    json.dump(drivers, open(f"{args.out_prefix}_drivers.json", "w"), indent=2, default=str)
    log.info("Adaptation mode: %d signatures; modes=%s", len(feats), drivers["mode_counts"])


if __name__ == "__main__":
    main()
