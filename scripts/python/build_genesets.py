#!/usr/bin/env python3
"""
build_genesets.py
-----------------
Build category memberships (gene sets) for one feature layer: which higher-level
functional category each feature belongs to, for each configured system. These
sets are the same regardless of contrast, so they are built once per layer and
intersected with the tested universe inside the enrichment step.

Sources, in order of defensibility:
  * deterministic (no external data): cazy_class (family prefix), cog_group (COG
    letter -> super-group), bgc_group (antiSMASH product -> broad class);
  * kegg_module: parsed from the KEGG module definitions already fetched;
  * file-based: kegg_pathway and any other system with a map file in canonical_dir
    (feature <tab> category_id [<tab> category_name]).
A system with no available source is skipped with a warning rather than guessed.

Output: genesets_<layer>.tsv  (system, category_id, category_name, feature)
"""
import argparse
import os
import re

import pandas as pd

from hgn_utils import load_config, get_logger
from kegg_module_completeness import tokenize  # reuse KO token extraction

log = get_logger("genesets")

COG_GROUP = {
    **{l: ("Information storage and processing", l) for l in "JAKLB"},
    **{l: ("Cellular processes and signaling", l) for l in "DYVTMNZWUO"},
    **{l: ("Metabolism", l) for l in "CGEFHIPQ"},
    **{l: ("Poorly characterized", l) for l in "RS"},
}
BGC_GROUP = {
    "NRPS": "NRPS", "NRPS-like": "NRPS",
    "T1PKS": "PKS", "T2PKS": "PKS", "T3PKS": "PKS", "transAT-PKS": "PKS", "hglE-KS": "PKS",
    "terpene": "Terpene",
    "RiPP-like": "RiPP", "lanthipeptide": "RiPP", "lassopeptide": "RiPP",
    "thiopeptide": "RiPP", "sactipeptide": "RiPP", "bacteriocin": "RiPP",
    "lanthipeptide-class-i": "RiPP", "lanthipeptide-class-ii": "RiPP",
    "siderophore": "Siderophore", "betalactone": "Betalactone",
    "arylpolyene": "Arylpolyene", "ectoine": "Ectoine", "phosphonate": "Phosphonate",
}


def deterministic(layer, system, features):
    rows = []
    if system == "cazy_class":
        for f in features:
            m = re.match(r"^([A-Za-z]+)", str(f))
            if m:
                rows.append((m.group(1), m.group(1), f))
    elif system == "cog_group":
        for f in features:
            g = COG_GROUP.get(str(f))
            if g:
                rows.append((g[0], g[0], f))
    elif system == "bgc_group":
        for f in features:
            rows.append((BGC_GROUP.get(str(f), "other"), BGC_GROUP.get(str(f), "other"), f))
    return rows


def kegg_module_sets(module_def_path, features):
    """KO -> module membership from module definitions (KOs cited in each def)."""
    rows = []
    if not os.path.exists(module_def_path):
        return rows
    mods = pd.read_csv(module_def_path, sep="\t", dtype=str, keep_default_na=False)
    mcol, ncol, dcol = mods.columns[0], mods.columns[1], mods.columns[2]
    feat = set(map(str, features))
    for _, m in mods.iterrows():
        kos = {t for t in tokenize(m[dcol]) if t.startswith("K")}
        for ko in kos & feat:
            rows.append((m[mcol], m[ncol], ko))
    return rows


def file_sets(path, features):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    feat = set(map(str, features))
    cid = d.columns[1]
    cname = d.columns[2] if d.shape[1] > 2 else d.columns[1]
    for _, r in d.iterrows():
        if r[d.columns[0]] in feat:
            rows.append((r[cid], r[cname], r[d.columns[0]]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--layer", required=True)
    ap.add_argument("--prevalence", required=True, help="prevalence_<layer>.parquet (feature universe)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    en = cfg["enrichment"]
    systems = en["systems"].get(args.layer, [])
    cdir = en["canonical_dir"]
    features = pd.read_parquet(args.prevalence)["feature"].astype(str).unique().tolist()

    out_rows = []
    for system in systems:
        if system in ("cazy_class", "cog_group", "bgc_group"):
            rows = deterministic(args.layer, system, features)
        elif system == "kegg_module":
            rows = kegg_module_sets(cfg["references"]["kegg_module_def"], features)
        elif system == "kegg_pathway":
            rows = file_sets(os.path.join(cdir, "kegg_pathway.tsv"), features)
        elif system == "cazy_substrate":
            rows = file_sets(en.get("cazy_substrate_map"), features)
        elif system == "amr_class":
            rows = file_sets(en.get("amr_class_map"), features)
        elif system == "pfam_clan":
            rows = file_sets(en.get("pfam_clan_map"), features)
        else:
            rows = file_sets(os.path.join(cdir, f"{system}.tsv"), features)
        if not rows:
            log.warning("Layer %s system %s: no membership source found, skipped.",
                        args.layer, system)
            continue
        for cid, cname, feat in rows:
            out_rows.append({"system": system, "category_id": cid,
                             "category_name": cname, "feature": feat})
        log.info("Layer %s system %s: %d memberships, %d categories",
                 args.layer, system, len(rows), len({r[0] for r in rows}))

    pd.DataFrame(out_rows, columns=["system", "category_id", "category_name", "feature"]).to_csv(
        args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
