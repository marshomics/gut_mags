#!/usr/bin/env python3
"""
parse_annotations.py
--------------------
Convert every annotation format into ONE uniform long table so the rest of the
pipeline never has to know which tool produced a feature:

    columns:  genome <TAB> layer <TAB> feature <TAB> count

`count` is the per-genome copy number (number of genes carrying that feature).
Presence/absence is derived later by binarising; keeping counts here means the
abundance-based analyses (CLR ordination, CAZyme/BGC/AMR load) and the
presence-based analyses both read from the same source.

A chunk of genome ids is processed per invocation (array-job friendly). File
paths come from annotation_paths.tsv, which resolve_annotation_paths.py builds
from either a manifest of explicit paths or the {genome} templates in the config,
so nothing about your directory layout is assumed in code.

Supported --type values and the layers they emit:
    kofam      -> ko
    eggnog     -> ko, pfam, cog, cazyme, ec      (cross-check layers)
    dbcan      -> cazyme                          (primary CAZyme source)
    antismash  -> bgc
    amrfinder  -> amr
"""
import argparse
import json
import os
import re

import pandas as pd

from hgn_utils import load_config, get_logger, load_annotation_paths

log = get_logger("parse")


# ----------------------------------------------------------------------------- parsers
def parse_kofam(path, fmt):
    """Return list of KO ids (significant assignments only)."""
    kos = []
    if not os.path.exists(path):
        return kos
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            # detail-tsv: first field is '*' for significant hits, KO in col 3
            if parts[0].strip() == "*" and len(parts) >= 3:
                ko = parts[2].strip()
                if ko.startswith("K"):
                    kos.append(ko)
            elif fmt == "mapper" or (len(parts) == 2):
                # mapper: gene<TAB>KO (only assigned genes listed)
                ko = parts[-1].strip()
                if ko.startswith("K"):
                    kos.append(ko)
            else:
                # space-delimited detail variant
                toks = line.split()
                if toks and toks[0] == "*" and len(toks) >= 3 and toks[2].startswith("K"):
                    kos.append(toks[2])
    return kos


def _split_field(val):
    if val in ("", "-", "NA"):
        return []
    return [x.strip() for x in val.split(",") if x.strip() not in ("", "-")]


def parse_eggnog(path):
    """Yield (layer, feature) tuples from an eggNOG-mapper .annotations file."""
    out = []
    if not os.path.exists(path):
        return out
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            for ko in _split_field(row.get("KEGG_ko", "")):
                out.append(("ko", ko.replace("ko:", "")))
            for pf in _split_field(row.get("PFAMs", "")):
                out.append(("pfam", pf))
            cog = row.get("COG_category", "")
            if cog not in ("", "-"):
                for ch in cog.strip():            # multi-letter categories -> each letter
                    if ch.isalpha():
                        out.append(("cog", ch))
            for cz in _split_field(row.get("CAZy", "")):
                out.append(("cazyme", cz))
            for ec in _split_field(row.get("EC", "")):
                out.append(("ec", ec))
    return out


def parse_dbcan(path):
    """Parse run_dbcan overview.txt; keep families supported by >=2 tools
    (HMMER, dbCAN_sub, DIAMOND). Returns list of CAZy family strings (one per
    gene call, so counts reflect copy number)."""
    fams = []
    if not os.path.exists(path):
        return fams
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    cols = {c.lower(): c for c in df.columns}
    tool_cols = [cols[c] for c in ("hmmer", "dbcan_sub", "diamond") if c in cols]
    fam_re = re.compile(r"([A-Za-z]+\d+)")     # GH13, GT2, PL9, CE1, CBM50, AA10 ...
    for _, r in df.iterrows():
        calls = []
        for tc in tool_cols:
            v = r[tc]
            if v in ("", "-", "N", "NA"):
                continue
            m = fam_re.findall(v)
            if m:
                calls.append(m[0])           # base family of this tool's top call
        # consensus: family found by >=2 tools
        for fam in set(calls):
            if calls.count(fam) >= 2:
                fams.append(fam)
    return fams


def parse_antismash(path):
    """Parse antiSMASH JSON; return one entry per detected BGC region with its
    product class."""
    prods = []
    if not os.path.exists(path):
        return prods
    with open(path) as fh:
        data = json.load(fh)
    for rec in data.get("records", []):
        for feat in rec.get("features", []):
            if feat.get("type") == "region":
                q = feat.get("qualifiers", {})
                for p in q.get("product", []):
                    prods.append(p)
    return prods


def parse_amrfinder(path):
    """Parse AMRFinderPlus TSV; return AMR gene symbols (Element type == AMR)."""
    genes = []
    if not os.path.exists(path):
        return genes
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    et = next((c for c in df.columns if c.lower().replace(" ", "") == "elementtype"), None)
    gs = next((c for c in df.columns if c.lower().replace(" ", "") in
               ("genesymbol", "gene")), None)
    if gs is None:
        return genes
    for _, r in df.iterrows():
        if et is None or r[et].strip().upper() == "AMR":
            genes.append(r[gs].strip())
    return genes


# ----------------------------------------------------------------------------- driver
def emit(genome, items):
    """items: list of (layer, feature). Aggregate to counts."""
    from collections import Counter
    c = Counter(items)
    return [{"genome": genome, "layer": l, "feature": f, "count": n}
            for (l, f), n in c.items()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--type", required=True,
                    choices=["kofam", "eggnog", "dbcan", "antismash", "amrfinder"])
    ap.add_argument("--genome-list", required=True,
                    help="text file, one genome id per line (a chunk)")
    ap.add_argument("--annotation-paths", required=True,
                    help="annotation_paths.tsv written by resolve_annotation_paths.py")
    ap.add_argument("--out", required=True, help="parquet long table for this chunk")
    args = ap.parse_args()

    cfg = load_config(args.config)
    kofam_fmt = cfg["inputs"].get("kofam_format", "detail-tsv")
    # Paths were resolved once, from a manifest or from {genome} templates; this
    # script no longer knows or cares which.
    paths = load_annotation_paths(args.annotation_paths)
    kind = {"kofam": "kofam", "eggnog": "eggnog", "dbcan": "dbcan_overview",
            "antismash": "antismash_json", "amrfinder": "amrfinder_tsv"}[args.type]

    with open(args.genome_list) as fh:
        genomes = [g.strip() for g in fh if g.strip()]

    absent = [g for g in genomes if g not in paths]
    if absent:
        raise SystemExit(f"{len(absent)} genomes in {args.genome_list} are missing from "
                         f"{args.annotation_paths} (e.g. {absent[:3]})")

    rows = []
    for g in genomes:
        p = paths[g].get(kind)
        if not p:
            raise SystemExit(f"no {kind!r} path for genome {g!r} in {args.annotation_paths}")
        if args.type == "kofam":
            rows += emit(g, [("ko", k) for k in parse_kofam(p, kofam_fmt)])
        elif args.type == "eggnog":
            rows += emit(g, parse_eggnog(p))
        elif args.type == "dbcan":
            rows += emit(g, [("cazyme", f) for f in parse_dbcan(p)])
        elif args.type == "antismash":
            rows += emit(g, [("bgc", f) for f in parse_antismash(p)])
        elif args.type == "amrfinder":
            rows += emit(g, [("amr", f) for f in parse_amrfinder(p)])

    out = pd.DataFrame(rows, columns=["genome", "layer", "feature", "count"])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_parquet(args.out, index=False)
    log.info("[%s] %d genomes -> %d feature rows", args.type, len(genomes), len(out))


if __name__ == "__main__":
    main()
