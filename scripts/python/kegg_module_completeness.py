#!/usr/bin/env python3
"""
kegg_module_completeness.py
---------------------------
Compute KEGG module completeness per species from species-level KO presence.

KEGG module DEFINITION strings encode pathway logic:
    space  = AND between steps (each space-separated block is one reaction step)
    ' , '  = OR between alternatives
    '+'    = AND (essential complex subunit / required combination)
    '-'    = the following KO/component is OPTIONAL (non-essential)
    (...)  = grouping
Completeness is the fraction of top-level steps that are satisfiable from the
species' KO set. This is the same "stepwise block completeness" used by anvi'o
and MicrobeAnnotator, so results are comparable to published practice.

A module is called "present" when completeness >= module_completeness_threshold
(config). Continuous completeness is also kept for the comparative tests, which
are more powerful than a binary call.

Inputs:
  --ko-prevalence  prevalence_ko.parquet (species, feature=KO, present)
  --module-def     TSV: module_id <TAB> name <TAB> definition   (see resources/)
Outputs:
  module_completeness.parquet  species, module, name, completeness, present
"""
import argparse
import re

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("modules")
_TOKEN = re.compile(r"K\d{5}|M\d{5}|[()+,\- ]")


def tokenize(defn):
    # keep only KO/module ids and structural operators
    return [t for t in _TOKEN.findall(defn)]


def split_top(tokens, sep):
    """Split a token list on a separator that is at parenthesis depth 0."""
    depth, cur, parts = 0, [], []
    for t in tokens:
        if t == "(":
            depth += 1; cur.append(t)
        elif t == ")":
            depth -= 1; cur.append(t)
        elif t == sep and depth == 0:
            parts.append(cur); cur = []
        else:
            cur.append(t)
    parts.append(cur)
    return parts


def strip_parens(tokens):
    if tokens and tokens[0] == "(" and tokens[-1] == ")":
        depth = 0
        for i, t in enumerate(tokens):
            depth += (t == "(") - (t == ")")
            if depth == 0 and i < len(tokens) - 1:
                return tokens
        return tokens[1:-1]
    return tokens


def evaluate(tokens, kos):
    """Return True if the (sub)expression is satisfied by KO set `kos`."""
    tokens = [t for t in tokens if t != " "]
    tokens = strip_parens(tokens)
    tokens = [t for t in tokens if t != " "]
    if not tokens:
        return True
    # OR (lowest precedence)
    ors = split_top(tokens, ",")
    if len(ors) > 1:
        return any(evaluate(o, kos) for o in ors)
    # AND via '+'
    plus = split_top(tokens, "+")
    if len(plus) > 1:
        return all(evaluate(p, kos) for p in plus)
    # optional component '-': required part is the token(s) before '-'
    minus = split_top(tokens, "-")
    if len(minus) > 1:
        return evaluate(minus[0], kos)         # optional tail ignored
    # single KO / module reference
    t = tokens[0]
    if t.startswith("K") or t.startswith("M"):
        return t in kos
    return evaluate(tokens[1:], kos) if len(tokens) > 1 else False


def module_completeness(defn, kos):
    toks = [t for t in tokenize(defn)]
    toks = strip_parens(toks)
    steps = [s for s in split_top(toks, " ") if any(x.strip() for x in s)]
    steps = [s for s in steps if [x for x in s if x not in ("(", ")")]]
    if not steps:
        return 0.0
    done = sum(1 for s in steps if evaluate(s, kos))
    return done / len(steps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--ko-prevalence", required=True)
    ap.add_argument("--module-def", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    thr = cfg["functional_annotation"]["module_completeness_threshold"]

    ko = pd.read_parquet(args.ko_prevalence)
    ko = ko[ko["present"] == 1]
    sp_kos = ko.groupby("species")["feature"].apply(set)

    mods = pd.read_csv(args.module_def, sep="\t", dtype=str, keep_default_na=False)
    # expected columns: module_id, name, definition
    mcol = mods.columns[0]; ncol = mods.columns[1]; dcol = mods.columns[2]

    rows = []
    for sp, kos in sp_kos.items():
        for _, m in mods.iterrows():
            c = module_completeness(m[dcol], kos)
            if c > 0:
                # uniform prevalence-style schema so this layer slots into the
                # same downstream code as ko/pfam/cazyme/... :
                #   prevalence == module completeness (continuous, [0,1])
                #   present    == completeness >= threshold
                rows.append({"species": sp, "feature": m[mcol],
                             "name": m[ncol], "prevalence": round(c, 4),
                             "present": int(c >= thr),
                             "mean_copies": round(c, 4), "n_genomes": -1})
    out = pd.DataFrame(rows, columns=["species", "feature", "name", "prevalence",
                                      "present", "mean_copies", "n_genomes"])
    # The prevalence parquets are concatenated across layers (combine_presence),
    # so every layer must share one schema. Module names go to a sidecar TSV.
    names = out[["feature", "name"]].drop_duplicates()
    names.to_csv(args.out.replace(".parquet", "_names.tsv"), sep="\t", index=False)
    out = out[["species", "feature", "prevalence", "present", "mean_copies", "n_genomes"]]
    out.to_parquet(args.out, index=False)
    log.info("Module completeness: %d species, %d modules evaluated",
             sp_kos.shape[0], mods.shape[0])


if __name__ == "__main__":
    main()
