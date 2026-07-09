#!/usr/bin/env python3
"""
cmh_stratified.py  --  METHOD B: clade-stratified enrichment
------------------------------------------------------------
Cochran-Mantel-Haenszel test of feature presence vs group (focal niche vs
comparator), stratified by clade (genus by default). Stratifying means human
species are only ever compared with non-human species *in the same genus*, so a
difference cannot be an artefact of comparing distantly related lineages. This
is an orthogonal phylogenetic control to phyloglm: phyloglm models the
continuous covariance; CMH matches like-with-like discretely. Agreement between
the two is strong evidence.

Per feature it returns the Mantel-Haenszel common odds ratio, its CI, and the
CMH chi-square p-value, using only strata that contain both groups and reach
the minimum size. Strata where the feature is invariant contribute no
information and are dropped.

A feature with no informative stratum is UNTESTABLE by CMH, not non-significant:
the association is completely confounded with clade at this rank, so no
within-clade statement can be made either way. That case is reported as
status="no_informative_strata" (p = NA) so downstream consensus can tell the two
apart instead of treating an undefined test as a negative one.

Output TSV: feature, mh_log_or, ci_lo, ci_hi, p, n_strata, n_species_used, status
"""
import argparse

import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable

from hgn_utils import load_config, get_logger

log = get_logger("cmh")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--presence", required=True)
    ap.add_argument("--features", required=False)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rank = cfg["stats"]["cmh"]["stratify_rank"].replace("gtdb_", "")
    min_str = cfg["stats"]["cmh"]["min_stratum_size"]

    meta = pd.read_csv(args.analysis, sep="\t")
    strat_col = "genus" if "genus" in meta.columns else rank
    meta[strat_col] = meta[strat_col].fillna("unknown").replace("", "unknown")

    pres = pd.read_parquet(args.presence)
    if args.features:
        feats = [f.strip() for f in open(args.features) if f.strip()]
        pres = pres[pres["feature"].isin(feats)]

    # species x feature presence (dense, tested features only)
    wide = (pres.pivot_table(index="species", columns="feature",
                             values="present", fill_value=0)
                .astype(np.int8))
    meta = meta[meta["species"].isin(wide.index)].copy()
    wide = wide.loc[meta["species"]]
    group = meta["group"].to_numpy()
    strata = meta[strat_col].to_numpy()

    rows = []
    for feat in wide.columns:
        y = wide[feat].to_numpy()
        tables, n_used = [], 0
        for s in np.unique(strata):
            m = strata == s
            if m.sum() < min_str:
                continue
            g, yy = group[m], y[m]
            if (g == 1).sum() == 0 or (g == 0).sum() == 0:
                continue
            a = int(((g == 1) & (yy == 1)).sum())   # focal present
            b = int(((g == 1) & (yy == 0)).sum())   # focal absent
            c = int(((g == 0) & (yy == 1)).sum())   # comp present
            d = int(((g == 0) & (yy == 0)).sum())   # comp absent
            if (a + c) == 0 or (b + d) == 0:        # feature invariant in stratum
                continue
            tables.append([[a, b], [c, d]])
            n_used += (a + b + c + d)
        rec = {"feature": feat, "mh_log_or": np.nan, "ci_lo": np.nan,
               "ci_hi": np.nan, "p": np.nan, "n_strata": len(tables),
               "n_species_used": n_used, "status": "no_informative_strata"}
        if len(tables) >= 1:
            try:
                st = StratifiedTable([np.array(t) for t in tables])
                orr = st.oddsratio_pooled
                lo, hi = st.oddsratio_pooled_confint()
                # Complete separation within every stratum gives OR 0 or inf: the
                # test is still valid (p is finite), the point estimate is not.
                # Keep the p-value and record the direction via a Haldane-corrected
                # OR so the effect is bounded but signed, never dropped.
                p = float(st.test_null_odds().pvalue)
                if np.isfinite(orr) and orr > 0:
                    rec["mh_log_or"] = float(np.log(orr))
                else:
                    num = sum(t[0][0] + 0.5 for t in tables) * sum(t[1][1] + 0.5 for t in tables)
                    den = sum(t[0][1] + 0.5 for t in tables) * sum(t[1][0] + 0.5 for t in tables)
                    rec["mh_log_or"] = float(np.log(num / den))
                rec["ci_lo"] = float(np.log(lo)) if np.isfinite(lo) and lo > 0 else np.nan
                rec["ci_hi"] = float(np.log(hi)) if np.isfinite(hi) and hi > 0 else np.nan
                rec["p"] = p
                rec["status"] = "ok" if np.isfinite(p) else "test_failed"
            except Exception:
                rec["status"] = "test_failed"
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    n_ok = int((out["status"] == "ok").sum())
    log.info("CMH: %d features, stratified by %s; %d testable, %d with no informative "
             "stratum (completely confounded with clade), %d test failures",
             len(out), strat_col, n_ok,
             int((out["status"] == "no_informative_strata").sum()),
             int((out["status"] == "test_failed").sum()))
    if n_ok == 0 and len(out):
        log.warning("No feature had an informative stratum at rank '%s'. Niche is "
                    "fully nested within clade at this rank; consider a coarser "
                    "stats.cmh.stratify_rank. Consensus will fall back to the "
                    "other methods and flag these features as CMH-untestable.",
                    strat_col)


if __name__ == "__main__":
    main()
