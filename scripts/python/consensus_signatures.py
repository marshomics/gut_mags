#!/usr/bin/env python3
"""
consensus_signatures.py
-----------------------
Combine the differential methods into a single, conservative call. Effect sizes
from every method are expressed as log2 odds ratios so they are directly
comparable; a positive value means enriched in the focal (human) niche.

Method support is defined per feature as:
  * phyloglm   : BH q < alpha and |log2 OR| >= min_abs_log2_or
  * cmh        : BH q < alpha and |log2 OR| >= min_abs_log2_or
  * resampling : bootstrap CI of the prevalence difference excludes 0
                 and sign_consistency >= 0.95
  * scoary     : Fisher q < alpha, empirical p < alpha, enough supporting pairs

Requiring agreement of a phylogenetic-covariance model (phyloglm), a clade-
matched stratified test (CMH), a balance-corrected bootstrap (resampling) and a
pan-GWAS with the pairwise-comparisons correction (Scoary2) means a call
survives every confounder control simultaneously, not just one.

APPLICABILITY IS NOT SIGNIFICANCE. A method can be undefined for a feature: CMH
has no informative stratum when niche is completely confounded with clade at the
stratifying rank; Scoary drops invariant features; phyloglm can fail to
converge. Treating "undefined" as "disagrees" would silently veto exactly the
features confined to niche-specific clades - the strongest biology - and leave
no trace of having done so. A feature is therefore a CONSENSUS signature when
  (i)   every REQUIRED method that is APPLICABLE to it is significant,
  (ii)  at least stats.consensus.min_applicable_methods are applicable,
  (iii) stats.consensus.anchor_method is applicable and significant,
  (iv)  (if require_same_direction) all supporting methods point the same way.
tier1_consensus = all required methods applicable and significant.
tier2_consensus_partial = passes (i)-(iv) but >=1 required method was untestable;
the methods_applicable / methods_untestable columns say which, per feature.
Set stats.consensus.strict_all_methods=true to restrict calls to tier1.

Outputs:
  signatures_<layer>_<contrast>.tsv   per-feature stats + consensus call + tier
  signatures_<layer>_<contrast>.json  summary incl. positive-control check
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from hgn_utils import load_config, get_logger

log = get_logger("consensus")
LN2 = np.log(2.0)


def bh(p):
    p = pd.to_numeric(p, errors="coerce")
    q = np.full(len(p), np.nan)
    m = p.notna().to_numpy()
    if m.sum():
        q[m] = multipletests(p[m].to_numpy(), method="fdr_bh")[1]
    return q


def positive_label(contrast):
    """The label that is the positive (group-1, 'enriched in') side of a contrast."""
    if contrast == "host_vs_free":
        return "host"
    return contrast.split("_vs_")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--phyloglm", required=True)
    ap.add_argument("--cmh", required=True)
    ap.add_argument("--resampling", required=True)
    ap.add_argument("--scoary", default=None, help="optional Scoary2 parsed table (4th method)")
    ap.add_argument("--layer", required=True)
    ap.add_argument("--contrast", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    alpha = cfg["stats"]["fdr_alpha"]
    eff = cfg["stats"]["min_abs_log2_or"]
    ccfg = cfg["stats"]["consensus"]
    req = list(ccfg["require_methods"])
    same_dir = ccfg["require_same_direction"]
    min_app = int(ccfg.get("min_applicable_methods", 2))
    anchor = ccfg.get("anchor_method", "phyloglm")
    strict = bool(ccfg.get("strict_all_methods", False))
    pos = positive_label(args.contrast)
    focal = cfg["inputs"]["focal_niche"]
    # sensitivity: with phylogenetic control OFF, drop the phylogeny-aware
    # methods so the report shows how many more calls a naive analysis would make.
    if cfg.get("flags", {}).get("phylogenetic_control") == "off":
        req = [m for m in req if m not in ("phyloglm", "cmh")] or ["resampling"]
        log.warning("phylogenetic_control=off: consensus uses %s only", req)

    pg = pd.read_csv(args.phyloglm, sep="\t")
    cm = pd.read_csv(args.cmh, sep="\t")
    rs = pd.read_csv(args.resampling, sep="\t")
    scoary = None
    if args.scoary and os.path.exists(args.scoary) and os.path.getsize(args.scoary):
        sc = pd.read_csv(args.scoary, sep="\t")
        if len(sc):
            sc = sc.copy()
            sc["scoary_tested"] = True   # a feature absent from this table was not tested
            scoary = sc

    # Express every method's effect as a log2 OR and clip to a finite,
    # reportable range (OR in ~[1e-3, 1e3]). Complete separation (feature in all
    # focal, none comparator) yields +/-inf, which clips to +/-CLIP: the signal
    # is kept at the bound and its DIRECTION is preserved, rather than being lost
    # to NaN. Genuine missing values (no test possible) stay NaN.
    CLIP = 10.0

    def clean_log2(x):
        x = pd.to_numeric(x, errors="coerce") / LN2
        return x.clip(-CLIP, CLIP)   # np.clip maps +/-inf to +/-CLIP

    pg["pg_log2or"] = clean_log2(pg["estimate_log_or"])
    pg["pg_q"] = bh(pg["p"])
    cm["cmh_log2or"] = clean_log2(cm["mh_log_or"])
    cm["cmh_q"] = bh(cm["p"])
    rs["rs_log2or"] = clean_log2(rs["median_log_or"])
    if "status" not in cm.columns:      # tolerate older cmh outputs
        cm["status"] = np.where(cm["p"].notna(), "ok", "no_informative_strata")

    df = (pg[["feature", "pg_log2or", "p", "pg_q", "status"]]
          .rename(columns={"p": "pg_p", "status": "pg_status"})
          .merge(cm[["feature", "cmh_log2or", "p", "cmh_q", "n_strata", "status"]]
                 .rename(columns={"p": "cmh_p", "status": "cmh_status",
                                  "n_strata": "cmh_n_strata"}),
                 on="feature", how="outer")
          .merge(rs[["feature", "median_prev_diff", "ci_lo", "ci_hi",
                     "rs_log2or", "sign_consistency"]], on="feature", how="outer"))
    if scoary is not None:
        df = df.merge(scoary[["feature", "scoary_log2or", "scoary_fisher_q",
                              "scoary_empirical_p", "scoary_supporting",
                              "scoary_best_p", "scoary_sig", "scoary_dir",
                              "scoary_tested"]],
                      on="feature", how="outer")

    # ---- applicability: could this method be evaluated for this feature? -----
    df["pg_applicable"] = (df["pg_status"].astype(str) == "ok") & df["pg_p"].notna()
    df["cmh_applicable"] = (df["cmh_status"].astype(str) == "ok") & df["cmh_p"].notna()
    df["rs_applicable"] = df["ci_lo"].notna() & df["ci_hi"].notna() & df["sign_consistency"].notna()
    if scoary is not None:
        df["scoary_applicable"] = df["scoary_tested"].fillna(False).astype(bool)

    # ---- significance (only meaningful where applicable) ---------------------
    df["pg_sig"] = df["pg_applicable"] & (df["pg_q"] < alpha) & (df["pg_log2or"].abs() >= eff)
    df["cmh_sig"] = df["cmh_applicable"] & (df["cmh_q"] < alpha) & (df["cmh_log2or"].abs() >= eff)
    df["rs_sig"] = (df["rs_applicable"] & ((df["ci_lo"] > 0) | (df["ci_hi"] < 0)) &
                    (df["sign_consistency"] >= 0.95))
    if scoary is not None:
        df["scoary_sig"] = (df["scoary_applicable"] &
                            df["scoary_sig"].fillna(False).astype(bool))

    def sgn(x):
        return np.sign(x).replace(0, np.nan)
    df["pg_dir"] = sgn(df["pg_log2or"])
    df["cmh_dir"] = sgn(df["cmh_log2or"])
    df["rs_dir"] = sgn(df["median_prev_diff"])
    if scoary is not None:
        df["scoary_dir"] = sgn(df["scoary_log2or"])

    flag = {"phyloglm": "pg_sig", "cmh": "cmh_sig", "resampling": "rs_sig",
            "scoary": "scoary_sig"}
    appc = {"phyloglm": "pg_applicable", "cmh": "cmh_applicable",
            "resampling": "rs_applicable", "scoary": "scoary_applicable"}
    dcol = {"phyloglm": "pg_dir", "cmh": "cmh_dir", "resampling": "rs_dir",
            "scoary": "scoary_dir"}
    # a method is "required" only if it ran for this contrast (input present)
    available = ["phyloglm", "cmh", "resampling"] + (["scoary"] if scoary is not None else [])
    required = [m for m in req if m in available]
    if anchor not in required:          # anchor must exist; else use the first required
        anchor = required[0]

    A = df[[appc[m] for m in required]].fillna(False).to_numpy(bool)   # applicable
    S = df[[flag[m] for m in required]].fillna(False).to_numpy(bool)   # significant
    n_app = A.sum(axis=1)
    n_sup = S.sum(axis=1)
    ai = required.index(anchor)

    # (i) every applicable required method is significant; (ii) enough applicable;
    # (iii) the anchor is applicable and significant.
    pass_all = (~A | S).all(axis=1) & (n_app >= min_app) & A[:, ai] & S[:, ai]

    if same_dir:
        # directions must agree among the methods that BOTH ran and were significant
        dirs = df[[dcol[m] for m in required]].to_numpy(float)
        d_sig = np.where(S, dirs, np.nan)
        with np.errstate(invalid="ignore"):
            consistent = np.array([
                (np.nansum(np.abs(r)) > 0) and (len(set(r[~np.isnan(r)])) == 1)
                for r in d_sig])
        pass_all &= consistent

    df["n_methods_applicable"] = n_app
    df["n_methods_support"] = n_sup
    df["methods_applicable"] = [",".join(m for m, a in zip(required, row) if a) for row in A]
    df["methods_untestable"] = [",".join(m for m, a in zip(required, row) if not a) for row in A]
    all_app = n_app == len(required)
    if strict:
        pass_all &= all_app
    df["consensus_signature"] = pass_all

    # direction labelled by the contrast's positive side (e.g. human, host, free);
    # taken from the applicable methods only, and "ns" when there is no direction.
    any_dir = pd.Series(np.nan, index=df.index)
    for m in required:
        d = df[dcol[m]].where(df[appc[m]].fillna(False))
        any_dir = any_dir.fillna(d)
    df["direction"] = np.where(any_dir > 0, f"{pos}_enriched",
                        np.where(any_dir < 0, f"{pos}_depleted", "ns"))
    log2_cols = [c for c in ["pg_log2or", "cmh_log2or", "rs_log2or", "scoary_log2or"]
                 if c in df.columns]
    df["consensus_log2or"] = df[log2_cols].mean(axis=1, skipna=True)

    df["tier"] = np.where(pass_all & all_app, "tier1_consensus",
                  np.where(pass_all, "tier2_consensus_partial",
                  np.where(n_sup >= 1, "tier3_suggestive", "ns")))

    df.sort_values(["consensus_signature", "consensus_log2or"],
                   ascending=[False, False]).to_csv(
        f"{args.out_prefix}.tsv", sep="\t", index=False)

    # positive-control check only applies to the human-enriched contrasts
    pc = cfg["stats"]["positive_controls"].get(args.layer, []) if pos == focal else []
    pc_hits = {f: bool(df.loc[df["feature"] == f, "consensus_signature"].any()
                       and (df.loc[df["feature"] == f, "direction"] == f"{pos}_enriched").any())
               for f in pc}

    n_enr = int(((df["consensus_signature"]) & (df["direction"] == f"{pos}_enriched")).sum())
    n_dep = int(((df["consensus_signature"]) & (df["direction"] == f"{pos}_depleted")).sum())
    untestable = {m: int((~df[appc[m]].fillna(False)).sum()) for m in required}
    summary = {
        "layer": args.layer, "contrast": args.contrast, "positive_label": pos,
        "n_tested": int(df["feature"].nunique()),
        "n_consensus": int(df["consensus_signature"].sum()),
        "n_tier1_all_methods": int((df["tier"] == "tier1_consensus").sum()),
        "n_tier2_partial": int((df["tier"] == "tier2_consensus_partial").sum()),
        "n_enriched": n_enr, "n_depleted": n_dep,
        "n_human_enriched": n_enr, "n_human_depleted": n_dep,   # back-compat keys
        "methods_used": required,
        "n_untestable_per_method": untestable,
        "median_methods_applicable": float(df["n_methods_applicable"].median()),
        "scoary_included": scoary is not None,
        "n_scoary_sig": int(df["scoary_sig"].sum()) if "scoary_sig" in df.columns else None,
        "tier_counts": df["tier"].value_counts().to_dict(),
        "positive_controls": pc_hits,
        "positive_controls_missing": [f for f, v in pc_hits.items() if not v],
        "params": {"alpha": alpha, "min_abs_log2_or": eff,
                   "require_methods": req, "required_for_contrast": required,
                   "require_same_direction": same_dir,
                   "min_applicable_methods": min_app, "anchor_method": anchor,
                   "strict_all_methods": strict},
    }
    with open(f"{args.out_prefix}.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("%s/%s: %d consensus signatures (%d enriched, %d depleted); "
             "tier1=%d tier2(partial)=%d; untestable per method: %s",
             args.layer, args.contrast, summary["n_consensus"],
             n_enr, n_dep, summary["n_tier1_all_methods"],
             summary["n_tier2_partial"], untestable)
    if summary["positive_controls_missing"]:
        log.warning("Positive controls NOT recovered: %s",
                    summary["positive_controls_missing"])


if __name__ == "__main__":
    main()
