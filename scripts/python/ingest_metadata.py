#!/usr/bin/env python3
"""
ingest_metadata.py
------------------
Parse the master metadata table, validate it, apply explicit QC, and emit a
clean, typed sample sheet plus a QC report. This is the single gate every
downstream step passes through, so the analysis set is defined here and
nowhere else.

Decisions made defensible:
  * QC thresholds come from config (qc:), not from upstream filtering.
  * The `qc` column is independently verified to equal Completeness-5*Contam.
  * Rows failing QC are written to a separate file with the reason, never
    silently dropped.
  * Niche values are validated against the declared niche_levels.

Outputs:
  samples.parquet        clean per-genome table (typed)
  samples.tsv            same, human-readable
  qc_failures.tsv        excluded genomes + reason
  qc_report.json         counts before/after each filter, per niche
"""
import argparse
import json

import numpy as np
import pandas as pd

from hgn_utils import (load_config, get_logger, set_global_seed, numeric,
                       is_placeholder_species, provenance_stamp)

log = get_logger("ingest")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    cols = cfg["inputs"]["columns"]
    qc = cfg["qc"]

    log.info("Reading metadata: %s", cfg["inputs"]["metadata"])
    df = pd.read_csv(cfg["inputs"]["metadata"], sep="\t",
                     dtype=str, keep_default_na=False)
    n0 = len(df)
    log.info("Loaded %d genomes, %d columns", n0, df.shape[1])

    # --- column presence check -------------------------------------------------
    required = [cols["genome_id"], cols["niche"], cols["species"],
                cols["quality_completeness"], cols["quality_contamination"]]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"FATAL: required columns absent: {missing}")

    # --- type coercion ---------------------------------------------------------
    numcols = [cols["quality_completeness"], cols["quality_contamination"],
               cols["quality_score"], cols["gc"], cols["cds_number"],
               cols["genome_size"], cols["genome_size_corrected"],
               cols["n_contigs"], cols["n50"]]
    df = numeric(df, [c for c in numcols if c in df.columns])

    # --- niche validation ------------------------------------------------------
    valid_niches = set(cfg["inputs"]["niche_levels"])
    bad_niche = ~df[cols["niche"]].isin(valid_niches)
    if bad_niche.any():
        log.warning("%d rows have unexpected niche values (will fail QC): %s",
                    int(bad_niche.sum()),
                    sorted(df.loc[bad_niche, cols["niche"]].unique())[:10])

    # --- verify qc == Completeness - 5*Contamination ---------------------------
    if cols["quality_score"] in df.columns:
        expected = df[cols["quality_completeness"]] - 5.0 * df[cols["quality_contamination"]]
        delta = (df[cols["quality_score"]] - expected).abs()
        n_bad = int((delta > 0.05).sum())
        log.info("qc formula check: %d / %d rows deviate >0.05 from "
                 "Completeness-5*Contamination", n_bad, n0)
    else:
        df[cols["quality_score"]] = (df[cols["quality_completeness"]]
                                     - 5.0 * df[cols["quality_contamination"]])
        log.info("qc column absent; computed it as Completeness-5*Contamination")

    # --- placeholder species flag ---------------------------------------------
    df["species_is_placeholder"] = df[cols["species"]].map(is_placeholder_species)

    # --- QC filters, recorded one at a time ------------------------------------
    reasons = pd.Series("", index=df.index)

    def fail(mask, label):
        nonlocal reasons
        newly = mask & (reasons == "")
        reasons.loc[newly] = label

    fail(~df[cols["niche"]].isin(valid_niches), "invalid_niche")
    # optional host exclusion (sensitivity: drop the dominant host, e.g. Mouse)
    drop_hosts = cfg.get("filters", {}).get("drop_hosts", [])
    if drop_hosts and cols["host_common_name"] in df.columns:
        fail(df[cols["host_common_name"]].isin(drop_hosts), "dropped_host")
    fail(df[cols["quality_completeness"]].isna() |
         (df[cols["quality_completeness"]] < qc["min_completeness"]),
         f"completeness<{qc['min_completeness']}")
    fail(df[cols["quality_contamination"]].isna() |
         (df[cols["quality_contamination"]] > qc["max_contamination"]),
         f"contamination>{qc['max_contamination']}")
    fail(df[cols["quality_score"]] < qc["min_quality_score"],
         f"qc<{qc['min_quality_score']}")
    if qc.get("min_n50") and cols["n50"] in df.columns:
        fail(df[cols["n50"]].isna() | (df[cols["n50"]] < qc["min_n50"]),
             f"n50<{qc['min_n50']}")
    if qc.get("max_contigs") and cols["n_contigs"] in df.columns:
        fail(df[cols["n_contigs"]] > qc["max_contigs"],
             f"contigs>{qc['max_contigs']}")
    if qc.get("drop_missing_genome_stats"):
        fail(df[cols["genome_size"]].isna() | (df[cols["genome_size"]] <= 0),
             "missing_genome_size")

    passed = reasons == ""
    clean = df.loc[passed].copy()
    failed = df.loc[~passed, [cols["genome_id"], cols["niche"], cols["species"]]].copy()
    failed["reason"] = reasons.loc[~passed].values

    # --- derived convenience columns ------------------------------------------
    clean["niche"] = clean[cols["niche"]]
    clean["species"] = clean[cols["species"]]
    clean["completeness"] = clean[cols["quality_completeness"]]
    clean["contamination"] = clean[cols["quality_contamination"]]
    clean["quality_score"] = clean[cols["quality_score"]]
    clean["genome_size"] = clean[cols["genome_size"]]
    clean["log10_genome_size"] = np.log10(clean["genome_size"].clip(lower=1))
    clean["gc"] = clean[cols["gc"]]
    clean["cds_number"] = clean[cols["cds_number"]]
    clean["n50"] = clean[cols["n50"]]
    clean["domain"] = clean[cols["gtdb_ranks"][0]]
    clean["host_common_name"] = clean.get(cols["host_common_name"], "")
    # population label (Western vs non-Western), human genomes only; normalised.
    # Tolerant of the column being absent (the metadata may not have it yet).
    pop_col = cols.get("western_nonwestern")
    pcfg = cfg.get("population", {})
    clean["population"] = ""
    if pop_col and pop_col in clean.columns:
        vmap = {}
        for canon, variants in (pcfg.get("values", {}) or {}).items():
            for v in variants:
                vmap[str(v).strip().lower()] = canon
        raw = clean[pop_col].astype(str).str.strip().str.lower()
        norm = raw.map(lambda x: vmap.get(x, "")).fillna("")
        # only human genomes carry a population label
        clean["population"] = norm.where(clean["niche"] == "human", "")
        log.info("population labels (human): %s",
                 clean.loc[clean["niche"] == "human", "population"].value_counts().to_dict())
    else:
        log.info("western_nonwestern column not present; population analyses will be inert.")
    # HQ flag for the sensitivity subset
    clean["is_hq"] = ((clean["completeness"] >= qc["hq_completeness"]) &
                      (clean["contamination"] <= qc["hq_contamination"]))

    # --- write -----------------------------------------------------------------
    clean.to_parquet(f"{args.out_prefix}.parquet", index=False)
    clean.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    failed.to_csv(f"{args.out_prefix}_qc_failures.tsv", sep="\t", index=False)

    report = {
        "provenance": provenance_stamp(cfg, {"script": "ingest_metadata"}),
        "n_input": int(n0),
        "n_passed": int(passed.sum()),
        "n_failed": int((~passed).sum()),
        "fail_reasons": failed["reason"].value_counts().to_dict(),
        "passed_per_niche": clean["niche"].value_counts().to_dict(),
        "passed_per_domain": clean["domain"].value_counts().to_dict(),
        "hq_per_niche": clean.loc[clean["is_hq"], "niche"].value_counts().to_dict(),
        "placeholder_species_genomes": int(clean["species_is_placeholder"].sum()),
        "qc_thresholds": qc,
    }
    with open(f"{args.out_prefix}_qc_report.json", "w") as fh:
        json.dump(report, fh, indent=2, default=str)

    log.info("Passed QC: %d / %d genomes", passed.sum(), n0)
    log.info("Per niche: %s", report["passed_per_niche"])


if __name__ == "__main__":
    main()
