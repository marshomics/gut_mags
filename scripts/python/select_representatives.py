#!/usr/bin/env python3
"""
select_representatives.py
-------------------------
Choose one representative genome per species for (a) de-novo annotation that is
too costly to run on all 581k genomes (antiSMASH) and (b) mapping species onto
the GTDB phylogenetic scaffold.

Representative = highest quality_score (Completeness - 5*Contamination), tie-
broken by N50 then genome size, per config species.representative.rank_by.
Ties beyond that are broken by genome id for determinism.

Outputs:
  representatives.tsv      species -> representative genome id (+ its stats)
  representatives.txt      bare list of representative genome ids (for array jobs)
"""
import argparse

import pandas as pd

from hgn_utils import load_config, get_logger, set_global_seed

log = get_logger("reps")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    gid = cfg["inputs"]["columns"]["genome_id"]
    rank_by = cfg["species"]["representative"]["rank_by"]

    df = pd.read_parquet(args.samples)
    # ensure deterministic, documented ordering
    sort_cols = rank_by + [gid]
    ascending = [False] * len(rank_by) + [True]
    df_sorted = df.sort_values(sort_cols, ascending=ascending, kind="mergesort")
    reps = df_sorted.groupby("species", as_index=False).first()

    keep = [gid, "species", "niche", "domain", "quality_score",
            "completeness", "contamination", "n50", "genome_size"]
    keep = [c for c in keep if c in reps.columns]
    reps[keep].to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    reps[gid].to_csv(f"{args.out_prefix}.txt", index=False, header=False)

    log.info("Selected %d representatives (one per species)", len(reps))


if __name__ == "__main__":
    main()
