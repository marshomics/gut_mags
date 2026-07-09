#!/usr/bin/env python3
"""selection_aggregate.py -- combine per-clade HyPhy selection tables into one,
and summarise the sequence-level-selection driver: how many human-associated
ortholog families show episodic positive selection (BUSTED) and intensified
selection on the human branches (RELAX K>1)."""
import argparse
import glob

import pandas as pd

from hgn_utils import get_logger

log = get_logger("sel-agg")

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--out-prefix", required=True)
a = ap.parse_args()

frames = []
for f in a.inputs:
    try:
        d = pd.read_csv(f, sep="\t")
        if len(d):
            frames.append(d)
    except Exception:
        continue
alld = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
    columns=["clade", "family", "busted_p", "busted_sig", "relax_K", "relax_p", "relax_class"])
alld.to_csv(f"{a.out_prefix}_all.tsv", sep="\t", index=False)

summary = {
    "n_families_tested": int(len(alld)),
    "n_positive_selection_busted": int(alld.get("busted_sig", pd.Series(dtype=bool))
                                        .astype(str).isin(["True", "TRUE", "1"]).sum()),
    "n_intensified_relax": int((alld.get("relax_class") == "intensified").sum()),
    "n_relaxed_relax": int((alld.get("relax_class") == "relaxed").sum()),
    "n_clades": int(alld["clade"].nunique()) if len(alld) else 0,
}
pd.DataFrame([summary]).to_csv(f"{a.out_prefix}_summary.tsv", sep="\t", index=False)
log.info("Selection aggregate: %s", summary)
