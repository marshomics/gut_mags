#!/usr/bin/env python3
"""combine_tsv.py -- concatenate TSV shards (e.g. phyloglm feature chunks)
into one table, keeping a single header. Used to merge parallel chunk outputs."""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
pd.concat([pd.read_csv(f, sep="\t") for f in a.inputs], ignore_index=True) \
  .to_csv(a.out, sep="\t", index=False)
print(f"Combined {len(a.inputs)} shards -> {a.out}")
