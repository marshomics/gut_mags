#!/usr/bin/env python3
"""combine_signatures.py -- concatenate per-layer consensus signature tables,
adding a `layer` column, into one table per contrast (used by the tree,
signal, ancestral and heatmap figures)."""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--layers", nargs="+", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

frames = []
for path, layer in zip(a.inputs, a.layers):
    d = pd.read_csv(path, sep="\t")
    d["layer"] = layer
    frames.append(d)
out = pd.concat(frames, ignore_index=True)
out.to_csv(a.out, sep="\t", index=False)
print(f"Combined {len(frames)} layers -> {a.out} ({len(out)} rows)")
