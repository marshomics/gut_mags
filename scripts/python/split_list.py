#!/usr/bin/env python3
"""
split_list.py -- split the genome id list (all genomes, and species
representatives) into fixed-size chunk files for array-style parallel
annotation. Deterministic ordering so chunk membership is reproducible.
"""
import argparse
import os

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--samples", required=True)
ap.add_argument("--representatives", required=True)
ap.add_argument("--chunk-size", type=int, required=True)
ap.add_argument("--out-dir", required=True)
a = ap.parse_args()

os.makedirs(a.out_dir, exist_ok=True)
allg = sorted(pd.read_parquet(a.samples)["genome"].astype(str).tolist())
reps = sorted(l.strip() for l in open(a.representatives) if l.strip())


def write_chunks(ids, prefix):
    n = 0
    for i in range(0, len(ids), a.chunk_size):
        with open(os.path.join(a.out_dir, f"{prefix}_{i // a.chunk_size:04d}.txt"), "w") as fh:
            fh.write("\n".join(ids[i:i + a.chunk_size]) + "\n")
        n += 1
    return n


na = write_chunks(allg, "all")
nr = write_chunks(reps, "reps")
print(f"Wrote {na} all-genome chunks and {nr} representative chunks "
      f"({len(allg)} genomes, {len(reps)} reps)")
