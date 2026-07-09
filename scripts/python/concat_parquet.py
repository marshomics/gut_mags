#!/usr/bin/env python3
"""concat_parquet.py -- concatenate parquet shards into one file, streaming row
groups so memory stays bounded for the large KO/eggNOG tables.

Schema safety: the inputs may come from different producers (e.g. the per-layer
prevalence files). We take the INTERSECTION of column names across all inputs, in
the order of the first file, and project every batch onto it. That prevents a
silent crash when one producer adds an extra column, and it never silently drops
a column that all inputs share.
"""
import argparse

import pyarrow as pa
import pyarrow.parquet as pq

ap = argparse.ArgumentParser()
ap.add_argument("--inputs", nargs="+", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

if not a.inputs:
    pq.write_table(pa.table({"genome": [], "layer": [], "feature": [], "count": []}), a.out)
    raise SystemExit(0)

# common columns, ordered as in the first input
first_cols = list(pq.ParquetFile(a.inputs[0]).schema_arrow.names)
common = set(first_cols)
for p in a.inputs[1:]:
    common &= set(pq.ParquetFile(p).schema_arrow.names)
cols = [c for c in first_cols if c in common]
dropped = [c for c in first_cols if c not in common]
if dropped:
    print(f"NOTE: columns not present in every input, dropped: {dropped}")
if not cols:
    raise SystemExit("FATAL: inputs share no common columns")

writer, n = None, 0
for path in a.inputs:
    pf = pq.ParquetFile(path)
    for batch in pf.iter_batches(batch_size=1_000_000, columns=cols):
        tbl = pa.Table.from_batches([batch]).select(cols)
        if writer is None:
            writer = pq.ParquetWriter(a.out, tbl.schema)
        writer.write_table(tbl)
        n += tbl.num_rows
if writer is not None:
    writer.close()
else:
    pq.write_table(pa.table({c: [] for c in cols}), a.out)
print(f"Concatenated {len(a.inputs)} shards -> {a.out} ({n} rows, {len(cols)} cols)")
