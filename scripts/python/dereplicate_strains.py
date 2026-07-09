#!/usr/bin/env python3
"""
dereplicate_strains.py
----------------------
Collapse near-identical strains so clonal/epidemic oversampling does not fake low
diversity or a spurious clade. Critically, dereplication is done WITHIN each
niche only: near-identical strains that span two niches are the signature of very
recent transfer and must be kept, so they are never collapsed across niches.

Single-linkage clustering on Mash distances at `max_distance` (default ~99.99%
ANI); the highest-quality genome represents each cluster, and cluster sizes are
recorded (so the original sampling is auditable and can re-weight if needed). The
outgroup genome is always retained.

Inputs:
  --candidates  candidate_genomes.tsv (genome, niche, role)
  --mash-dist   mash dist output: g1 <tab> g2 <tab> dist <tab> p <tab> shared
  --samples     samples.parquet (for quality_score)
Outputs:
  --out         dereplicated_genomes.tsv (genome, niche, role)
  --clusters    clusters.tsv (genome, niche, cluster_id, cluster_size, kept)
"""
import argparse

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("derep")


def components(nodes, edges):
    """Union-find connected components."""
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in edges:
        parent[find(a)] = find(b)
    comp = {}
    for n in nodes:
        comp.setdefault(find(n), []).append(n)
    return list(comp.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--mash-dist", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--clusters", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    dc = cfg["transition"]["dereplicate"]
    thr = dc["max_distance"]

    cand = pd.read_csv(args.candidates, sep="\t")
    qual = (pd.read_parquet(args.samples)[["genome", "quality_score"]]
            .set_index("genome")["quality_score"].to_dict())
    niche_of = dict(zip(cand["genome"], cand["niche"]))
    role_of = dict(zip(cand["genome"], cand["role"]))

    if not dc["enabled"]:
        cand.to_csv(args.out, sep="\t", index=False)
        cand.assign(cluster_id=range(len(cand)), cluster_size=1, kept=True).to_csv(
            args.clusters, sep="\t", index=False); return

    md = pd.read_csv(args.mash_dist, sep="\t", header=None,
                     names=["g1", "g2", "dist", "p", "shared"])
    # normalise ids to basenames if paths were used
    for c in ("g1", "g2"):
        md[c] = md[c].astype(str).map(lambda x: x.split("/")[-1]
                                      .replace(".fna", "").replace(".fasta", ""))

    cluster_rows, kept = [], []
    cid = 0
    for niche, d in cand.groupby("niche"):
        gs = set(d["genome"])
        if niche == "outgroup":
            for g in gs:
                kept.append(g); cluster_rows.append((g, niche, cid, 1, True)); cid += 1
            continue
        # within-niche edges below threshold
        sub = md[(md["g1"].isin(gs)) & (md["g2"].isin(gs)) & (md["dist"] <= thr)]
        edges = list(zip(sub["g1"], sub["g2"]))
        for comp in components(gs, edges):
            rep = max(comp, key=lambda g: qual.get(g, 0))
            kept.append(rep)
            for g in comp:
                cluster_rows.append((g, niche, cid, len(comp), g == rep))
            cid += 1

    out = pd.DataFrame([{"genome": g, "niche": niche_of[g], "role": role_of[g]}
                        for g in kept])
    out.to_csv(args.out, sep="\t", index=False)
    pd.DataFrame(cluster_rows, columns=["genome", "niche", "cluster_id",
                                        "cluster_size", "kept"]).to_csv(
        args.clusters, sep="\t", index=False)
    n_in = (cand["role"] == "focal").sum()
    n_out = (out["role"] == "focal").sum()
    log.info("Dereplicated focal genomes %d -> %d (within-niche, dist<=%.5f)",
             n_in, n_out, thr)


if __name__ == "__main__":
    main()
