#!/usr/bin/env python3
"""
overlap_nullmodel.py
--------------------
Is the cross-niche species overlap lower than expected by chance? "Niche-
specific with minor overlap" is only meaningful relative to a null. The species x
niche occupancy matrix is randomised while preserving margins, and the observed
shared-species counts are compared to the null.

Null models (config taxonomy.overlap_null.model):
  * curveball  : preserves BOTH margins exactly (each species' number of occupied
                 niches and each niche's richness). The conservative, standard
                 choice (Strona et al. 2014).
  * fixed_degree: preserves each species' degree exactly and niche richness in
                 expectation (each species' niches drawn without replacement with
                 probability proportional to niche size). Faster.

Reports, for each niche pair and the triple intersection: observed overlap,
null mean/sd, standardised effect size (SES) and a two-sided empirical p. A
strongly negative SES with small p is the statistical statement of niche
specificity.

Output: overlap_nullmodel.tsv, overlap_observed.tsv
"""
import argparse
import itertools

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("overlapnull")


def overlap_stats(M, niches):
    """M: bool array species x niche. Return dict of pair/triple overlaps."""
    out = {}
    K = M.shape[1]
    for i, j in itertools.combinations(range(K), 2):
        out[f"{niches[i]}&{niches[j]}"] = int((M[:, i] & M[:, j]).sum())
    if K >= 3:
        out["&".join(niches)] = int(M.all(axis=1).sum())
    return out


def curveball(M, n_swaps, rng):
    """In-place fixed-fixed randomisation via curveball trades (margins preserved)."""
    rows = [set(np.where(r)[0]) for r in M]
    n = len(rows)
    done = 0
    while done < n_swaps:
        a, b = rng.integers(0, n), rng.integers(0, n)
        if a == b:
            continue
        ra, rb = rows[a], rows[b]
        only_a = list(ra - rb)
        only_b = list(rb - ra)
        if not only_a or not only_b:
            continue
        # swap one differing membership each way (preserves both row and col sums)
        x = only_a[rng.integers(0, len(only_a))]
        y = only_b[rng.integers(0, len(only_b))]
        ra.discard(x); ra.add(y)
        rb.discard(y); rb.add(x)
        done += 1
    out = np.zeros_like(M)
    for i, s in enumerate(rows):
        for c in s:
            out[i, c] = True
    return out


def fixed_degree(M, rng, niche_sizes):
    """Preserve each species' degree; draw its niches without replacement weighted
    by niche richness."""
    K = M.shape[1]
    deg = M.sum(axis=1)
    p = niche_sizes / niche_sizes.sum()
    out = np.zeros_like(M)
    for i, d in enumerate(deg):
        choice = rng.choice(K, size=int(d), replace=False, p=p)
        out[i, choice] = True
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    niches = cfg["inputs"]["niche_levels"]
    onc = cfg["taxonomy"]["overlap_null"]
    N = onc["iterations"]
    model = onc["model"]

    import os
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    samples = pd.read_parquet(args.samples)
    occ = (samples[["species", "niche"]].drop_duplicates()
           .assign(v=1).pivot_table(index="species", columns="niche", values="v",
                                    fill_value=0))
    occ = occ.reindex(columns=niches, fill_value=0)
    M = occ.to_numpy().astype(bool)
    niche_sizes = M.sum(axis=0).astype(float)

    obs = overlap_stats(M, niches)
    pd.Series(obs, name="observed").rename_axis("set").reset_index().to_csv(
        f"{args.out_prefix}_observed.tsv", sep="\t", index=False)

    keys = list(obs.keys())
    null = {k: np.empty(N) for k in keys}
    rng = np.random.default_rng(cfg["seed"])
    n_swaps = 5 * int(M.sum())
    # config uses "fixed_fixed" (both margins preserved) as the canonical name;
    # "curveball" is the algorithm that implements it. Anything else falls back to
    # the weaker fixed-degree null, which is stated explicitly in the output.
    use_curveball = model in ("curveball", "fixed_fixed")
    if not use_curveball:
        log.warning("overlap null model '%s' preserves species degree only "
                    "(niche richness in expectation), not both margins.", model)
    for it in range(N):
        if use_curveball:
            Mr = curveball(M, n_swaps, rng)
        else:
            Mr = fixed_degree(M, rng, niche_sizes)
        st = overlap_stats(Mr, niches)
        for k in keys:
            null[k][it] = st[k]
        if (it + 1) % max(1, N // 10) == 0:
            log.info("null %d/%d", it + 1, N)

    rows = []
    for k in keys:
        nd = null[k]; o = obs[k]
        mu, sd = nd.mean(), nd.std(ddof=1)
        ses = (o - mu) / sd if sd > 0 else np.nan
        # two-sided empirical p (Davison-Hinkley: (1 + #as-extreme) / (N + 1))
        p = min(1.0, 2 * (1 + min((nd <= o).sum(), (nd >= o).sum())) / (N + 1))
        rows.append({"set": k, "observed": o, "null_mean": round(mu, 2),
                     "null_sd": round(sd, 2), "SES": round(ses, 3) if ses == ses else np.nan,
                     "p_two_sided": p,
                     "null_model": "curveball(fixed-fixed)" if use_curveball else "fixed-degree",
                     "interpretation": ("undefined (no variance in null)" if not (sd > 0)
                                        else "less overlap than chance" if ses < 0
                                        else "more overlap than chance")})
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}_nullmodel.tsv", sep="\t", index=False)
    log.info("Overlap null (%s) done: %s", model,
             {r["set"]: r["SES"] for r in rows})


if __name__ == "__main__":
    main()
