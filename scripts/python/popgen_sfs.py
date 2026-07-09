#!/usr/bin/env python3
"""
popgen_sfs.py
-------------
Within-niche population genetics from the recombination-masked SNP alignment
(Gubbins filtered_polymorphic_sites.fasta), the direct test for a founder
bottleneck and post-colonisation expansion that mark a recently acquired niche.

Per niche population (outgroup dropped), after subsampling to equal strain number
each bootstrap so the comparison is not driven by sample size:
  S (segregating sites), pi (total pairwise differences), Watterson theta,
  Tajima's D, folded SFS, and the singleton fraction (expansion indicator).
A recently colonised niche should show lower pi/theta and a Tajima's D pushed
negative by an excess of rare variants.

Between niches: Hudson Fst, dxy, net divergence da, and a split estimate (da in
substitutions/site; converted to time if a mutation rate is supplied, else left
in mutational units). Statistics use SNP sites only, which is valid for S, total
pi, theta and Tajima's D; per-bp values are reported when the core-alignment
length is provided.

Formulas are implemented directly (Watterson 1975; Tajima 1989; Hudson Fst) so
each step is auditable. Gaps/N are treated as missing per site.
"""
import argparse
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("popgen")
BASES = set("ACGT")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def a1a2(n):
    i = np.arange(1, n)
    return (1.0 / i).sum(), (1.0 / i ** 2).sum()


def tajima_constants(n):
    a1, a2 = a1a2(n)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n ** 2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / a1 ** 2
    e1 = c1 / a1
    e2 = c2 / (a1 ** 2 + a2)
    return a1, e1, e2


def pop_stats(mat):
    """mat: list of equal-length strings (one per sequence). Returns dict."""
    n = len(mat)
    if n < 2:
        return None
    arr = np.array([list(s) for s in mat])           # n x L
    L = arr.shape[1]
    a1, e1, e2 = tajima_constants(n)
    S = 0
    pair_diff_total = 0.0
    sfs = np.zeros(n // 2 + 1, dtype=float)           # folded
    singletons = 0
    for j in range(L):
        col = arr[:, j]
        keep = np.array([c in BASES for c in col])
        if keep.sum() < 2:
            continue
        c = col[keep]
        alleles, counts = np.unique(c, return_counts=True)
        m = len(c)
        if len(alleles) < 2:
            continue
        S += 1
        # per-site pairwise differences (unbiased over present samples)
        same = sum(cc * (cc - 1) / 2 for cc in counts)
        diff = m * (m - 1) / 2 - same
        pair_diff_total += diff / (m * (m - 1) / 2)    # per-site pi contribution
        minor = counts.min()
        fold = min(minor, m - minor)
        if fold <= n // 2:
            sfs[fold] += 1
        if fold == 1:
            singletons += 1
    pi_total = pair_diff_total                          # summed per-site pi (per-bp units already)
    thetaW = S / a1
    var = e1 * S + e2 * S * (S - 1)
    D = (pi_total - thetaW) / np.sqrt(var) if var > 0 and S > 0 else np.nan
    return {"n": n, "S": S, "pi": pi_total, "thetaW": thetaW, "tajimaD": D,
            "singleton_frac": singletons / S if S else np.nan, "sfs": sfs, "L": L}


def hudson(matX, matY):
    """Hudson Fst, dxy, da between two populations (per-site averages)."""
    aX = np.array([list(s) for s in matX]); aY = np.array([list(s) for s in matY])
    L = aX.shape[1]
    dxy_sum = piX_sum = piY_sum = 0.0; nsite = 0
    for j in range(L):
        cx = [c for c in aX[:, j] if c in BASES]
        cy = [c for c in aY[:, j] if c in BASES]
        if len(cx) < 1 or len(cy) < 1:
            continue
        nsite += 1
        # dxy: prob two seqs (one per pop) differ
        ax, axc = np.unique(cx, return_counts=True)
        ay, ayc = np.unique(cy, return_counts=True)
        px = dict(zip(ax, axc / len(cx))); py = dict(zip(ay, ayc / len(cy)))
        same = sum(px.get(b, 0) * py.get(b, 0) for b in BASES)
        dxy_sum += 1 - same
        mx = len(cx); my = len(cy)
        if mx > 1:                                   # unbiased within-pop pi per site
            piX_sum += (mx / (mx - 1)) * (1 - sum((c / mx) ** 2 for c in axc))
        if my > 1:
            piY_sum += (my / (my - 1)) * (1 - sum((c / my) ** 2 for c in ayc))
    dxy = dxy_sum / nsite if nsite else np.nan
    piX = piX_sum / nsite if nsite else np.nan
    piY = piY_sum / nsite if nsite else np.nan
    Hw = (piX + piY) / 2
    da = dxy - Hw
    fst = (dxy - Hw) / dxy if dxy and dxy > 0 else np.nan
    return dxy, da, fst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--niche-map", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--core-length", type=int, default=0)
    args = ap.parse_args()

    cfg = load_config(args.config)
    pc = cfg["transition"]["popgen"]
    dc = cfg["transition"]["demography"]
    os.makedirs(args.out_dir, exist_ok=True)

    seqs = read_fasta(args.alignment)
    nm = pd.read_csv(args.niche_map, sep="\t")
    nm = nm[nm["role"] == "focal"]
    pops = {}
    for niche, d in nm.groupby("niche"):
        ids = [g for g in d["genome"] if g in seqs]
        if len(ids) >= 2:
            pops[niche] = ids
    niches = sorted(pops)
    if len(niches) < 2:
        log.warning("Fewer than 2 niches with sequences; writing empty outputs.")
        pd.DataFrame().to_csv(f"{args.out_dir}/popgen_diversity.tsv", sep="\t")
        pd.DataFrame().to_csv(f"{args.out_dir}/divergence.tsv", sep="\t")
        return

    n_min = min(len(v) for v in pops.values())
    B = pc["bootstrap"]
    L = args.core_length

    div_rows = []
    sfs_store = {}
    for niche in niches:
        ids = pops[niche]
        boots = {k: [] for k in ["S", "pi", "thetaW", "tajimaD", "singleton_frac"]}
        sfs_acc = None
        for b in range(B):
            rng = np.random.default_rng(derive_seed(cfg["seed"], "popgen", niche, b))
            pick = rng.choice(ids, size=n_min, replace=len(ids) < n_min)
            st = pop_stats([seqs[g] for g in pick])
            if st is None:
                continue
            for k in boots:
                boots[k].append(st[k])
            s = st["sfs"]
            sfs_acc = s if sfs_acc is None else sfs_acc + s
        row = {"niche": niche, "n_pop": len(ids), "n_subsample": n_min}
        for k in boots:
            v = np.array(boots[k], float)
            row[f"{k}_mean"] = np.nanmean(v)
            row[f"{k}_lo"] = np.nanpercentile(v, 2.5)
            row[f"{k}_hi"] = np.nanpercentile(v, 97.5)
        if L:
            row["pi_per_bp"] = row["pi_mean"] / L
            row["thetaW_per_bp"] = row["thetaW_mean"] / L
        div_rows.append(row)
        sfs_store[niche] = (sfs_acc / B) if sfs_acc is not None else None

    pd.DataFrame(div_rows).to_csv(f"{args.out_dir}/popgen_diversity.tsv", sep="\t", index=False)
    for niche, s in sfs_store.items():
        if s is not None:
            pd.DataFrame({"minor_allele_count": np.arange(len(s)), "mean_sites": s}).to_csv(
                f"{args.out_dir}/sfs_{niche}.tsv", sep="\t", index=False)

    # pairwise divergence + split estimate
    mu = dc.get("mutation_rate_per_site_per_year")
    drows = []
    for i in range(len(niches)):
        for j in range(i + 1, len(niches)):
            X, Y = niches[i], niches[j]
            dxy, da, fst = hudson([seqs[g] for g in pops[X]], [seqs[g] for g in pops[Y]])
            rec = {"pair": f"{X}|{Y}", "dxy": dxy, "da": da, "Fst": fst}
            if mu:
                rec["split_years"] = da / (2 * mu) if da and da > 0 else np.nan
            else:
                rec["split_subs_per_site"] = da
            drows.append(rec)
    pd.DataFrame(drows).to_csv(f"{args.out_dir}/divergence.tsv", sep="\t", index=False)
    log.info("Popgen: niches=%s, subsample n=%d; divergence pairs=%d",
             niches, n_min, len(drows))


if __name__ == "__main__":
    main()
