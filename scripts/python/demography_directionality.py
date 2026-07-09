#!/usr/bin/env python3
"""
demography_directionality.py
----------------------------
Which niche population is derived (recently founded) and which is the source,
inferred from the joint polymorphism pattern; plus the 2D folded SFS for
optional downstream coalescent modelling.

A population founded from another carries a SUBSET of the source's standing
variation: colonisation samples a few lineages, so ancestral alleles are lost,
few new ones have had time to arise, and diversity is lower. Three consequences
are measured, each at EQUAL SAMPLE SIZE (bootstrap subsampling to the smaller
population, because segregating-site counts and private-allele counts both grow
with n and would otherwise simply track sample size):

  L1 nestedness   fraction of a population's segregating sites at which its
                  FULL allele set is contained in the other population. The
                  derived population is the more nested one.
  L2 private alleles  number of alleles seen in one population and not the
                  other. The derived population has fewer (its variation is
                  inherited; the source retains alleles lost at the founding).
  L3 diversity    nucleotide diversity pi. The derived population has less.

Uncertainty and significance are separate things here, and are computed
separately. A BOOTSTRAP over subsamples of the two fixed populations gives the
confidence interval of each estimate; it says nothing about the null, because
resampling within fixed populations reproduces whatever asymmetry those samples
happen to contain. Significance therefore comes from a LABEL PERMUTATION: niche
labels are shuffled across the pooled genomes, the statistic recomputed, and the
observed value compared with that distribution. All three statistics are
antisymmetric under swapping the two labels, so the permutation null is centred
on zero by construction and the two-sided empirical p is well defined.

The call is made by L1 when its permutation p is significant, otherwise by L2
and L3 when they agree; when nothing resolves, derived_call is "unresolved"
rather than a coin flip. n_lines_agree records the corroboration, and
transition_verdict.py consumes derived_call as one of five evidence lines.

Definition note: nestedness is set containment, {alleles in Y at site j} subset
of {alleles in X at site j}. An earlier formulation asked only whether Y's MINOR
allele occurred in X; which allele is minor depends on frequencies rather than
on ancestry, so that statistic was not a nestedness measure.

Polarisation is not required for any of the above, so no outgroup is assumed and
no derived-allele mis-polarisation can bias the call.

The 2D folded SFS is written per pair for optional fastsimcoal2 / dadi modelling
(see resources/fastsimcoal2_recipe.md); this script does not run those tools, so
no untested demographic templates are shipped.

Outputs: directionality.tsv, joint_sfs_<X>_<Y>.tsv
"""
import argparse
import itertools
import os

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, derive_seed

log = get_logger("demog-dir")
BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def encode(seqs):
    """(n x L) int8 array; ACGT -> 0..3, anything else (N, -, IUPAC) -> -1."""
    arr = np.full((len(seqs), len(seqs[0])), -1, dtype=np.int8)
    for i, s in enumerate(seqs):
        v = np.frombuffer(s.encode(), dtype=np.uint8)
        for b, k in CODE.items():
            arr[i, v == ord(b)] = k
    return arr


def allele_presence(a):
    """(L x 4) bool: which alleles are observed in this population at each site,
    and (L,) int: number of called sequences per site."""
    counts = np.zeros((a.shape[1], 4), dtype=np.int32)
    for k in range(4):
        counts[:, k] = (a == k).sum(axis=0)
    return counts > 0, counts, counts.sum(axis=1)


def pi_from_counts(counts, n_called):
    """Mean per-site nucleotide diversity, unbiased: (n^2 - sum c^2)/(n(n-1))."""
    ok = n_called >= 2
    if not ok.any():
        return np.nan
    n = n_called[ok].astype(float)
    ss = (counts[ok].astype(float) ** 2).sum(axis=1)
    return float((((n ** 2 - ss) / (n * (n - 1)))).mean())


def pair_directionality(aX, aY):
    """One replicate: X and Y are (n x L) encoded arrays with equal n."""
    pX, cX, nX = allele_presence(aX)
    pY, cY, nY = allele_presence(aY)
    ok = (nX >= 1) & (nY >= 1)                       # both populations observed
    pX, pY, cX, cY = pX[ok], pY[ok], cX[ok], cY[ok]
    nXo, nYo = nX[ok], nY[ok]

    kX, kY = pX.sum(axis=1), pY.sum(axis=1)
    polyX, polyY = kX > 1, kY > 1                    # segregating in that population

    # set containment: every allele of Y at this site is also seen in X
    Y_in_X = ~(pY & ~pX).any(axis=1)
    X_in_Y = ~(pX & ~pY).any(axis=1)

    nestY_in_X = float((Y_in_X & polyY).sum() / polyY.sum()) if polyY.sum() else np.nan
    nestX_in_Y = float((X_in_Y & polyX).sum() / polyX.sum()) if polyX.sum() else np.nan

    return {
        "S_X": int(polyX.sum()), "S_Y": int(polyY.sum()),
        "private_X": int((polyX & ~polyY).sum()),      # privately segregating sites
        "private_Y": int((polyY & ~polyX).sum()),
        "private_alleles_X": int((pX & ~pY).sum()),    # alleles absent from the other
        "private_alleles_Y": int((pY & ~pX).sum()),
        "shared": int((polyX & polyY).sum()),
        "polyX": int(polyX.sum()), "polyY": int(polyY.sum()),
        "nestY_in_X": nestY_in_X, "nestX_in_Y": nestX_in_Y,
        "pi_X": pi_from_counts(cX, nXo), "pi_Y": pi_from_counts(cY, nYo),
    }


def folded_joint_sfs(aX, aY):
    """2D folded SFS over biallelic sites called in both populations."""
    nX, nY = aX.shape[0], aY.shape[0]
    M = np.zeros((nX + 1, nY + 1))
    pX, cX, nXc = allele_presence(aX)
    pY, cY, nYc = allele_presence(aY)
    biallelic = ((pX | pY).sum(axis=1) == 2) & (nXc == nX) & (nYc == nY)
    for j in np.flatnonzero(biallelic):
        alleles = np.flatnonzero(pX[j] | pY[j])
        a0 = alleles[0]                       # lexicographically first of ACGT
        i, k = int(cX[j, a0]), int(cY[j, a0])
        if i + k > (nX + nY) / 2:             # fold by the global minor allele
            i, k = nX - i, nY - k
        M[i, k] += 1
    return M


def boot_ci(diffs, ci=0.95):
    """Percentile CI of an estimate from its bootstrap replicates."""
    d = np.asarray([x for x in diffs if np.isfinite(x)], float)
    if d.size < 2:
        return np.nan, np.nan
    lo, hi = np.percentile(d, [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100])
    return float(lo), float(hi)


def perm_p(obs, null):
    """Two-sided empirical p from a permutation null (Phipson & Smyth: add one)."""
    n = np.asarray([x for x in null if np.isfinite(x)], float)
    if not np.isfinite(obs) or n.size < 10:
        return np.nan
    return float((1 + (np.abs(n) >= abs(obs)).sum()) / (n.size + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--niche-map", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    B = cfg["transition"]["popgen"]["bootstrap"]
    min_poly = int(cfg["transition"]["popgen"].get("min_segregating_sites", 10))
    alpha = float(cfg["stats"]["fdr_alpha"])
    os.makedirs(args.out_dir, exist_ok=True)

    seqs = read_fasta(args.alignment)
    nm = pd.read_csv(args.niche_map, sep="\t")
    nm = nm[nm["role"] == "focal"]
    pops = {n: [g for g in d["genome"] if g in seqs]
            for n, d in nm.groupby("niche") if (d["genome"].isin(seqs)).sum() >= 2}
    niches = sorted(pops)
    if len(niches) < 2:
        pd.DataFrame().to_csv(f"{args.out_dir}/directionality.tsv", sep="\t", index=False)
        log.warning("Fewer than two niches with >=2 genomes; no directionality."); return
    n_min = min(len(v) for v in pops.values())
    enc = {n: encode([seqs[g] for g in gs]) for n, gs in pops.items()}

    rows = []
    for X, Y in itertools.combinations(niches, 2):
        # ---- bootstrap: estimates + CIs, subsampling each population to n_min --
        boots = []
        for b in range(B):
            rng = np.random.default_rng(derive_seed(cfg["seed"], "demog", X, Y, b))
            ix = rng.choice(len(pops[X]), size=n_min, replace=len(pops[X]) < n_min)
            iy = rng.choice(len(pops[Y]), size=n_min, replace=len(pops[Y]) < n_min)
            boots.append(pair_directionality(enc[X][ix], enc[Y][iy]))
        bd = pd.DataFrame(boots)

        nyx, nxy = bd["nestY_in_X"].mean(), bd["nestX_in_Y"].mean()
        pax, pay = bd["private_alleles_X"].mean(), bd["private_alleles_Y"].mean()
        pix, piy = bd["pi_X"].mean(), bd["pi_Y"].mean()
        obs_nest = nyx - nxy                 # >0 => Y more nested => Y derived
        obs_priv = pax - pay                 # >0 => Y has fewer private => Y derived
        obs_pi = pix - piy                   # >0 => Y less diverse   => Y derived
        lo, hi = boot_ci(bd["nestY_in_X"] - bd["nestX_in_Y"])
        enough = (bd["polyY"].median() >= min_poly) and (bd["polyX"].median() >= min_poly)

        # ---- permutation null: shuffle niche labels over the pooled genomes ----
        pool = np.vstack([enc[X], enc[Y]])
        nX = len(pops[X])
        null = {"nest": [], "priv": [], "pi": []}
        for b in range(B):
            rng = np.random.default_rng(derive_seed(cfg["seed"], "demog_perm", X, Y, b))
            perm = rng.permutation(pool.shape[0])
            px, py = pool[perm[:nX]], pool[perm[nX:]]
            ix = rng.choice(px.shape[0], size=n_min, replace=px.shape[0] < n_min)
            iy = rng.choice(py.shape[0], size=n_min, replace=py.shape[0] < n_min)
            s = pair_directionality(px[ix], py[iy])
            null["nest"].append(s["nestY_in_X"] - s["nestX_in_Y"])
            null["priv"].append(s["private_alleles_X"] - s["private_alleles_Y"])
            null["pi"].append(s["pi_X"] - s["pi_Y"])
        p_nest = perm_p(obs_nest, null["nest"])
        p_priv = perm_p(obs_priv, null["priv"])
        p_pi = perm_p(obs_pi, null["pi"])

        # L1 nestedness (primary): the more-nested population is derived
        nest_call = (Y if obs_nest > 0 else X) if (enough and np.isfinite(p_nest)
                                                   and p_nest < alpha and obs_nest != 0) else None
        # L2 private alleles: the derived population has fewer
        priv_call = (Y if obs_priv > 0 else X) if (np.isfinite(p_priv)
                                                   and p_priv < alpha and obs_priv != 0) else None
        # L3 diversity: the derived population has lower pi
        pi_call = (Y if obs_pi > 0 else X) if (np.isfinite(p_pi)
                                               and p_pi < alpha and obs_pi != 0) else None

        if nest_call is not None:
            derived, status = nest_call, "resolved_nestedness"
        elif priv_call is not None and priv_call == pi_call:
            derived, status = priv_call, "resolved_private_and_pi"
        else:
            derived, status = None, "unresolved"
        calls = [c for c in (nest_call, priv_call, pi_call) if c is not None]
        agree = sum(c == derived for c in calls) if derived else 0

        rows.append({
            "pair": f"{X}|{Y}", "X": X, "Y": Y,
            "source_call": (X if derived == Y else Y) if derived else "unresolved",
            "derived_call": derived if derived else "unresolved",
            "call_status": status, "n_lines_agree": agree, "n_lines_evaluated": len(calls),
            "nest_Y_in_X": round(nyx, 4), "nest_X_in_Y": round(nxy, 4),
            "nest_diff": round(obs_nest, 4),
            "nest_diff_ci_lo": round(lo, 4) if np.isfinite(lo) else np.nan,
            "nest_diff_ci_hi": round(hi, 4) if np.isfinite(hi) else np.nan,
            "nest_p": round(p_nest, 4) if np.isfinite(p_nest) else np.nan,
            "nestedness_call": nest_call or "unresolved",
            "private_alleles_X": round(pax, 1), "private_alleles_Y": round(pay, 1),
            "private_alleles_p": round(p_priv, 4) if np.isfinite(p_priv) else np.nan,
            "private_alleles_call": priv_call or "unresolved",
            "pi_X": round(pix, 5) if np.isfinite(pix) else np.nan,
            "pi_Y": round(piy, 5) if np.isfinite(piy) else np.nan,
            "pi_ratio_YoverX": round(piy / pix, 4) if np.isfinite(pix) and pix > 0 else np.nan,
            "pi_p": round(p_pi, 4) if np.isfinite(p_pi) else np.nan,
            "pi_call": pi_call or "unresolved",
            "private_X": int(bd["private_X"].mean()), "private_Y": int(bd["private_Y"].mean()),
            "shared": int(bd["shared"].mean()),
            "S_X": int(bd["polyX"].mean()), "S_Y": int(bd["polyY"].mean()),
            "n_subsample": n_min, "n_bootstrap": B, "n_permutations": B,
            "sufficient_polymorphism": bool(enough)})
        # 2D folded SFS (one full-data draw at equal n)
        M = folded_joint_sfs(enc[X][:n_min], enc[Y][:n_min])
        pd.DataFrame(M).to_csv(f"{args.out_dir}/joint_sfs_{X}_{Y}.tsv", sep="\t")

    pd.DataFrame(rows).to_csv(f"{args.out_dir}/directionality.tsv", sep="\t", index=False)
    log.info("Directionality (n=%d per population, %d bootstraps, %d permutations): %s",
             n_min, B, B, [(r["pair"], r["derived_call"], r["call_status"]) for r in rows])


if __name__ == "__main__":
    main()
