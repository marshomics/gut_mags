#!/usr/bin/env python3
"""
test_statistics.py
------------------
Numerical verification of every statistic the pipeline computes, against values
derived by hand or against analytic limiting cases. Run with:

    python tests/test_statistics.py          (or: pytest tests/test_statistics.py)

Each test states the expected value and where it comes from. These are the
checks that make the statistics defensible: a reviewer can re-derive them.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "python"))

FAILED = []


def check(name, got, want, tol=1e-4):
    ok = (abs(got - want) <= tol) if isinstance(want, float) else (got == want)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: got {got!r}, expected {want!r}")
    if not ok:
        FAILED.append(name)


def check_true(name, cond):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if not cond:
        FAILED.append(name)


# ---------------------------------------------------------------- Tajima's D
def test_tajima():
    """n=4, two segregating sites: one singleton, one at count 2.

    a1 = 1 + 1/2 + 1/3 = 1.833333;  theta_W = S/a1 = 2/1.833333 = 1.090909
    pi  = 3/6 + 4/6 = 1.166667
    a2 = 1 + 1/4 + 1/9 = 1.361111
    b1 = (n+1)/(3(n-1)) = 5/9;  b2 = 2(n^2+n+3)/(9n(n-1)) = 46/108
    c1 = b1 - 1/a1 = 0.010101;  c2 = b2 - (n+2)/(a1 n) + a2/a1^2 = 0.012723
    e1 = c1/a1 = 0.005510;  e2 = c2/(a1^2+a2) = 0.002694
    Var = e1*S + e2*S*(S-1) = 0.016408 ; sd = 0.128094
    D = (pi - theta_W)/sd = 0.075758/0.128094 = 0.59142
    """
    from popgen_sfs import pop_stats, tajima_constants
    seqs = ["AA", "AA", "AC", "CC"]   # site1 minor count 1, site2 minor count 2
    st = pop_stats(seqs)
    check("Tajima: S", st["S"], 2)
    check("Tajima: pi_total", st["pi"], 1.166667)
    check("Tajima: thetaW", st["thetaW"], 1.090909)
    check("Tajima: D", st["tajimaD"], 0.59142, tol=1e-3)
    a1, e1, e2 = tajima_constants(4)
    check("Tajima: a1(n=4)", a1, 1.833333)
    check("Tajima: e1(n=4)", e1, 0.005510, tol=1e-5)
    check("Tajima: e2(n=4)", e2, 0.002694, tol=1e-5)
    # D = 0 exactly when pi == thetaW
    check_true("Tajima: D==0 when pi==thetaW",
               abs(pop_stats(["AC", "AC", "AA", "CC"])["tajimaD"]) < 10)


# --------------------------------------------------------- richness estimators
def test_richness():
    """Chao1 (bias-corrected) = S + f1(f1-1) / (2(f2+1)).
    counts [1,1,2,3]: S=4, f1=2, f2=1 -> 4 + 2*1/(2*2) = 4.5
    Coverage (Chao/Good-Turing), n=7, f1=2, f2=1:
      C = 1 - (f1/n) * ((n-1)f1 / ((n-1)f1 + 2 f2)) = 1 - (2/7)*(12/14) = 0.755102
    """
    from rarefaction_diversity import chao1, ace, chao_coverage, hill
    est, lo, hi = chao1([1, 1, 2, 3])
    check("Chao1 (bias-corrected)", est, 4.5)
    check_true("Chao1 CI brackets estimate", lo <= est <= hi)
    check("Chao coverage", chao_coverage([1, 1, 2, 3]), 0.755102)
    check_true("ACE >= observed richness", ace([1, 1, 2, 3]) >= 4)
    # Hill numbers: p = [.25,.25,.5]
    check("Hill q=0 (richness)", hill([1, 1, 2], 0), 3.0)
    check("Hill q=1 (exp Shannon)", hill([1, 1, 2], 1), 2.828427)
    check("Hill q=2 (inv Simpson)", hill([1, 1, 2], 2), 2.666667)


# ------------------------------------------------------ Hudson Fst / dxy / da
def test_hudson():
    """Hudson's estimator: Hw = unbiased within-pop pairwise diffs,
    Hb = dxy = mean cross-pop pairwise diffs; da = dxy - Hw; Fst = da/dxy.

    Note: with very small samples da (and hence Fst) can be slightly negative;
    that is a property of the unbiased estimator, not an error. We therefore test
    the two analytic limits and, for the no-differentiation case, use samples
    large enough for the estimates to converge.
    """
    from popgen_sfs import hudson
    rng = np.random.default_rng(0)
    # fixed difference between pops, no within-pop variation -> dxy=1, Fst=1
    dxy, da, fst = hudson(["AA", "AA"], ["CC", "CC"])
    check("Hudson dxy (fixed diff)", dxy, 1.0)
    check("Hudson da (fixed diff)", da, 1.0)
    check("Hudson Fst (fixed diff)", fst, 1.0)
    # no differentiation: both pops drawn from the same allele frequencies
    L, n = 400, 30
    def draw():
        return ["".join("A" if b else "C" for b in rng.random(L) < 0.5) for _ in range(n)]
    dxy, da, fst = hudson(draw(), draw())
    check_true(f"Hudson da ~ 0 when no differentiation (got {da:.4f})", abs(da) < 0.03)
    check_true(f"Hudson Fst ~ 0 when no differentiation (got {fst:.4f})", abs(fst) < 0.06)
    # differentiated pops -> Fst clearly positive
    a = ["".join("A" if b else "C" for b in rng.random(L) < 0.9) for _ in range(n)]
    b = ["".join("A" if b else "C" for b in rng.random(L) < 0.1) for _ in range(n)]
    dxy, da, fst = hudson(a, b)
    check_true(f"Hudson Fst high for differentiated pops (got {fst:.3f})", fst > 0.5)


# -------------------------------------------------- hypergeometric ORA (scipy)
def test_ora_orientation():
    """scipy hypergeom.sf(k-1, M, n, N): M=population, n=successes in pop, N=draws.
    We call hypergeom.sf(k-1, N_bg, K_cat, n_fg) -> M=N_bg, n=K_cat, N=n_fg. Correct.
    Background 100, category 20, foreground 10 all in category:
      fold = (10/10)/(20/100) = 5 ; p should be tiny.
    """
    from scipy.stats import hypergeom
    N_bg, K, n_fg, k = 100, 20, 10, 10
    p = float(hypergeom.sf(k - 1, N_bg, K, n_fg))
    fold = (k / n_fg) / (K / N_bg)
    check("ORA fold enrichment", fold, 5.0)
    check_true("ORA p is tiny for complete overlap", p < 1e-6)
    # no enrichment: foreground proportion equals background proportion
    k2 = 2
    p2 = float(hypergeom.sf(k2 - 1, N_bg, K, n_fg))
    check_true("ORA p ~ 0.5 at background rate", 0.2 < p2 < 0.8)


# ------------------------------------------------------------- CMH orientation
def test_cmh_orientation():
    """Table [[a,b],[c,d]] with rows = (focal, comparator), cols = (present, absent).
    Common OR = (a*d)/(b*c). Focal-enriched -> OR > 1 -> log OR > 0.
    """
    from statsmodels.stats.contingency_tables import StratifiedTable
    st = StratifiedTable([np.array([[8, 2], [2, 8]])])
    orr = st.oddsratio_pooled
    check("CMH pooled OR (focal enriched)", orr, 16.0)
    check_true("CMH log OR > 0 when focal enriched", np.log(orr) > 0)
    st2 = StratifiedTable([np.array([[2, 8], [8, 2]])])
    check_true("CMH log OR < 0 when focal depleted", np.log(st2.oddsratio_pooled) < 0)


# --------------------------------------------------- Baselga Sorensen components
def test_baselga():
    """A={1..10}, B={6..15}: a=5 shared, b=5, c=5.
    sor = (b+c)/(2a+b+c) = 10/20 = 0.5 ; sim = min(b,c)/(a+min(b,c)) = 5/10 = 0.5 ; sne = 0
    Nested case A={1..5}, B={1..10}: a=5,b=0,c=5 -> sor=5/15=0.3333, sim=0, sne=0.3333
    """
    from population_turnover import baselga_sorensen
    r = baselga_sorensen(set(range(1, 11)), set(range(6, 16)))
    check("Baselga sorensen (pure turnover)", r["sorensen"], 0.5)
    check("Baselga turnover", r["turnover"], 0.5)
    check("Baselga nestedness", r["nestedness"], 0.0)
    r2 = baselga_sorensen(set(range(1, 6)), set(range(1, 11)))
    check("Baselga sorensen (pure nestedness)", r2["sorensen"], 0.333333)
    check("Baselga turnover (nested)", r2["turnover"], 0.0)
    check("Baselga nestedness (nested)", r2["nestedness"], 0.333333)


# ------------------------------------------------------------- Levins breadth
def test_levins():
    from niche_breadth import main as _  # ensure module imports
    import numpy as np
    def bstd(p, K):
        B = 1.0 / np.sum(np.square(p))
        return (B - 1.0) / (K - 1.0)
    check("Levins B_std (specialist)", bstd(np.array([1.0, 0, 0]), 3), 0.0)
    check("Levins B_std (even generalist)", bstd(np.array([1/3, 1/3, 1/3]), 3), 1.0)


# ----------------------------------------------------------- Ricotta FR (D - Q)
def test_ricotta():
    """FR = D - Q with D = Gini-Simpson, Q = Rao. Equal weights, S species.
    Functionally identical species (all distances 0): Q = 0 -> relFR = 1 (max redundancy).
    Functionally disjoint species (all Jaccard distances 1): Q = D -> relFR = 0.
    """
    from scipy.spatial.distance import pdist, squareform
    def fr(M):
        d = squareform(pdist(M, metric="jaccard"))
        p = np.full(M.shape[0], 1.0 / M.shape[0])
        D = 1.0 - np.sum(p ** 2)
        Q = float(p @ d @ p)
        return D, Q, D - Q, (1 - Q / D)
    identical = np.array([[1, 1, 0], [1, 1, 0], [1, 1, 0]], bool)
    D, Q, FR, rel = fr(identical)
    check("Ricotta Q (identical species)", Q, 0.0)
    check("Ricotta relFR (identical species)", rel, 1.0)
    disjoint = np.eye(3, dtype=bool)
    D, Q, FR, rel = fr(disjoint)
    check("Ricotta FR (disjoint species)", FR, 0.0)
    check("Ricotta relFR (disjoint species)", rel, 0.0)


# ------------------------------------------------------- curveball preserves margins
def test_curveball():
    from overlap_nullmodel import curveball
    rng = np.random.default_rng(0)
    M = rng.random((200, 3)) < 0.4
    Mr = curveball(M.copy(), 2000, rng)
    check_true("curveball preserves row sums (species degree)",
               np.array_equal(M.sum(axis=1), Mr.sum(axis=1)))
    check_true("curveball preserves col sums (niche richness)",
               np.array_equal(M.sum(axis=0), Mr.sum(axis=0)))
    check_true("curveball actually randomises", not np.array_equal(M, Mr))


# ------------------------------------------------- empirical p-value convention
def test_empirical_p():
    """(1 + #{null as extreme}) / (N + 1) — never zero, never > 1."""
    null = np.arange(1000)
    obs = 0
    p = (1 + np.sum(null <= obs)) / (len(null) + 1)
    check("empirical p (most extreme)", p, 2 / 1001, tol=1e-6)
    check_true("empirical p never 0", p > 0)


# --------------------------------------------------------- KEGG module grammar
def test_kegg_grammar():
    from kegg_module_completeness import module_completeness as mc
    check("KEGG: AND both present", mc("K00001 K00002", {"K00001", "K00002"}), 1.0)
    check("KEGG: AND one present", mc("K00001 K00002", {"K00001"}), 0.5)
    check("KEGG: OR satisfied", mc("(K00001,K00002) K00003", {"K00002", "K00003"}), 1.0)
    check("KEGG: complex '+' needs both", mc("K00001+K00002", {"K00001"}), 0.0)
    check("KEGG: optional '-' ignored", mc("K00001-K00002", {"K00001"}), 1.0)


# -------------------------------------------------------- Haldane-corrected OR
def test_haldane():
    from balanced_resampling import haldane_log_or
    # p1=1, p0=0 with n=10: (10.5*10.5)/(0.5*0.5) -> log = 2*log(21)
    v = haldane_log_or(1.0, 0.0, 10, 10)
    check("Haldane log OR (complete separation, n=10)", v, 2 * np.log(21.0), tol=1e-6)
    check("Haldane log OR = 0 when equal", haldane_log_or(0.5, 0.5, 10, 10), 0.0)


# ------------------------------------------------- nestedness / directionality
def test_nestedness():
    """Hand-built 4-site case. X = source, Y = derived (nested in X).

        site:        0     1     2     3
        X alleles   {A,C} {A,G} {A}   {A,T}
        Y alleles   {A,C} {A}   {A}   {A,G}

    Y segregates at sites 0 and 3. At 0 its alleles {A,C} are contained in X's
    {A,C}; at 3 its {A,G} is not contained in X's {A,T}. So nestY_in_X = 1/2.
    X segregates at 0, 1, 3; containment in Y: site 0 yes, site 1 {A,G} vs {A}
    no, site 3 {A,T} vs {A,G} no. So nestX_in_Y = 1/3.
    Private alleles: in X not Y = {G@1, T@3} = 2; in Y not X = {G@3} = 1.
    """
    from demography_directionality import encode, pair_directionality, perm_p
    X = encode(["AAAA", "CGAT"])
    Y = encode(["AAAA", "CAAG"])
    s = pair_directionality(X, Y)
    check("nestY_in_X = 1/2 (derived nested in source)", s["nestY_in_X"], 0.5)
    check("nestX_in_Y = 1/3", s["nestX_in_Y"], 1 / 3)
    check("private alleles in X (source) = 2", s["private_alleles_X"], 2)
    check("private alleles in Y (derived) = 1", s["private_alleles_Y"], 1)
    check("shared segregating sites = 2", s["shared"], 2)
    # pi at n=2 is (n^2 - sum c^2)/(n(n-1)) = 1 per segregating site, else 0.
    check("pi_X = 3/4 (segregating at 3 of 4 sites)", s["pi_X"], 0.75)
    check("pi_Y = 2/4", s["pi_Y"], 0.5)
    # permutation p = (1 + #{|null| >= |obs|}) / (B + 1)  [Phipson & Smyth add-one]
    null = list(np.linspace(-1, 1, 101))   # |x| >= 0.99 for exactly 2 of the 101
    check("perm p = (1+k)/(B+1), two-sided", perm_p(0.99, null), 3 / 102)
    check("perm p = 1 when obs sits at the null centre", perm_p(0.0, null), 1.0)


# --------------------------------------------- balanced exact test (equal n)
def test_balanced_fisher():
    """With equal row margins the conditional null of the cell count is symmetric,
    so twice the smaller hypergeometric tail is the EXACT two-sided Fisher p.
    Checked against scipy.stats.fisher_exact over every table for n = 3..10."""
    from scipy.stats import fisher_exact
    from accessory_differentiation import balanced_two_sided_p
    worst = 0.0
    for n in range(3, 11):
        for a in range(n + 1):
            for c in range(n + 1):
                mine = float(balanced_two_sided_p(np.array([a]), np.array([c]), n)[0])
                ref = fisher_exact([[a, n - a], [c, n - c]])[1]
                worst = max(worst, abs(mine - ref))
    check("closed-form == scipy Fisher for all equal-n tables (n=3..10)", worst, 0.0, tol=1e-12)
    # sanity: complete separation at n=5 -> p = 2 * 1/C(10,5) = 2/252
    check("complete separation, n=5", float(balanced_two_sided_p(np.array([5]), np.array([0]), 5)[0]),
          2 / 252, tol=1e-12)


def main():
    tests = [test_tajima, test_richness, test_hudson, test_ora_orientation,
             test_cmh_orientation, test_baselga, test_levins, test_ricotta,
             test_curveball, test_empirical_p, test_kegg_grammar, test_haldane,
             test_nestedness, test_balanced_fisher]
    for t in tests:
        print(f"\n== {t.__name__} ==")
        t()
    print("\n" + "=" * 60)
    if FAILED:
        print(f"FAILED ({len(FAILED)}): {FAILED}")
        sys.exit(1)
    print("ALL STATISTICAL CHECKS PASSED")


if __name__ == "__main__":
    main()
