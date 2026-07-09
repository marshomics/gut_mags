#!/usr/bin/env python3
"""
test_pipeline_smoke.py
----------------------
End-to-end smoke tests on small synthetic datasets with PLANTED signals. Each
chain must not only run, it must recover the signal that was planted, which is
what makes this a test rather than a syntax check.

    python tests/test_pipeline_smoke.py

Covers: ingest -> species units; the differential consensus; taxonomy specificity /
enrichment / overlap null; population genetics (founder bottleneck); functional
redundancy; Western vs non-Western turnover; the community metabolic module; and
ORA enrichment.
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(__file__))
from _shim import run  # noqa: E402  (installs the parquet shim)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE_CFG = os.path.join(ROOT, "config", "config.yaml")
FAILED = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{(' — ' + detail) if detail else ''}")
    if not cond:
        FAILED.append(name)


def write_cfg(tmp, **overrides):
    cfg = yaml.safe_load(open(BASE_CFG))
    for k, v in overrides.items():
        cur = cfg
        parts = k.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    p = os.path.join(tmp, "cfg.yaml")
    yaml.safe_dump(cfg, open(p, "w"))
    return p


# --------------------------------------------------------------------------
def make_metadata(tmp, western=False):
    rows = []
    def add(g, sp, niche, fam, gen, pop=""):
        r = dict(genome=g, gtdb_domain="d__Bacteria", gtdb_phylum="p__P", gtdb_class="c__C",
                 gtdb_order="o__O", gtdb_family=fam, gtdb_genus=gen, gtdb_species=sp,
                 source=niche, Completeness="97", Contamination="1", gc_avg="0.5",
                 cds_number="2000", genome_size="2000000", corrected_genome_size="2000000",
                 common_name_x=("Mouse" if niche == "animal" else ""),
                 n_contigs="30", ctg_L50="50000")
        for c in ["class", "domain", "family", "kingdom", "order", "phylum"]:
            r[c] = ""
        if western:
            r["western_nonwestern"] = pop
        rows.append(r)
    gi = 0
    plan = [("human", 6, "f__Fh", "g__Gh"), ("animal", 6, "f__Fa", "g__Ga"),
            ("free", 6, "f__Ff", "g__Gf")]
    for niche, nsp, fam, gen in plan:
        for s in range(nsp):
            sp = f"s__{niche}_sp{s}"
            for st in range(3):
                gi += 1
                pop = ("western" if s < 3 else "non_western") if niche == "human" else ""
                add(f"G{gi:04d}", sp, niche, fam, gen, pop)
    df = pd.DataFrame(rows)
    df["qc"] = "92.0"
    p = os.path.join(tmp, "meta.tsv")
    df.to_csv(p, sep="\t", index=False)
    return p


def prevalence(tmp, name, mapping):
    """mapping: feature -> list of species carrying it."""
    rows = []
    for feat, spp in mapping.items():
        for s in spp:
            rows.append(dict(species=s, feature=feat, prevalence=1.0, present=1,
                             mean_copies=1.0, n_genomes=3))
    p = os.path.join(tmp, f"prevalence_{name}.parquet")
    pd.DataFrame(rows).to_parquet(p)
    return p


# --------------------------------------------------------------------------
def test_ingest_and_units(tmp):
    print("\n== ingest + species units (+ population labels) ==")
    meta = make_metadata(tmp, western=True)
    cfg = write_cfg(tmp, **{"inputs.metadata": meta, "population.enabled": True,
                            "population.min_species_per_population": 2})
    run("ingest_metadata", ["--config", cfg, "--out-prefix", f"{tmp}/samples"])
    run("define_species_units", ["--config", cfg, "--samples", f"{tmp}/samples.parquet",
                                 "--out-prefix", f"{tmp}/su"])
    st = pd.read_csv(f"{tmp}/su_species_table.tsv", sep="\t")
    qc = json.load(open(f"{tmp}/samples_qc_report.json"))
    check("all 54 genomes pass QC", qc["n_passed"] == 54, f"{qc['n_passed']}")
    check("18 species defined", len(st) == 18, f"{len(st)}")
    check("6 human specialists", (st["specialist_niche"] == "human").sum() == 6)
    hum = st[st["niche_primary"] == "human"]
    check("population split 3 western / 3 non_western",
          (hum["population"] == "western").sum() == 3 and (hum["population"] == "non_western").sum() == 3)
    return cfg


def test_differential_consensus(tmp, cfg):
    print("\n== differential prep + CMH + resampling + consensus (planted signal) ==")
    st = pd.read_csv(f"{tmp}/su_species_table.tsv", sep="\t")
    human = st.loc[st["specialist_niche"] == "human", "species"].tolist()
    others = st.loc[st["specialist_niche"].isin(["animal", "free"]), "species"].tolist()
    # SIG present only in human; BG present everywhere
    prevalence(tmp, "ko", {"SIG": human, "BG": human + others})
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "population.enabled": True,
                             "stats.cmh.min_stratum_size": 1,
                             "stats.min_abs_log2_or": 0.5,
                             "functional_annotation.min_prevalence_any_niche": 0.01,
                             "stats.resampling.iterations": 100})
    for contrast, g1, g0 in [("human_vs_rest", 6, 12), ("free_vs_rest", 6, 12),
                             ("host_vs_free", 12, 6), ("western_vs_nonwestern", 3, 3)]:
        run("differential_prep", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                                  "--species-table", f"{tmp}/su_species_table.tsv",
                                  "--prevalence", f"{tmp}/prevalence_ko.parquet",
                                  "--layer", "ko", "--contrast", contrast,
                                  "--out-prefix", f"{tmp}/d_{contrast}"])
        a = pd.read_csv(f"{tmp}/d_{contrast}_analysis_species.tsv", sep="\t")
        check(f"contrast {contrast} groups", (a.group == 1).sum() == g1 and (a.group == 0).sum() == g0,
              f"{(a.group==1).sum()} vs {(a.group==0).sum()}")

    run("cmh_stratified", ["--config", cfg2, "--analysis", f"{tmp}/d_human_vs_rest_analysis_species.tsv",
                           "--presence", f"{tmp}/d_human_vs_rest_presence.parquet",
                           "--out", f"{tmp}/cmh.tsv"])
    run("balanced_resampling", ["--config", cfg2, "--analysis", f"{tmp}/d_human_vs_rest_analysis_species.tsv",
                                "--presence", f"{tmp}/d_human_vs_rest_presence.parquet",
                                "--species-table", f"{tmp}/su_species_table.tsv",
                                "--contrast", "human_vs_rest", "--out", f"{tmp}/rs.tsv"])
    feats = [l.strip() for l in open(f"{tmp}/d_human_vs_rest_tested_features.txt") if l.strip()]
    pd.DataFrame({"feature": feats,
                  "estimate_log_or": [4.0 if f == "SIG" else 0.0 for f in feats],
                  "se": 0.5, "z": 5, "p": [1e-6 if f == "SIG" else 0.9 for f in feats],
                  "n_species": 18, "n_present": 6, "status": "ok"}).to_csv(f"{tmp}/pg.tsv", sep="\t", index=False)
    pd.DataFrame({"feature": feats,
                  "scoary_log2or": [3.0 if f == "SIG" else 0.0 for f in feats],
                  "scoary_fisher_q": [1e-4 if f == "SIG" else 0.9 for f in feats],
                  "scoary_empirical_p": [0.001 if f == "SIG" else 0.8 for f in feats],
                  "scoary_supporting": [5 if f == "SIG" else 1 for f in feats],
                  "scoary_best_p": [0.01 if f == "SIG" else 0.5 for f in feats],
                  "scoary_sig": [f == "SIG" for f in feats],
                  "scoary_dir": [1.0 if f == "SIG" else np.nan for f in feats]}).to_csv(
        f"{tmp}/scoary.tsv", sep="\t", index=False)
    run("consensus_signatures", ["--config", cfg2, "--phyloglm", f"{tmp}/pg.tsv",
                                 "--cmh", f"{tmp}/cmh.tsv", "--resampling", f"{tmp}/rs.tsv",
                                 "--scoary", f"{tmp}/scoary.tsv", "--layer", "ko",
                                 "--contrast", "human_vs_rest", "--out-prefix", f"{tmp}/sig"])
    d = pd.read_csv(f"{tmp}/sig.tsv", sep="\t")
    js = json.load(open(f"{tmp}/sig.json"))
    sig = d[d.feature == "SIG"].iloc[0]
    bg = d[d.feature == "BG"].iloc[0]
    check("4 methods used", js["methods_used"] == ["phyloglm", "cmh", "resampling", "scoary"])
    # every species is niche-specific by construction, so no genus contains both
    # groups: CMH is UNTESTABLE, not negative, and must not veto the call.
    check("CMH reported untestable (niche fully nested in clade)",
          sig["methods_untestable"] == "cmh" and js["n_untestable_per_method"]["cmh"] == 2)
    check("planted SIG is consensus, human_enriched",
          bool(sig["consensus_signature"]) and sig["direction"] == "human_enriched")
    check("SIG is tier2 (one required method untestable)",
          sig["tier"] == "tier2_consensus_partial", str(sig["tier"]))
    check("background BG is not a signature", not bool(bg["consensus_signature"]))
    check("BG direction is 'ns', not a spurious depletion", bg["direction"] == "ns")
    check("consensus effect size finite", np.isfinite(sig["consensus_log2or"]))

    # An APPLICABLE but non-significant method must veto; an untestable ANCHOR
    # must block the call outright.
    cmh = pd.read_csv(f"{tmp}/cmh.tsv", sep="\t")
    cmh.loc[cmh.feature == "SIG", ["mh_log_or", "p", "n_strata", "status"]] = [0.0, 0.9, 3, "ok"]
    cmh.to_csv(f"{tmp}/cmh_veto.tsv", sep="\t", index=False)
    run("consensus_signatures", ["--config", cfg2, "--phyloglm", f"{tmp}/pg.tsv",
                                 "--cmh", f"{tmp}/cmh_veto.tsv", "--resampling", f"{tmp}/rs.tsv",
                                 "--scoary", f"{tmp}/scoary.tsv", "--layer", "ko",
                                 "--contrast", "human_vs_rest", "--out-prefix", f"{tmp}/sig_veto"])
    dv = pd.read_csv(f"{tmp}/sig_veto.tsv", sep="\t")
    check("an applicable-but-null CMH vetoes the call",
          not bool(dv[dv.feature == "SIG"].iloc[0]["consensus_signature"]))
    pgf = pd.read_csv(f"{tmp}/pg.tsv", sep="\t")
    pgf["status"] = "fit_failed"; pgf["p"] = np.nan
    pgf.to_csv(f"{tmp}/pg_fail.tsv", sep="\t", index=False)
    run("consensus_signatures", ["--config", cfg2, "--phyloglm", f"{tmp}/pg_fail.tsv",
                                 "--cmh", f"{tmp}/cmh.tsv", "--resampling", f"{tmp}/rs.tsv",
                                 "--scoary", f"{tmp}/scoary.tsv", "--layer", "ko",
                                 "--contrast", "human_vs_rest", "--out-prefix", f"{tmp}/sig_anchor"])
    da = pd.read_csv(f"{tmp}/sig_anchor.tsv", sep="\t")
    check("an untestable anchor method blocks every call",
          not da["consensus_signature"].any())
    return cfg2


def test_taxonomy(tmp, cfg):
    print("\n== taxonomy: specificity null, enrichment, overlap null, novelty ==")
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "taxonomy.null_iterations": 99,
                             "taxonomy.overlap_null.iterations": 50,
                             "taxonomy.enrichment.min_taxon_species": 3})
    os.makedirs(f"{tmp}/spec", exist_ok=True)
    run("taxon_niche_specificity", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                                    "--species-table", f"{tmp}/su_species_table.tsv",
                                    "--out-dir", f"{tmp}/spec"])
    g = pd.read_csv(f"{tmp}/spec/taxon_specificity_genus.tsv", sep="\t")
    check("each genus is niche-specific (Levins B_std == 0)", (g["levins_B_std"] == 0).all())
    check("specificity p is significant for a 6-species specialist genus",
          (g["q_specificity"] < 0.05).any(), str(g["q_specificity"].tolist()))
    run("taxon_enrichment", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                             "--species-table", f"{tmp}/su_species_table.tsv",
                             "--out-dir", f"{tmp}/enr"])
    e = pd.read_csv(f"{tmp}/enr/enrichment_genus.tsv", sep="\t")
    gh = e[(e.taxon == "g__Gh") & (e.niche == "human")].iloc[0]
    check("human genus enriched in human (log2 OR > 0, q<0.05)", gh.log2_or > 0 and gh.q < 0.05)
    run("overlap_nullmodel", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                              "--out-prefix", f"{tmp}/ov"])
    ov = pd.read_csv(f"{tmp}/ov_nullmodel.tsv", sep="\t")
    check("overlap null reports the model used", "null_model" in ov.columns and
          ov["null_model"].iloc[0].startswith("curveball"))
    check("no cross-niche overlap observed (by construction)", (ov["observed"] == 0).all())
    run("taxonomic_novelty", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                              "--species-table", f"{tmp}/su_species_table.tsv",
                              "--out-prefix", f"{tmp}/nov", "--n-boot", "20"])
    check("novelty output written", os.path.exists(f"{tmp}/nov_by_niche.tsv"))


def test_popgen(tmp, cfg):
    print("\n== population genetics: founder bottleneck signature ==")
    # Source population (human): 20 haplotypes segregating at 300 of 600 sites.
    # Derived population (animal): founded by 3 of those haplotypes, expanded to
    # 16 with a little drift and 3 novel mutations. Its variation is therefore a
    # near-subset of the source's, with lower pi and fewer private alleles.
    rng = np.random.default_rng(7)
    L, nH, nA, nfound = 600, 20, 16, 3
    hum = np.full((nH, L), "A", dtype="<U1")
    seg = rng.choice(L, 300, replace=False)
    for s in seg:
        f = rng.uniform(0.1, 0.9)
        hum[rng.random(nH) < f, s] = "C"
    founders = hum[rng.choice(nH, nfound, replace=False)]
    ani = founders[rng.integers(0, nfound, nA)].copy()
    for s in rng.choice(L, 3, replace=False):          # novel mutations post-founding
        ani[rng.integers(nA), s] = "T"
    names = [f"H{i}" for i in range(nH)] + [f"A{i}" for i in range(nA)]
    mat = np.vstack([hum, ani])
    mat = mat[:, [j for j in range(L) if len(set(mat[:, j])) > 1]]
    with open(f"{tmp}/snps.fa", "w") as fh:
        for n, row in zip(names, mat):
            fh.write(f">{n}\n{''.join(row)}\n")
    pd.DataFrame([{"genome": n, "niche": "human" if n.startswith("H") else "animal",
                   "role": "focal"} for n in names]).to_csv(f"{tmp}/derep.tsv", sep="\t", index=False)
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "transition.popgen.bootstrap": 100})
    run("popgen_sfs", ["--config", cfg2, "--alignment", f"{tmp}/snps.fa",
                       "--niche-map", f"{tmp}/derep.tsv", "--out-dir", f"{tmp}/pg"])
    d = pd.read_csv(f"{tmp}/pg/popgen_diversity.tsv", sep="\t").set_index("niche")
    check("derived pop has lower pi (founder effect)",
          d.loc["animal", "pi_mean"] < d.loc["human", "pi_mean"],
          f"{d.loc['animal','pi_mean']:.3f} < {d.loc['human','pi_mean']:.3f}")
    check("derived pop has fewer segregating sites at equal n",
          d.loc["animal", "S_mean"] < d.loc["human", "S_mean"],
          f"{d.loc['animal','S_mean']:.0f} < {d.loc['human','S_mean']:.0f}")
    run("demography_directionality", ["--config", cfg2, "--alignment", f"{tmp}/snps.fa",
                                      "--niche-map", f"{tmp}/derep.tsv", "--out-dir", f"{tmp}/dm"])
    dd = pd.read_csv(f"{tmp}/dm/directionality.tsv", sep="\t").iloc[0]
    check("nestedness calls animal the derived population", dd["derived_call"] == "animal",
          f"{dd['derived_call']} ({dd['call_status']})")
    check("nestedness resolved with CI excluding 0", dd["call_status"] == "resolved_nestedness",
          f"nest_diff={dd['nest_diff']} CI=[{dd['nest_diff_ci_lo']},{dd['nest_diff_ci_hi']}]")
    check("all three directionality lines agree", dd["n_lines_agree"] == 3,
          f"{dd['n_lines_agree']}/{dd['n_lines_evaluated']}")
    check("derived pop has fewer private alleles",
          dd["private_alleles_X"] < dd["private_alleles_Y"],
          f"animal={dd['private_alleles_X']} human={dd['private_alleles_Y']}")

    # Accessory gene gain/loss with a deliberately imbalanced design: 20 human
    # strains vs 16 animal. One gene (GAIN) is fixed in animal, absent in human.
    # 40 neutral genes segregate identically in both. Under the unbalanced Fisher
    # the larger population wins more calls; the balanced draws must not.
    rng2 = np.random.default_rng(11)
    genes, rowsA = [], []
    for g in range(40):
        f = rng2.uniform(0.3, 0.7)
        rowsA.append((rng2.random(nH + nA) < f).astype(int)); genes.append(f"neutral{g}")
    gain = np.r_[np.zeros(nH, int), np.ones(nA, int)]
    rowsA.append(gain); genes.append("GAIN")
    rtab = pd.DataFrame(np.vstack(rowsA), index=genes, columns=names)
    rtab.index.name = "Gene"
    rtab.to_csv(f"{tmp}/gpa.Rtab", sep="\t")
    cfg3 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "transition.accessory.balanced_iterations": 50})
    run("accessory_differentiation", ["--config", cfg3, "--rtab", f"{tmp}/gpa.Rtab",
                                      "--niche-map", f"{tmp}/derep.tsv",
                                      "--out", f"{tmp}/acc.tsv", "--summary", f"{tmp}/acc_sum.tsv"])
    ac = pd.read_csv(f"{tmp}/acc.tsv", sep="\t").set_index("gene")
    asum = pd.read_csv(f"{tmp}/acc_sum.tsv", sep="\t").iloc[0]
    check("planted gene gain called in the derived niche",
          ac.loc["GAIN", "balanced_call"] == "animal", str(ac.loc["GAIN", "balanced_call"]))
    check("no neutral gene is called", (ac.drop("GAIN")["balanced_call"] == "ns").all())
    check("balanced counts: exactly one gene gained in animal",
          asum["enriched_in_animal"] == 1 and asum["enriched_in_human"] == 0)
    check("balanced subsample equalises strain counts", asum["n_balanced_subsample"] == nA)

    # Two populations drawn from the same pool must NOT yield a directional call.
    perm = rng.permutation(len(names))
    pd.DataFrame([{"genome": names[i], "role": "focal",
                   "niche": "human" if k < len(names) // 2 else "animal"}
                  for k, i in enumerate(perm)]).to_csv(f"{tmp}/derep_null.tsv", sep="\t", index=False)
    run("demography_directionality", ["--config", cfg2, "--alignment", f"{tmp}/snps.fa",
                                      "--niche-map", f"{tmp}/derep_null.tsv", "--out-dir", f"{tmp}/dmnull"])
    dn = pd.read_csv(f"{tmp}/dmnull/directionality.tsv", sep="\t").iloc[0]
    check("shuffled labels give no significant nestedness asymmetry",
          dn["nestedness_call"] == "unresolved" or dn["nest_p"] > 0.05,
          f"p={dn['nest_p']}")


def test_redundancy_and_population(tmp, cfg):
    print("\n== functional redundancy + Western/non-Western turnover ==")
    st = pd.read_csv(f"{tmp}/su_species_table.tsv", sep="\t")
    hum = st.loc[st["niche_primary"] == "human", "species"].tolist()
    west = st.loc[st["population"] == "western", "species"].tolist()
    nonw = st.loc[st["population"] == "non_western", "species"].tolist()
    prevalence(tmp, "ko", {f"CORE{i}": hum for i in range(5)} | {"RARE": hum[:1]})
    prevalence(tmp, "cazyme", {"GHW": west, "GHN": nonw, "GHS": hum})
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "population.enabled": True, "population.min_species_per_population": 2,
                             "redundancy.layers": ["ko", "cazyme"],
                             "redundancy.ricotta.max_species": 20,
                             "redundancy.accumulation.bootstrap": 10})
    run("functional_redundancy", ["--config", cfg2, "--samples", f"{tmp}/samples.parquet",
                                  "--species-table", f"{tmp}/su_species_table.tsv",
                                  "--profiles-dir", tmp, "--out-dir", f"{tmp}/red"])
    r = pd.read_csv(f"{tmp}/red/redundancy_summary.tsv", sep="\t")
    ko = r[(r.layer == "ko") & (r.community == "human")].iloc[0]
    check("KO layer is highly redundant (relFR high)", ko["ricotta_relFR"] > 0.8, f"{ko['ricotta_relFR']}")
    check("redundancy reports pool and annotated species counts",
          {"n_species_pool", "n_species"} <= set(r.columns))
    run("population_turnover", ["--config", cfg2, "--species-table", f"{tmp}/su_species_table.tsv",
                                "--profiles-dir", tmp, "--out-dir", f"{tmp}/pop"])
    t = pd.read_csv(f"{tmp}/pop/population_turnover.tsv", sep="\t").set_index("layer")
    check("taxonomic turnover complete between populations", t.loc["ko", "taxonomic_sorensen"] == 1.0)
    check("KO functions conserved despite species turnover",
          t.loc["ko", "functional_sorensen"] < 0.2, f"{t.loc['ko','functional_sorensen']}")
    check("CAZymes diverge between populations",
          t.loc["cazyme", "functional_sorensen"] > t.loc["ko", "functional_sorensen"])
    check("carrier substitution detected (shared functions, different species)",
          t.loc["ko", "median_carrier_jaccard_shared_fn"] == 0.0)


def test_community(tmp, cfg):
    print("\n== community: auxotrophy, cross-feeding, complementarity, traits ==")
    samples = pd.read_parquet(f"{tmp}/samples.parquet")
    AAs = yaml.safe_load(open(BASE_CFG))["community"]["amino_acids"]
    rng = np.random.default_rng(2)
    aa, cb, md, tr = [], [], [], []
    for _, r in samples.iterrows():
        g, niche = r["genome"], r["niche"]
        can = 0.2 if niche == "human" else (0.5 if niche == "animal" else 0.9)
        aa.append({"genome": g, **{a: int(rng.random() < can) for a in AAs}, "chorismate": 1})
        row = {"genome": g, "glucose": 1, "acetate": 0, "L-lactate": 0, "succinate": 0,
               "propionate": 0, "D-lactate": 0, "ethanol": 0, "pyruvate": 0, "fumarate": 0}
        if niche == "human":
            row.update({"acetate": 1, "L-lactate": 1})
        cb.append(row)
        m = {"file": g, "M00579": int(niche == "human"), "M00001": 1, "M00002": 0, "M00003": 0}
        if niche == "human":
            m["M00002" if g[-1] in "0123456789" and int(g[-1]) % 2 else "M00003"] = 1
        md.append(m)
        t = {"": g, "Anaerobe": int(niche == "human"), "Aerobe": int(niche == "free"),
             "Facultative": int(niche == "animal"), "Spore formation": int(niche == "human"),
             "Glucose fermenter": int(niche != "free"), "Bile-susceptible": 0, "Motile": int(niche == "free"),
             "Gram positive": 0, "Gram negative": 0}
        tr.append(t)
    for name, rows in [("aa", aa), ("carbon", cb), ("modules", md), ("traits", tr)]:
        pd.DataFrame(rows).to_csv(f"{tmp}/{name}.tsv", sep="\t", index=False)
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "community.rarefaction.bootstrap": 20})
    prof = f"{tmp}/cprof"; os.makedirs(prof, exist_ok=True)
    for ct, f in [("auxotrophy", "aa"), ("carbon", "carbon"), ("modulep", "modules"), ("trait", "traits")]:
        run("community_ingest", ["--config", cfg2, "--type", ct, "--input", f"{tmp}/{f}.tsv",
                                 "--samples", f"{tmp}/samples.parquet",
                                 "--out", f"{prof}/prevalence_{ct}.parquet"])
    # species_table must have >=10 species per niche pool for the community step
    st = pd.read_csv(f"{tmp}/su_species_table.tsv", sep="\t")
    check("community ingest produced 4 layers",
          all(os.path.exists(f"{prof}/prevalence_{c}.parquet")
              for c in ["auxotrophy", "carbon", "modulep", "trait"]))
    aux = pd.read_parquet(f"{prof}/prevalence_auxotrophy.parquet")
    hum = set(st.loc[st.niche_primary == "human", "species"])
    free = set(st.loc[st.niche_primary == "free", "species"])
    hrate = aux[aux.species.isin(hum)]["present"].mean()
    frate = aux[aux.species.isin(free)]["present"].mean()
    check("human species more auxotrophic than free-living", hrate > frate, f"{hrate:.2f} vs {frate:.2f}")


def test_ora(tmp, cfg):
    print("\n== ORA enrichment (proper background, per direction) ==")
    # 20 tested features: 10 GH, 10 GT. All 8 human-enriched signatures are GH.
    # Hypergeometric: k=8 of K=10 GH drawn in n=8 from N=20 -> p = 6.2e-05.
    ghs = [f"GH{i}" for i in range(10)]
    gts = [f"GT{i}" for i in range(10)]
    rows = []
    for i, f in enumerate(ghs):
        hit = i < 8
        rows.append(dict(feature=f, layer="cazyme", consensus_signature=hit,
                         direction="human_enriched" if hit else "ns",
                         consensus_log2or=3.0 - 0.1 * i if hit else 0.05,
                         pg_q=1e-5 if hit else 0.5))
    for i, f in enumerate(gts):
        rows.append(dict(feature=f, layer="cazyme", consensus_signature=False,
                         direction="ns", consensus_log2or=-0.05, pg_q=0.6))
    pd.DataFrame(rows).to_csv(f"{tmp}/sigs.tsv", sep="\t", index=False)
    gs = pd.DataFrame([{"system": "cazy_class", "category_id": "GH",
                        "category_name": "Glycoside hydrolases", "feature": f} for f in ghs] +
                      [{"system": "cazy_class", "category_id": "GT",
                        "category_name": "Glycosyltransferases", "feature": f} for f in gts])
    gs.to_csv(f"{tmp}/gs.tsv", sep="\t", index=False)
    cfg2 = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                             "enrichment.min_set_size": 1})
    run("ora_enrichment", ["--config", cfg2, "--signatures", f"{tmp}/sigs.tsv",
                           "--genesets", f"{tmp}/gs.tsv", "--layer", "cazyme",
                           "--contrast", "human_vs_rest", "--out", f"{tmp}/ora.tsv"])
    o = pd.read_csv(f"{tmp}/ora.tsv", sep="\t")
    up = o[(o.direction == "up") & (o.category_id == "GH")].iloc[0]
    check("ORA background = tested annotatable features (N=20)", up["N"] == 20, f"N={up['N']}")
    check("GH over-represented among human signatures (fold=2)",
          abs(up["fold_enrichment"] - 2.0) < 1e-6, f"{up['fold_enrichment']}")
    # P(X>=8) = C(10,8)C(10,0)/C(20,8) = 45/125970
    check("GH p matches the hypergeometric by hand (45/125970)",
          abs(up["p"] - 45 / 125970) < 1e-12, f"p={up['p']:.4g}")
    check("GH significant after BH", up["q"] < 0.05, f"q={up['q']:.4g}")
    gt_up = o[(o.direction == "up") & (o.category_id == "GT")]
    check("GT not enriched among human signatures",
          len(gt_up) == 0 or gt_up.iloc[0]["q"] > 0.05)


def test_annotation_paths(tmp):
    print("\n== annotation paths: manifest mode, template mode, validation ==")
    import subprocess
    d = os.path.join(tmp, "ap"); os.makedirs(d + "/raw", exist_ok=True)
    G = [f"P{i:03d}" for i in range(4)]
    pd.DataFrame({"genome": G}).to_parquet(f"{d}/samples.parquet")

    # files on disk with irregular basenames, as a real manifest would have
    rows = []
    for i, g in enumerate(G):
        faa = f"{d}/raw/study{i}_orig.faa"; open(faa, "w").write(">x\nMK\n")
        kof = f"{d}/raw/weird-{i}.kofam.tsv"; open(kof, "w").write("")
        egg = f"{d}/raw/e{i}.annotations"; open(egg, "w").write("")
        rows.append({"genome": g, "faa": faa, "kofam": kof, "emapper": egg,
                     "unknown_column": "ignored"})
    pd.DataFrame(rows).to_csv(f"{d}/manifest.tsv", sep="\t", index=False)

    cfg_m = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                              "inputs.annotation_manifest": f"{d}/manifest.tsv",
                              "inputs.validate_paths": "all"})
    run("resolve_annotation_paths", ["--config", cfg_m,
                                     "--samples", f"{d}/samples.parquet",
                                     "--out", f"{d}/ap_manifest.tsv"])
    ap = pd.read_csv(f"{d}/ap_manifest.tsv", sep="\t", keep_default_na=False).set_index("genome")
    check("manifest: aliased columns resolved (faa/emapper/kofam)",
          ap.loc["P000", "prokka_faa"].endswith("study0_orig.faa") and
          ap.loc["P000", "eggnog"].endswith("e0.annotations") and
          ap.loc["P000", "kofam"].endswith("weird-0.kofam.tsv"))
    check("manifest: pipeline-generated kinds still come from templates",
          ap.loc["P001", "amrfinder_tsv"].endswith("amrfinder/P001.tsv"))
    check("manifest: unsupplied kind left empty, not guessed",
          ap.loc["P000", "assembly"] == "")

    # template mode reproduces the same schema
    cfg_t = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                              "inputs.annotation_manifest": None,
                              "inputs.validate_paths": "none",
                              "inputs.annotations.prokka_faa": f"{d}/raw/{{genome}}.faa",
                              "inputs.annotations.eggnog_annotations": f"{d}/raw/{{genome}}.egg",
                              "inputs.annotations.kofam_annotations": f"{d}/raw/{{genome}}.ko"})
    run("resolve_annotation_paths", ["--config", cfg_t,
                                     "--samples", f"{d}/samples.parquet",
                                     "--out", f"{d}/ap_tmpl.tsv"])
    at = pd.read_csv(f"{d}/ap_tmpl.tsv", sep="\t", keep_default_na=False).set_index("genome")
    check("template: {genome} expanded", at.loc["P002", "prokka_faa"].endswith("/P002.faa"))
    check("both modes give the same columns", list(ap.columns) == list(at.columns))

    # a genome absent from the manifest must be an error, never a silent skip
    pd.DataFrame(rows[:2]).to_csv(f"{d}/manifest_short.tsv", sep="\t", index=False)
    cfg_s = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                              "inputs.annotation_manifest": f"{d}/manifest_short.tsv",
                              "inputs.validate_paths": "none"})
    try:
        run("resolve_annotation_paths", ["--config", cfg_s,
                                         "--samples", f"{d}/samples.parquet",
                                         "--out", f"{d}/ap_short.tsv"])
        ok = False
    except SystemExit:
        ok = True
    check("a genome missing from the manifest is a hard error", ok)

    # a declared path that does not exist is caught before cluster time
    bad = [dict(r) for r in rows]
    bad[0]["faa"] = f"{d}/raw/does_not_exist.faa"
    pd.DataFrame(bad).to_csv(f"{d}/manifest_bad.tsv", sep="\t", index=False)
    cfg_b = write_cfg(tmp, **{"inputs.metadata": os.path.join(tmp, "meta.tsv"),
                              "inputs.annotation_manifest": f"{d}/manifest_bad.tsv",
                              "inputs.validate_paths": "all"})
    try:
        run("resolve_annotation_paths", ["--config", cfg_b,
                                         "--samples", f"{d}/samples.parquet",
                                         "--out", f"{d}/ap_bad.tsv"])
        ok = False
    except SystemExit:
        ok = True
    check("a nonexistent declared path is caught up front", ok)

    # the symlink farm renames irregular basenames to <genome>.<ext>, which is
    # what keeps Mash/Panaroo sample names equal to genome ids
    with open(f"{d}/faa_paths.tsv", "w") as fh:
        for r in rows:
            fh.write(f"{r['genome']}\t{r['faa']}\n")
    farm_sh = os.path.join(ROOT, "scripts", "sh", "symlink_farm.sh")
    subprocess.run(["bash", farm_sh, f"{d}/faa_paths.tsv", f"{d}/farm", "faa"], check=True,
                   capture_output=True)
    linked = sorted(f for f in os.listdir(f"{d}/farm") if f.endswith(".faa"))
    check("symlink farm names files <genome>.faa", linked == [f"{g}.faa" for g in G],
          str(linked))
    check("symlink targets resolve to the real files",
          os.path.realpath(f"{d}/farm/P000.faa") == os.path.realpath(rows[0]["faa"]))
    open(f"{d}/faa_missing.tsv", "w").write(f"P000\t{d}/raw/nope.faa\n")
    r = subprocess.run(["bash", farm_sh, f"{d}/faa_missing.tsv", f"{d}/farm3", "faa"],
                       capture_output=True)
    check("symlink farm fails on a missing source rather than dangling", r.returncode != 0)


def main():
    tmp = tempfile.mkdtemp(prefix="hgn_smoke_")
    print(f"workdir: {tmp}")
    try:
        cfg = test_ingest_and_units(tmp)
        test_annotation_paths(tmp)
        cfg2 = test_differential_consensus(tmp, cfg)
        test_taxonomy(tmp, cfg)
        test_popgen(tmp, cfg)
        test_redundancy_and_population(tmp, cfg)
        test_community(tmp, cfg)
        test_ora(tmp, cfg)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("\n" + "=" * 60)
    if FAILED:
        print(f"FAILED ({len(FAILED)}): {FAILED}")
        sys.exit(1)
    print("ALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()
