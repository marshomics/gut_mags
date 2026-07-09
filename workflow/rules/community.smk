# =============================================================================
# community.smk -- community & metabolic-interaction module (target
# `snakemake community_all`). Uses the per-genome tables in input/.
#
# Confounders handled consistently: community_ingest collapses strains to species
# (species-level prevalence); community_interactions rarefies every niche pool to
# equal species number. The microbial "which are human-specific" question runs the
# four new layers (auxotrophy/carbon/trait/modulep) through the differential +
# consensus stack automatically (they are added to DIFF_LAYERS in the Snakefile).
# =============================================================================

CM = f"{RESULTS}/12_community"
COMM = config["community"]
CTYPE_INPUT = {
    "auxotrophy": COMM["inputs"]["auxotrophy"],
    "carbon":     COMM["inputs"]["carbon"],
    "trait":      COMM["inputs"]["traits"],
    "modulep":    COMM["inputs"]["modules"],
}

rule community_all:
    input:
        f"{RESULTS}/report/community_report.html",

# ---- ingest each per-genome table to a species-level prevalence layer --------
rule community_ingest:
    input:
        table=lambda wc: CTYPE_INPUT[wc.ctype],
        samples=f"{RESULTS}/00_ingest/samples.parquet",
    output:
        f"{RESULTS}/03_profiles/prevalence_{{ctype}}.parquet",
    wildcard_constraints:
        ctype="auxotrophy|carbon|trait|modulep",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/community_ingest.py --config config/config.yaml --type {{wildcards.ctype}} "
        f"--input {{input.table}} --samples {{input.samples}} --out {{output}}"

# ---- community-level interaction metrics (rarefied) -------------------------
rule community_interactions:
    input:
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        prof=expand(f"{RESULTS}/03_profiles/prevalence_{{ctype}}.parquet",
                    ctype=["auxotrophy", "carbon", "trait", "modulep"]),
        mmap=COMM["cross_feeding"]["metabolite_map"],
    output:
        summ=f"{CM}/community_summary.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/community_interactions.py --config config/config.yaml "
        f"--species-table {{input.species}} --profiles-dir {RESULTS}/03_profiles "
        f"--metabolite-map {{input.mmap}} --out-dir {CM}"

rule fig_community:
    input: summ=rules.community_interactions.output.summ
    output: f"{RESULTS}/figures/community/community.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_community.py --config config/config.yaml "
        f"--community-dir {CM} --out {RESULTS}/figures/community/community"

# ---- report (pulls the microbial consensus on the community layers) ---------
rule community_report:
    input:
        summ=rules.community_interactions.output.summ,
        fig=rules.fig_community.output,
        micro=f"{RESULTS}/05_diff/combined/{config['synthesis']['focal_contrast']}_signatures_all.tsv",
    output:
        html=f"{RESULTS}/report/community_report.html",
        md=f"{RESULTS}/report/community_report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_community_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/community_report"
