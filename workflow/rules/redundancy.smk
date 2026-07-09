# =============================================================================
# redundancy.smk -- functional redundancy of the human gut, and the Western vs
# non-Western comparison. Targets:
#   snakemake redundancy_all   redundancy metrics + figure + report
#   snakemake western_all      adds the Western/non-Western turnover, the
#                              western_vs_nonwestern 4-method consensus + enrichment
#                              (these flow through automatically when
#                              population.enabled), and the population figure.
# The population analyses are inert until the western_nonwestern column is added
# and population.enabled is set true.
# =============================================================================

RED = f"{RESULTS}/10_redundancy"
POP = f"{RESULTS}/11_population"
POP_ENABLED = config.get("population", {}).get("enabled", False)
POP_CONTRAST = config.get("population", {}).get("contrast", "western_vs_nonwestern")

rule redundancy_all:
    input:
        f"{RESULTS}/report/redundancy_report.html",

rule western_all:
    input:
        f"{RESULTS}/report/redundancy_report.html",
        # when enabled, also build the western_vs_nonwestern functional divergence
        (expand(f"{RESULTS}/05_diff/{{layer}}/{POP_CONTRAST}_signatures.tsv", layer=DIFF_LAYERS)
         if POP_ENABLED else []),
        (f"{RESULTS}/figures/enrichment/dotplot_{POP_CONTRAST}.png" if POP_ENABLED else []),

# ---- functional redundancy (human gut) --------------------------------------
rule functional_redundancy:
    input:
        samples=f"{RESULTS}/00_ingest/samples.parquet",
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        # anchor on the profiles so this waits for the species profiles to exist
        anchor=f"{RESULTS}/03_profiles/prevalence_ko.parquet",
    output:
        summ=f"{RED}/redundancy_summary.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/functional_redundancy.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} "
        f"--profiles-dir {RESULTS}/03_profiles --out-dir {RED}"

rule fig_redundancy:
    input: summ=rules.functional_redundancy.output.summ
    output: f"{RESULTS}/figures/redundancy/redundancy.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_redundancy.py --config config/config.yaml "
        f"--redundancy-dir {RED} --out {RESULTS}/figures/redundancy/redundancy"

# ---- Western vs non-Western turnover ----------------------------------------
rule population_turnover:
    input:
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        anchor=f"{RESULTS}/03_profiles/prevalence_ko.parquet",
    output:
        turn=f"{POP}/population_turnover.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/population_turnover.py --config config/config.yaml "
        f"--species-table {{input.species}} --profiles-dir {RESULTS}/03_profiles --out-dir {POP}"

rule fig_population:
    input: turn=rules.population_turnover.output.turn
    output: f"{RESULTS}/figures/population/population.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_population.py --config config/config.yaml "
        f"--turnover {{input.turn}} --substitution-dir {POP} "
        f"--out {RESULTS}/figures/population/population"

# ---- combined report --------------------------------------------------------
rule redundancy_report:
    input:
        red=rules.functional_redundancy.output.summ,
        turn=rules.population_turnover.output.turn,
        fred=rules.fig_redundancy.output,
        fpop=rules.fig_population.output,
        # the western_vs_nonwestern divergent functions (only when enabled)
        div=(f"{RESULTS}/05_diff/combined/{POP_CONTRAST}_signatures_all.tsv"
             if POP_ENABLED else []),
    output:
        html=f"{RESULTS}/report/redundancy_report.html",
        md=f"{RESULTS}/report/redundancy_report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_redundancy_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/redundancy_report"
