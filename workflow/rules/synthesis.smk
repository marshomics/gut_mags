# =============================================================================
# synthesis.smk -- capstone integrating every stage to answer the six questions
# (target `snakemake synthesis_all`). Adds the driver analyses (adaptation mode,
# diversification rate, curated ecological pressures, HyPhy selection) and builds
# the integrated human-specific catalogues, master figure and report.
# Depends on taxonomy + functional + enrichment + Scoary2 clade outputs.
# =============================================================================

SY = f"{RESULTS}/09_synthesis"
SYN = config["synthesis"]
FOCAL_CONTRAST = SYN["focal_contrast"]
ECO_CONTRASTS = SYN["ecological_pressures"]["contrasts"]
COMBINED_SIG = f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv"

rule synthesis_all:
    input:
        f"{RESULTS}/report/synthesis_report.html",

# ---- driver: ecological pressures (per contrast) ----------------------------
rule ecological_pressure:
    input:
        sig=f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv",
        pmap=SYN["ecological_pressures"]["map"],
    output: f"{SY}/ecological/ecological_pressure_{{contrast}}.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/ecological_pressure_test.py --config config/config.yaml "
        f"--signatures {{input.sig}} --pressure-map {{input.pmap}} "
        f"--contrast {{wildcards.contrast}} --out {{output}}"

# ---- driver: diversification rate by niche ----------------------------------
rule diversification_rate:
    input:
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        tree=f"{RESULTS}/04_phylo/species_tree.nwk",
        tipmap=f"{RESULTS}/04_phylo/tip_map.tsv",
    output:
        dr=f"{SY}/diversification/diversification_dr.tsv",
        test=f"{SY}/diversification/diversification_test.json",
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/diversification_rate.R --config config/config.yaml "
        f"--species-table {{input.species}} --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --out-dir {SY}/diversification"

# ---- driver: sequence-level selection (HyPhy), per clade --------------------
def _sel_clade_ids(wildcards):
    ck = checkpoints.scoary_select_clades.get().output[0]
    return sorted(glob_wildcards(os.path.join(ck, "clades", "{cid}",
                                              "species_genomes.tsv")).cid)

rule selection_prepare:
    input:
        gpa=f"{RESULTS}/scoary/clade/clades/{{cid}}/panaroo/gene_presence_absence_roary.csv",
        scoary=f"{RESULTS}/scoary/clade/clades/{{cid}}/scoary/human_vs_rest.tsv",
        sg=f"{RESULTS}/scoary/clade/selection/clades/{{cid}}/species_genomes.tsv",
        paths=ANNOT_PATHS,
    output: directory(f"{SY}/selection/clades/{{cid}}/prep")
    conda: "../envs/selection.yaml"
    shell:
        f"{PY}/selection_prepare.py --config config/config.yaml "
        f"--gpa-csv {{input.gpa}} --scoary {{input.scoary}} "
        f"--species-genomes {{input.sg}} --annotation-paths {{input.paths}} "
        f"--out-dir {{output}}"

rule selection_hyphy:
    input: prep=rules.selection_prepare.output
    output: sentinel=f"{SY}/selection/clades/{{cid}}/hyphy.done"
    threads: 4
    conda: "../envs/selection.yaml"
    shell:
        r"""
        set -euo pipefail
        for fam in {input.prep}/families/*/; do
          aln="$fam/codon.fasta"; tre="$fam/tree.nwk"
          [ -s "$aln" ] && [ -s "$tre" ] || continue
          hyphy busted --alignment "$aln" --tree "$tre" --output "$fam/BUSTED.json" || true
          if grep -q "Foreground" "$tre"; then
            hyphy relax --alignment "$aln" --tree "$tre" --test Foreground \
                  --output "$fam/RELAX.json" || true
          fi
        done
        touch {output.sentinel}
        """

rule selection_parse:
    input:
        sentinel=rules.selection_hyphy.output.sentinel,
        prep=rules.selection_prepare.output,
    output: f"{SY}/selection/clades/{{cid}}/selection.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/selection_parse.py --config config/config.yaml "
        f"--sel-dir {{input.prep}} --clade {{wildcards.cid}} --out {{output}}"

def _selection_tables(wildcards):
    if not SYN["selection"]["enabled"]:
        return []
    return [f"{SY}/selection/clades/{c}/selection.tsv" for c in _sel_clade_ids(wildcards)]

rule selection_aggregate:
    input: _selection_tables
    output:
        allt=f"{SY}/selection/selection_all.tsv",
        summ=f"{SY}/selection/selection_summary.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/selection_aggregate.py --inputs {{input}} --out-prefix {SY}/selection/selection "
        f"|| (mkdir -p {SY}/selection && : > {{output.allt}} && : > {{output.summ}})"

# ---- driver: adaptation mode ------------------------------------------------
def _selection_for_adapt(wildcards):
    return rules.selection_aggregate.output.summ if SYN["selection"]["enabled"] else []

rule adaptation_mode:
    input:
        sig=f"{RESULTS}/05_diff/combined/{FOCAL_CONTRAST}_signatures_all.tsv",
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        pgls=f"{RESULTS}/06_comparative/pgls_results.tsv",
        selection=_selection_for_adapt,
    output:
        feats=f"{SY}/adaptation/adaptation_mode_features.tsv",
        summ=f"{SY}/adaptation/adaptation_mode_summary.tsv",
        drivers=f"{SY}/adaptation/adaptation_mode_drivers.json",
    params:
        sel_arg=lambda wc, input: f"--selection {input.selection}" if input.selection else "",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/adaptation_mode.py --config config/config.yaml "
        f"--signatures {{input.sig}} --species-table {{input.species}} "
        f"--profiles-dir {RESULTS}/03_profiles --pgls {{input.pgls}} {{params.sel_arg}} "
        f"--out-prefix {SY}/adaptation/adaptation_mode"

# ---- integrated catalogues --------------------------------------------------
rule synthesis_catalogues:
    input:
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        indicator=f"{RESULTS}/taxonomy/03_enrichment_indicator",
        sig=f"{RESULTS}/05_diff/combined/{FOCAL_CONTRAST}_signatures_all.tsv",
        adapt=rules.adaptation_mode.output.feats,
        enr=f"{RESULTS}/08_enrichment/enrichment_top_by_contrast.tsv",
        eco=f"{SY}/ecological/ecological_pressure_{FOCAL_CONTRAST}.tsv",
    output:
        counts=f"{SY}/catalogues/catalogue_counts.json",
        species_cat=f"{SY}/catalogues/human_specific_species.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/synthesis_catalogues.py --config config/config.yaml "
        f"--species-table {{input.species}} --indicator-dir {{input.indicator}} "
        f"--signatures {{input.sig}} --adaptation {{input.adapt}} "
        f"--enrichment-top {{input.enr}} --ecological {{input.eco}} "
        f"--out-dir {SY}/catalogues"

# ---- master figure + report -------------------------------------------------
rule fig_synthesis_master:
    input:
        counts=rules.synthesis_catalogues.output.counts,
        adapt=rules.adaptation_mode.output.summ,
        eco=f"{SY}/ecological/ecological_pressure_{FOCAL_CONTRAST}.tsv",
        dr=rules.diversification_rate.output.dr,
    output: f"{RESULTS}/figures/synthesis/master.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_synthesis_master.py --config config/config.yaml "
        f"--counts {{input.counts}} --adaptation-summary {{input.adapt}} "
        f"--ecological {{input.eco}} --diversification {{input.dr}} "
        f"--out {RESULTS}/figures/synthesis/master"

rule synthesis_report:
    input:
        counts=rules.synthesis_catalogues.output.counts,
        master=rules.fig_synthesis_master.output,
        seln=rules.selection_aggregate.output.summ,
        eco=expand(f"{SY}/ecological/ecological_pressure_{{contrast}}.tsv", contrast=ECO_CONTRASTS),
        divtest=rules.diversification_rate.output.test,
    output:
        html=f"{RESULTS}/report/synthesis_report.html",
        md=f"{RESULTS}/report/synthesis_report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_synthesis_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/synthesis_report"
