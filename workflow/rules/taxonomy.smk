# =============================================================================
# taxonomy.smk -- species/taxonomy-by-niche analysis (the current focus).
#
# Depends only on the metadata (ingest, species_units) and the GTDB scaffold;
# it does NOT need the functional annotations, so this whole stage runs now.
# Upstream files are referenced by PATH (not rules.*) so include order does not
# matter. The default target `rule all` (in the Snakefile) builds the taxonomy
# report, which pulls every rule below.
# =============================================================================

SAMPLES = f"{RESULTS}/00_ingest/samples.parquet"
SPTAB   = f"{RESULTS}/00_ingest/species_table.tsv"
TREE    = f"{RESULTS}/04_phylo/species_tree.nwk"
TIPMAP  = f"{RESULTS}/04_phylo/tip_map.tsv"
TX      = f"{RESULTS}/taxonomy"
TXFIG   = f"{RESULTS}/figures/taxonomy"
PLOTRANK = "family"   # rank shown in enrichment/indicator/specificity figures

# ---- analyses ---------------------------------------------------------------
rule taxon_specificity:
    input: samples=SAMPLES, species=SPTAB
    output: directory(f"{TX}/02_specificity")
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/taxon_niche_specificity.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} --out-dir {{output}}"

rule taxon_enrichment:
    input: samples=SAMPLES, species=SPTAB
    output: directory(f"{TX}/03_enrichment")
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/taxon_enrichment.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} --out-dir {{output}}"

rule indicator_species:
    input: species=SPTAB
    output: directory(f"{TX}/03_enrichment_indicator")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/indicator_species.R --config config/config.yaml "
        f"--species-table {{input.species}} --out-dir {{output}}"

rule overlap_nullmodel:
    input: samples=SAMPLES
    output:
        nm=f"{TX}/04_overlap/overlap_nullmodel.tsv",
        obs=f"{TX}/04_overlap/overlap_observed.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/overlap_nullmodel.py --config config/config.yaml "
        f"--samples {{input.samples}} --out-prefix {TX}/04_overlap/overlap"

rule taxonomic_novelty:
    input: samples=SAMPLES, species=SPTAB
    output:
        byniche=f"{TX}/05_novelty/novelty_by_niche.tsv",
        byrank=f"{TX}/05_novelty/novelty_by_rank.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/taxonomic_novelty.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} "
        f"--out-prefix {TX}/05_novelty/novelty"

rule phylo_community:
    input: species=SPTAB, tree=TREE, tipmap=TIPMAP
    output: directory(f"{TX}/06_phylo_community")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/phylo_community.R --config config/config.yaml "
        f"--species-table {{input.species}} --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --out-dir {{output}}"

rule beta_partition:
    input: species=SPTAB
    output: directory(f"{TX}/07_beta")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/beta_partition.R --config config/config.yaml "
        f"--species-table {{input.species}} --out-dir {{output}}"

rule host_resolved:
    input: samples=SAMPLES
    output:
        hs=f"{TX}/08_host/host_host_summary.tsv",
        hp=f"{TX}/08_host/host_host_phylum_composition.tsv",
        ww=f"{TX}/08_host/host_animal_with_without_mouse.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/host_resolved.py --config config/config.yaml "
        f"--samples {{input.samples}} --out-prefix {TX}/08_host/host"

# ---- taxonomy figures (PNG + editable-text SVG) -----------------------------
rule fig_specificity:
    input: spec=rules.taxon_specificity.output
    output: f"{TXFIG}/fig_specificity.png", f"{TXFIG}/fig_specificity.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_specificity.py --config config/config.yaml "
        f"--spec-dir {{input.spec}} --rank {PLOTRANK} --out {TXFIG}/fig_specificity"

rule fig_taxa_signatures:
    input:
        enrich=rules.taxon_enrichment.output,
        indic=rules.indicator_species.output,
    output: f"{TXFIG}/fig_taxa_signatures.png", f"{TXFIG}/fig_taxa_signatures.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_taxa_signatures.py --config config/config.yaml "
        f"--enrichment {{input.enrich}}/enrichment_{PLOTRANK}.tsv "
        f"--indicator {{input.indic}}/indicator_{PLOTRANK}.tsv "
        f"--rank {PLOTRANK} --out {TXFIG}/fig_taxa_signatures"

rule fig_phylo_community:
    input: pc=rules.phylo_community.output
    output: f"{TXFIG}/fig_phylo_community.png", f"{TXFIG}/fig_phylo_community.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_phylo_community.py --config config/config.yaml "
        f"--phylo-community {{input.pc}}/phylo_community.tsv --out {TXFIG}/fig_phylo_community"

rule fig_beta_overlap:
    input:
        beta=rules.beta_partition.output,
        ov=rules.overlap_nullmodel.output.nm,
    output: f"{TXFIG}/fig_beta_overlap.png", f"{TXFIG}/fig_beta_overlap.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_beta_overlap.py --config config/config.yaml "
        f"--beta-pairwise {{input.beta}}/beta_pairwise.tsv "
        f"--overlap-null {{input.ov}} --out {TXFIG}/fig_beta_overlap"

rule fig_novelty:
    input:
        byniche=rules.taxonomic_novelty.output.byniche,
        byrank=rules.taxonomic_novelty.output.byrank,
    output: f"{TXFIG}/fig_novelty.png", f"{TXFIG}/fig_novelty.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_novelty.py --config config/config.yaml "
        f"--by-niche {{input.byniche}} --by-rank {{input.byrank}} --out {TXFIG}/fig_novelty"

rule fig_host_resolved:
    input:
        hs=rules.host_resolved.output.hs,
        hp=rules.host_resolved.output.hp,
    output: f"{TXFIG}/fig_host_resolved.png", f"{TXFIG}/fig_host_resolved.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_host_resolved.py --config config/config.yaml "
        f"--host-summary {{input.hs}} --host-phylum {{input.hp}} --out {TXFIG}/fig_host_resolved"

rule fig_tree_niche:
    input: species=SPTAB, tree=TREE, tipmap=TIPMAP
    output: f"{TXFIG}/fig_tree_niche.png", f"{TXFIG}/fig_tree_niche.svg"
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/fig_tree_niche.R --config config/config.yaml --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --species-table {{input.species}} "
        f"--out {TXFIG}/fig_tree_niche"

# ---- taxonomy report (default target dependency) ----------------------------
rule taxonomy_report:
    input:
        # ecology figures (from figures.smk, routed under figures/taxonomy/)
        f"{TXFIG}/fig1_dataset_overview.png",
        f"{TXFIG}/fig2_diversity_composition.png",
        f"{TXFIG}/fig3_occupancy.png",
        # taxonomy figures
        rules.fig_specificity.output,
        rules.fig_taxa_signatures.output,
        rules.fig_phylo_community.output,
        rules.fig_beta_overlap.output,
        rules.fig_novelty.output,
        rules.fig_host_resolved.output,
        rules.fig_tree_niche.output,
        # tables for the report body
        f"{RESULTS}/01_ecology/div_richness_estimators.tsv",
        rules.taxon_specificity.output,
        rules.overlap_nullmodel.output.nm,
        rules.phylo_community.output,
        rules.beta_partition.output,
        rules.taxonomic_novelty.output.byniche,
        rules.host_resolved.output.ww,
        f"{RESULTS}/04_phylo/scaffold_coverage.json",
    output:
        html=f"{RESULTS}/report/taxonomy_report.html",
        md=f"{RESULTS}/report/taxonomy_report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_taxonomy_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/taxonomy_report"
