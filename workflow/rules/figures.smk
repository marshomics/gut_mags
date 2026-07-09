# =============================================================================
# figures.smk -- every figure as PNG (raster) + SVG (editable text).
# Python figures use plotting_theme.save(); R figures use save_pub()/save_base().
# =============================================================================

FIG = f"{RESULTS}/figures"

rule fig_dataset_overview:
    input:
        samples=rules.ingest.output.parquet,
        strains=rules.species_units.output.strains,
    output: f"{FIG}/taxonomy/fig1_dataset_overview.png", f"{FIG}/taxonomy/fig1_dataset_overview.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_dataset_overview.py --config config/config.yaml "
        f"--samples {{input.samples}} --strains-per-species {{input.strains}} "
        f"--out {FIG}/taxonomy/fig1_dataset_overview"

rule fig_diversity:
    input:
        rare=rules.diversity.output.rare,
        hill=rules.diversity.output.hill,
        comp=rules.taxonomy_composition.output,
    output: f"{FIG}/taxonomy/fig2_diversity_composition.png", f"{FIG}/taxonomy/fig2_diversity_composition.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_diversity_composition.py --config config/config.yaml "
        f"--rarefaction {{input.rare}} --hill {{input.hill}} "
        f"--composition-phylum {{input.comp}}/composition_phylum.tsv "
        f"--out {FIG}/taxonomy/fig2_diversity_composition"

rule fig_occupancy:
    input:
        upset=rules.niche_breadth.output.upset,
        breadth=rules.niche_breadth.output.breadth,
    output: f"{FIG}/taxonomy/fig3_occupancy.png", f"{FIG}/taxonomy/fig3_occupancy.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_occupancy_upset.py --config config/config.yaml "
        f"--upset {{input.upset}} --breadth {{input.breadth}} "
        f"--out {FIG}/taxonomy/fig3_occupancy"

rule fig_volcano:
    input:
        sig=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.tsv",
    output:
        f"{FIG}/volcano/{{layer}}_{{contrast}}.png",
        f"{FIG}/volcano/{{layer}}_{{contrast}}.svg",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_volcano_concordance.py --config config/config.yaml "
        f"--signatures {{input.sig}} --layer {{wildcards.layer}} "
        f"--contrast {{wildcards.contrast}} "
        f"--out {FIG}/volcano/{{wildcards.layer}}_{{wildcards.contrast}}"

rule fig_signature_heatmap:
    input:
        sig=f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv",
        species=rules.species_units.output.table,
        anchor=rules.species_profiles.output.ko,
    output:
        f"{FIG}/heatmap/{{contrast}}_signature_heatmap.png",
        f"{FIG}/heatmap/{{contrast}}_signature_heatmap.svg",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_signature_heatmap.py --config config/config.yaml "
        f"--signatures-combined {{input.sig}} --profiles-dir {RESULTS}/03_profiles "
        f"--species-table {{input.species}} "
        f"--out {FIG}/heatmap/{{wildcards.contrast}}_signature_heatmap"

rule fig_function_landscape:
    input:
        traits=rules.species_traits.output,
        pgls=rules.pgls.output,
    output: f"{FIG}/fig_function_landscape.png", f"{FIG}/fig_function_landscape.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_function_landscape.py --config config/config.yaml "
        f"--traits {{input.traits}} --pgls {{input.pgls}} "
        f"--out {FIG}/fig_function_landscape"

rule fig_ordination:
    input:
        ord=rules.ordination.output,
    output: f"{FIG}/fig_ordination_varpart.png", f"{FIG}/fig_ordination_varpart.svg"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_ordination_varpart.py --config config/config.yaml "
        f"--coords {{input.ord}}/pcoa_coords.tsv --varexp {{input.ord}}/pcoa_varexp.tsv "
        f"--permanova {{input.ord}}/permanova.tsv --varpart {{input.ord}}/varpart.tsv "
        f"--out {FIG}/fig_ordination_varpart"

rule fig_tree_annotated:
    input:
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
        species=rules.species_units.output.table,
        sig=f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv",
        presence=rules.combine_presence.output,
    output:
        f"{FIG}/tree/{{contrast}}_tree_annotated.png",
        f"{FIG}/tree/{{contrast}}_tree_annotated.svg",
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/fig_tree_annotated.R --config config/config.yaml --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --species-table {{input.species}} "
        f"--signatures {{input.sig}} --presence {{input.presence}} "
        f"--out {FIG}/tree/{{wildcards.contrast}}_tree_annotated"

rule fig_convergence:
    input:
        anc=f"{RESULTS}/06_comparative/{{contrast}}_ancestral",
    output:
        f"{FIG}/tree/{{contrast}}_convergence.png",
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/fig_convergence.R --config config/config.yaml "
        f"--convergence {{input.anc}}/convergence_summary.tsv "
        f"--simmap-dir {{input.anc}}/simmaps "
        f"--out {FIG}/tree/{{wildcards.contrast}}_convergence"
