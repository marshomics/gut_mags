# =============================================================================
# enrichment.smk -- functional enrichment of the niche-signature features
# (part of `functional_all`; convenience target `enrichment_all`).
#
# Per feature layer, category memberships (gene sets) are built once, then each
# layer x contrast consensus signature set is tested by over-representation (ORA,
# hypergeometric with the tested features as background) and preranked GSEA
# (fgsea). Results are combined (agreement flagged), aggregated across layers,
# and plotted as a dot plot per contrast and a category-by-contrast heatmap.
# =============================================================================

EN = config["enrichment"]
ENRICH_LAYERS = list(EN["systems"].keys())            # ko, pfam, cog, cazyme, bgc, amr
ENRICH_CONTRASTS = list(EN["contrasts"])
if config.get("population", {}).get("enabled"):        # add western_vs_nonwestern when on
    _pc = config["population"]["contrast"]
    if _pc not in ENRICH_CONTRASTS:
        ENRICH_CONTRASTS.append(_pc)
ENR = f"{RESULTS}/08_enrichment"

rule enrichment_all:
    input:
        f"{ENR}/enrichment_all.tsv",
        expand(f"{RESULTS}/figures/enrichment/dotplot_{{contrast}}.png",
               contrast=ENRICH_CONTRASTS),
        f"{RESULTS}/figures/enrichment/themes_heatmap.png",

rule build_genesets:
    input:
        prevalence=f"{RESULTS}/03_profiles/prevalence_{{layer}}.parquet",
    output:
        f"{ENR}/genesets/genesets_{{layer}}.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/build_genesets.py --config config/config.yaml "
        f"--layer {{wildcards.layer}} --prevalence {{input.prevalence}} --out {{output}}"

rule ora_enrichment:
    input:
        sig=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.tsv",
        gs=rules.build_genesets.output,
    output: f"{ENR}/{{layer}}/{{contrast}}_ora.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/ora_enrichment.py --config config/config.yaml "
        f"--signatures {{input.sig}} --genesets {{input.gs}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} --out {{output}}"

rule gsea_enrichment:
    input:
        sig=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.tsv",
        gs=rules.build_genesets.output,
    output: f"{ENR}/{{layer}}/{{contrast}}_gsea.tsv"
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/gsea_enrichment.R --config config/config.yaml "
        f"--signatures {{input.sig}} --genesets {{input.gs}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} --out {{output}}"

rule enrichment_combine:
    input:
        ora=rules.ora_enrichment.output,
        gsea=rules.gsea_enrichment.output,
    output: f"{ENR}/{{layer}}/{{contrast}}_combined.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/enrichment_combine.py --config config/config.yaml "
        f"--ora {{input.ora}} --gsea {{input.gsea}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} --out {{output}}"

rule enrichment_aggregate:
    input:
        expand(f"{ENR}/{{layer}}/{{contrast}}_combined.tsv",
               layer=ENRICH_LAYERS, contrast=ENRICH_CONTRASTS),
    output:
        allt=f"{ENR}/enrichment_all.tsv",
        top=f"{ENR}/enrichment_top_by_contrast.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/enrichment_aggregate.py --inputs {{input}} --out-prefix {ENR}/enrichment"

rule fig_enrichment_dotplot:
    input: enr=rules.enrichment_aggregate.output.allt
    output: f"{RESULTS}/figures/enrichment/dotplot_{{contrast}}.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_enrichment_dotplot.py --config config/config.yaml "
        f"--enrichment {{input.enr}} --contrast {{wildcards.contrast}} "
        f"--out {RESULTS}/figures/enrichment/dotplot_{{wildcards.contrast}}"

rule fig_enrichment_heatmap:
    input: enr=rules.enrichment_aggregate.output.allt
    output: f"{RESULTS}/figures/enrichment/themes_heatmap.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_enrichment_heatmap.py --config config/config.yaml "
        f"--enrichment {{input.enr}} --out {RESULTS}/figures/enrichment/themes_heatmap"
