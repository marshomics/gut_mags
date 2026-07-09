# =============================================================================
# report.smk -- the final target. Depends on every figure and summary so that
# requesting report.html builds the whole DAG, then assembles an HTML/markdown
# report with the headline numbers, confounder-control recap, signature tables,
# permutation FDR, positive-control check and a figure gallery.
# =============================================================================

def _report_inputs(wildcards):
    f = []
    # ecology + overview figures
    f += [f"{FIG}/taxonomy/fig1_dataset_overview.png",
          f"{FIG}/taxonomy/fig2_diversity_composition.png",
          f"{FIG}/taxonomy/fig3_occupancy.png",
          f"{FIG}/fig_function_landscape.png",
          f"{FIG}/fig_ordination_varpart.png"]
    # per layer x contrast
    for c in CONTRASTS:
        f += [f"{FIG}/heatmap/{c}_signature_heatmap.png",
              f"{FIG}/tree/{c}_tree_annotated.png",
              f"{FIG}/tree/{c}_convergence.png",
              f"{RESULTS}/05_diff/combined/{c}_signatures_all.tsv",
              f"{RESULTS}/06_comparative/{c}_signal",
              f"{RESULTS}/06_comparative/{c}_ancestral"]
        for l in DIFF_LAYERS:
            f += [f"{FIG}/volcano/{l}_{c}.png",
                  f"{RESULTS}/05_diff/{l}/{c}_signatures.json",
                  f"{RESULTS}/05_diff/{l}/{c}_permnull.json"]
    # Scoary2 (4th method): concordance figures + clade-stratified summary
    sc = config.get("scoary", {})
    if sc.get("enabled"):
        for c in sc.get("contrasts", []):
            for l in sc.get("layers", []):
                f.append(f"{FIG}/scoary/{l}_{c}_concordance.png")
        f += [f"{FIG}/scoary/clade_summary.png",
              f"{RESULTS}/scoary/clade/clade_scoary_summary.tsv"]
    # functional enrichment of niche signatures
    en = config.get("enrichment", {})
    if en.get("enabled"):
        f += [f"{RESULTS}/08_enrichment/enrichment_top_by_contrast.tsv",
              f"{FIG}/enrichment/themes_heatmap.png"]
        f += [f"{FIG}/enrichment/dotplot_{c}.png" for c in en.get("contrasts", [])]
    # summaries
    f += [rules.ingest.output.report,
          rules.scaffold.output.cov,
          rules.species_profiles.output.ko_conc,
          rules.pgls.output]
    return f

rule report:
    input:
        _report_inputs,
    output:
        html=f"{RESULTS}/report/report.html",
        md=f"{RESULTS}/report/report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/report"
