# =============================================================================
# scoary.smk -- Scoary2 pan-GWAS as a 4th association method (part of
# `functional_all`; convenience target `scoary_all`).
#
# (1) Species-level: reuse differential_prep's groups + the species-level
#     functional presence matrices; the GTDB species tree supplies the pairwise-
#     comparisons correction. Results fold into the consensus (see comparative.smk).
# (2) Clade-stratified: Panaroo per genus on species reps -> Scoary2 on real
#     ortholog families within the clade -> aggregate.
# =============================================================================

SCOARY = config["scoary"]
SCOARY_CONTRASTS = list(SCOARY["contrasts"])
if config.get("population", {}).get("enabled"):        # 4th method for western_vs_nonwestern too
    _pc = config["population"]["contrast"]
    if _pc not in SCOARY_CONTRASTS:
        SCOARY_CONTRASTS.append(_pc)
SCOARY_LAYERS = SCOARY["layers"]
SC = f"{RESULTS}/scoary"
TREE_S = f"{RESULTS}/04_phylo/species_tree.nwk"
TIPMAP_S = f"{RESULTS}/04_phylo/tip_map.tsv"

rule scoary_all:
    input:
        expand(f"{SC}/species/{{layer}}/{{contrast}}/scoary.tsv",
               layer=SCOARY_LAYERS, contrast=SCOARY_CONTRASTS),
        f"{SC}/clade/clade_scoary_summary.tsv",
        f"{RESULTS}/figures/scoary/clade_summary.png",

# ---- species-level ----------------------------------------------------------
rule scoary_prepare:
    input:
        analysis=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_analysis_species.tsv",
        presence=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_presence.parquet",
        tree=TREE_S, tipmap=TIPMAP_S,
    output:
        genes=f"{SC}/species/{{layer}}/{{contrast}}/input/genes.tsv",
        traits=f"{SC}/species/{{layer}}/{{contrast}}/input/traits.tsv",
        tree=f"{SC}/species/{{layer}}/{{contrast}}/input/tree.nwk",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_prepare.py --config config/config.yaml "
        f"--analysis {{input.analysis}} --presence {{input.presence}} "
        f"--tip-map {{input.tipmap}} --tree {{input.tree}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} "
        f"--out-dir {SC}/species/{{wildcards.layer}}/{{wildcards.contrast}}/input"

rule scoary_run:
    input:
        genes=rules.scoary_prepare.output.genes,
        traits=rules.scoary_prepare.output.traits,
        tree=rules.scoary_prepare.output.tree,
    output:
        results=f"{SC}/species/{{layer}}/{{contrast}}/out/traits/{{contrast}}/results.tsv",
    params:
        outdir=f"{SC}/species/{{layer}}/{{contrast}}/out",
        mt=SCOARY["multiple_testing"], np=SCOARY["n_permut"],
        wc=SCOARY["worst_cutoff"],
        mg=("--max_genes %d" % SCOARY["max_genes"]) if SCOARY["max_genes"] else "",
        seed=config["seed"],
    threads: 8
    conda: "../envs/scoary.yaml"
    shell:
        r"""
        set -euo pipefail
        TREEARG=""; [ -s {input.tree} ] && TREEARG="--newicktree {input.tree}"
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        scoary2 --genes {input.genes} --gene_data_type 'gene-count:\t' \
                --traits {input.traits} --trait_data_type 'binary:\t' \
                $TREEARG --outdir {params.outdir} \
                --multiple_testing {params.mt} --n_permut {params.np} \
                --n_cpus {threads} --worst_cutoff {params.wc} {params.mg} \
                --random_state {params.seed} \
        || echo "Scoary2 produced no significant traits for {wildcards.layer}/{wildcards.contrast}"
        # ensure the expected results file exists even if Scoary2 filtered the trait out
        mkdir -p $(dirname {output.results}); [ -f {output.results} ] || touch {output.results}
        """

rule scoary_parse:
    input: results=rules.scoary_run.output.results
    output: f"{SC}/species/{{layer}}/{{contrast}}/scoary.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_parse.py --config config/config.yaml --results {{input.results}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} --out {{output}}"

# ---- clade-stratified -------------------------------------------------------
checkpoint scoary_select_clades:
    input:
        species=f"{RESULTS}/00_ingest/species_table.tsv",
        reps=f"{RESULTS}/00_ingest/representatives.tsv",
    output: directory(f"{SC}/clade/selection")
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_select_clades.py --config config/config.yaml "
        f"--species-table {{input.species}} --representatives {{input.reps}} "
        f"--out-dir {SC}/clade/selection"

def _clade_ids(wildcards):
    ck = checkpoints.scoary_select_clades.get().output[0]
    return sorted(glob_wildcards(os.path.join(ck, "clades", "{cid}",
                                              "species_genomes.tsv")).cid)

rule clade_panaroo:
    input:
        sg=f"{SC}/clade/selection/clades/{{cid}}/species_genomes.tsv",
        paths=ANNOT_PATHS,
    output:
        rtab=f"{SC}/clade/clades/{{cid}}/panaroo/gene_presence_absence.Rtab",
        csv=f"{SC}/clade/clades/{{cid}}/panaroo/gene_presence_absence_roary.csv",
    params:
        wd=f"{SC}/clade/clades/{{cid}}/panaroo",
        ct=SCOARY["clade"]["core_threshold"],
        lookup=lambda wc, input: annot_lookup(f"{SC}/clade/clades/{wc.cid}/panaroo/ids.txt",
                                              "prokka_gff", input.paths),
    threads: 8
    conda: "../envs/transition.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        # genome ids are column 2 of species_genomes.tsv; the symlink farm of
        # <genome>.gff keeps Panaroo's sample names equal to the genome ids
        tail -n +2 {input.sg} | cut -f2 > {params.wd}/ids.txt
        {params.lookup} > {params.wd}/gff_paths.tsv
        bash {SCRIPTS}/sh/symlink_farm.sh {params.wd}/gff_paths.tsv {params.wd}/gff gff
        cp {params.wd}/gff/list.txt {params.wd}/gff_list.txt
        panaroo -i $(cat {params.wd}/gff_list.txt) -o {params.wd} \
          --clean-mode moderate --core_threshold {params.ct} -t {threads}
        """

rule clade_scoary_prepare:
    input:
        sg=f"{SC}/clade/selection/clades/{{cid}}/species_genomes.tsv",
        rtab=rules.clade_panaroo.output.rtab,
        tree=TREE_S, tipmap=TIPMAP_S,
    output:
        genes=f"{SC}/clade/clades/{{cid}}/input/{{ccontrast}}/genes.tsv",
        traits=f"{SC}/clade/clades/{{cid}}/input/{{ccontrast}}/traits.tsv",
        tree=f"{SC}/clade/clades/{{cid}}/input/{{ccontrast}}/tree.nwk",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_clade_prepare.py --config config/config.yaml "
        f"--species-genomes {{input.sg}} --rtab {{input.rtab}} "
        f"--tip-map {{input.tipmap}} --tree {{input.tree}} "
        f"--contrast {{wildcards.ccontrast}} "
        f"--out-dir {SC}/clade/clades/{{wildcards.cid}}/input/{{wildcards.ccontrast}}"

rule clade_scoary_run:
    input:
        genes=rules.clade_scoary_prepare.output.genes,
        traits=rules.clade_scoary_prepare.output.traits,
        tree=rules.clade_scoary_prepare.output.tree,
    output:
        results=f"{SC}/clade/clades/{{cid}}/out/{{ccontrast}}/traits/{{ccontrast}}/results.tsv",
    params:
        outdir=f"{SC}/clade/clades/{{cid}}/out/{{ccontrast}}",
        mt=SCOARY["multiple_testing"], np=SCOARY["n_permut"], seed=config["seed"],
    threads: 4
    conda: "../envs/scoary.yaml"
    shell:
        r"""
        set -euo pipefail
        TREEARG=""; [ -s {input.tree} ] && TREEARG="--newicktree {input.tree}"
        rm -rf {params.outdir}; mkdir -p {params.outdir}
        scoary2 --genes {input.genes} --gene_data_type 'gene-count:\t' \
                --traits {input.traits} --trait_data_type 'binary:\t' \
                $TREEARG --outdir {params.outdir} --multiple_testing {params.mt} \
                --n_permut {params.np} --n_cpus {threads} --random_state {params.seed} \
        || echo "no sig traits {wildcards.cid}/{wildcards.ccontrast}"
        mkdir -p $(dirname {output.results}); [ -f {output.results} ] || touch {output.results}
        """

rule clade_scoary_parse:
    input: results=rules.clade_scoary_run.output.results
    output: f"{SC}/clade/clades/{{cid}}/scoary/{{ccontrast}}.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_parse.py --config config/config.yaml --results {{input.results}} "
        f"--layer pangenome --contrast {{wildcards.ccontrast}} --out {{output}}"

def _clade_scoary_outputs(wildcards):
    cids = _clade_ids(wildcards)
    return [f"{SC}/clade/clades/{c}/scoary/{cc}.tsv"
            for c in cids for cc in SCOARY["clade"]["contrasts"]]

rule clade_scoary_aggregate:
    input: _clade_scoary_outputs
    output:
        allt=f"{SC}/clade/clade_scoary_all.tsv",
        summ=f"{SC}/clade/clade_scoary_summary.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/scoary_clade_aggregate.py --root {SC}/clade "
        f"--out-prefix {SC}/clade/clade_scoary"

# ---- figures ----------------------------------------------------------------
rule fig_scoary_concordance:
    input:
        sig=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.tsv",
    output: f"{RESULTS}/figures/scoary/{{layer}}_{{contrast}}_concordance.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_scoary_concordance.py --config config/config.yaml "
        f"--signatures {{input.sig}} --layer {{wildcards.layer}} "
        f"--contrast {{wildcards.contrast}} "
        f"--out {RESULTS}/figures/scoary/{{wildcards.layer}}_{{wildcards.contrast}}_concordance"

rule fig_scoary_clade:
    input: summ=rules.clade_scoary_aggregate.output.summ
    output: f"{RESULTS}/figures/scoary/clade_summary.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_scoary_clade.py --config config/config.yaml "
        f"--summary {{input.summ}} --out {RESULTS}/figures/scoary/clade_summary"
