# =============================================================================
# transition.smk -- within-species niche-transition analysis (separate target:
# `snakemake transition_all`). For each qualifying species (>=2 niches, enough
# near-complete strains each): dereplicate -> Panaroo core -> Gubbins
# recombination masking -> IQ-TREE -> ancestral-state polarization + per-niche
# SFS/diversity + directionality + accessory gene gains -> verdict; then a
# cross-species meta-analysis of transition direction.
#
# Heavy and per-species, so it is NOT in the default target. Needs the metadata,
# the GTDB scaffold (for species/outgroup taxonomy via species_table), the
# assemblies and the Prokka GFFs.
# =============================================================================

TR = f"{RESULTS}/transition"
TRFIG = f"{RESULTS}/figures/transition"
TC = config["transition"]

# Mash and Panaroo both name each sample after its file's basename, and every
# downstream step (dereplication, the tree's tip labels, the niche map) keys on
# the genome id. With {genome} templates those coincide; with a manifest of
# arbitrary paths they do not. So each rule first builds a symlink farm of
# canonically named <genome>.fna / <genome>.gff files (scripts/sh/symlink_farm.sh)
# and points the tool at that. Nothing downstream needs to know where the real
# file lives.
FARM = f"bash {SCRIPTS}/sh/symlink_farm.sh"

rule transition_all:
    input:
        f"{RESULTS}/report/transition_report.html",

# ---- checkpoint: choose species + assemble inputs ---------------------------
checkpoint select_transition_species:
    input:
        samples=f"{RESULTS}/00_ingest/samples.parquet",
        species=f"{RESULTS}/00_ingest/species_table.tsv",
    output:
        directory(f"{TR}/selection"),
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/transition_select.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} "
        f"--out-dir {TR}/selection"


def _species_ids(wildcards):
    ck = checkpoints.select_transition_species.get().output[0]
    return sorted(glob_wildcards(os.path.join(ck, "species", "{sid}",
                                              "candidate_genomes.tsv")).sid)

# ---- per-species steps ------------------------------------------------------
rule mash_dereplicate:
    input:
        cand=f"{TR}/selection/species/{{sid}}/candidate_genomes.tsv",
        samples=f"{RESULTS}/00_ingest/samples.parquet",
        paths=ANNOT_PATHS,
    output:
        derep=f"{TR}/work/{{sid}}/derep/dereplicated_genomes.tsv",
        clusters=f"{TR}/work/{{sid}}/derep/clusters.tsv",
    params:
        wd=f"{TR}/work/{{sid}}/derep",
        lookup=lambda wc, input: annot_lookup(f"{TR}/work/{wc.sid}/derep/ids.txt",
                                              "assembly", input.paths),
    threads: 4
    conda: "../envs/transition.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        # genome -> assembly path (from annotation_paths.tsv), then a symlink farm
        # of <genome>.fna so mash's sample names come back as genome ids
        tail -n +2 {input.cand} | cut -f1 > {params.wd}/ids.txt
        {params.lookup} > {params.wd}/assembly_paths.tsv
        {FARM} {params.wd}/assembly_paths.tsv {params.wd}/fa fna
        cp {params.wd}/fa/list.txt {params.wd}/fasta_list.txt
        mash sketch -o {params.wd}/ref -l {params.wd}/fasta_list.txt -k {TC[dereplicate][mash_kmer]} -s {TC[dereplicate][mash_sketch]} -p {threads}
        mash dist {params.wd}/ref.msh {params.wd}/ref.msh > {params.wd}/mash_dist.tsv
        python {SCRIPTS}/python/dereplicate_strains.py --config config/config.yaml \
          --candidates {input.cand} --mash-dist {params.wd}/mash_dist.tsv \
          --samples {input.samples} --out {output.derep} --clusters {output.clusters}
        """

rule panaroo_species:
    input:
        derep=rules.mash_dereplicate.output.derep,
        paths=ANNOT_PATHS,
    output:
        aln=f"{TR}/work/{{sid}}/panaroo/core_gene_alignment.aln",
        rtab=f"{TR}/work/{{sid}}/panaroo/gene_presence_absence.Rtab",
    params:
        wd=f"{TR}/work/{{sid}}/panaroo",
        lookup=lambda wc, input: annot_lookup(f"{TR}/work/{wc.sid}/panaroo/ids.txt",
                                              "prokka_gff", input.paths),
    threads: 8
    conda: "../envs/transition.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        # symlink farm of <genome>.gff so Panaroo's sample names, and hence the
        # tree's tip labels, are the genome ids
        tail -n +2 {input.derep} | cut -f1 > {params.wd}/ids.txt
        {params.lookup} > {params.wd}/gff_paths.tsv
        {FARM} {params.wd}/gff_paths.tsv {params.wd}/gff gff
        cp {params.wd}/gff/list.txt {params.wd}/gff_list.txt
        panaroo -i $(cat {params.wd}/gff_list.txt) -o {params.wd} \
          --clean-mode {TC[pangenome][clean_mode]} -a core \
          --aligner {TC[pangenome][aligner]} --core_threshold {TC[pangenome][core_threshold]} \
          -t {threads}
        """

rule gubbins_species:
    input:
        aln=rules.panaroo_species.output.aln,
    output:
        snps=f"{TR}/work/{{sid}}/gubbins/{{sid}}.filtered_polymorphic_sites.fasta",
    params:
        wd=f"{TR}/work/{{sid}}/gubbins", tool=TC["recombination"]["tool"]
    threads: 8
    conda: "../envs/transition.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}; cd {params.wd}
        if [ "{params.tool}" = "gubbins" ]; then
          run_gubbins.py --prefix {wildcards.sid} --threads {threads} \
            --iterations {TC[recombination][iterations]} {input.aln} \
          || snp-sites -o {wildcards.sid}.filtered_polymorphic_sites.fasta {input.aln}
        else
          snp-sites -o {wildcards.sid}.filtered_polymorphic_sites.fasta {input.aln}
        fi
        """

rule iqtree_species:
    input:
        snps=rules.gubbins_species.output.snps,
    output:
        tree=f"{TR}/work/{{sid}}/tree/{{sid}}.treefile",
    params:
        wd=f"{TR}/work/{{sid}}/tree",
        asc=("+ASC" if TC["tree"]["asc"] else "")
    threads: 8
    conda: "../envs/transition.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.wd}
        iqtree2 -s {input.snps} -m {TC[tree][model]}{params.asc} \
          -B {TC[tree][ufboot]} -alrt {TC[tree][alrt]} -T AUTO --threads-max {threads} \
          --prefix {params.wd}/{wildcards.sid} -redo \
        || iqtree2 -s {input.snps} -m GTR+G -B {TC[tree][ufboot]} -T AUTO \
          --prefix {params.wd}/{wildcards.sid} -redo
        """

rule ancestral_niche:
    input:
        tree=rules.iqtree_species.output.tree,
        derep=rules.mash_dereplicate.output.derep,
    output:
        summary=f"{TR}/work/{{sid}}/ancestral/ancestral_summary.tsv",
        nodes=f"{TR}/work/{{sid}}/ancestral/ancestral_nodes.rds",
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/ancestral_niche.R --config config/config.yaml --tree {{input.tree}} "
        f"--niche-map {{input.derep}} --species-id {{wildcards.sid}} "
        f"--out-dir {TR}/work/{{wildcards.sid}}/ancestral"

rule popgen_species:
    input:
        snps=rules.gubbins_species.output.snps,
        derep=rules.mash_dereplicate.output.derep,
        aln=rules.panaroo_species.output.aln,
    output:
        div=f"{TR}/work/{{sid}}/popgen/popgen_diversity.tsv",
    params: wd=f"{TR}/work/{{sid}}/popgen"
    conda: "../envs/python.yaml"
    shell:
        r"""
        L=$(awk '/^>/{{next}}{{print length($0); exit}}' {input.aln})
        python {SCRIPTS}/python/popgen_sfs.py --config config/config.yaml \
          --alignment {input.snps} --niche-map {input.derep} \
          --out-dir {params.wd} --core-length ${{L:-0}}
        """

rule demography_species:
    input:
        snps=rules.gubbins_species.output.snps,
        derep=rules.mash_dereplicate.output.derep,
    output:
        dirn=f"{TR}/work/{{sid}}/demography/directionality.tsv",
    params: wd=f"{TR}/work/{{sid}}/demography"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/demography_directionality.py --config config/config.yaml "
        f"--alignment {{input.snps}} --niche-map {{input.derep}} --out-dir {{params.wd}}"

rule accessory_species:
    input:
        rtab=rules.panaroo_species.output.rtab,
        derep=rules.mash_dereplicate.output.derep,
    output:
        summ=f"{TR}/work/{{sid}}/accessory/accessory_summary.tsv",
        full=f"{TR}/work/{{sid}}/accessory/accessory_differentiation.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/accessory_differentiation.py --config config/config.yaml "
        f"--rtab {{input.rtab}} --niche-map {{input.derep}} "
        f"--out {{output.full}} --summary {{output.summ}}"

rule transition_verdict:
    input:
        anc=rules.ancestral_niche.output.summary,
        pg=rules.popgen_species.output.div,
        dirn=rules.demography_species.output.dirn,
        acc=rules.accessory_species.output.summ,
    output:
        verdict=f"{TR}/work/{{sid}}/verdict/transition_verdict.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/transition_verdict.py --config config/config.yaml "
        f"--species-id {{wildcards.sid}} --ancestral {{input.anc}} --popgen {{input.pg}} "
        f"--divergence {{input.dirn}} --directionality {{input.dirn}} "
        f"--accessory {{input.acc}} --out {{output.verdict}}"

# ---- per-species figures ----------------------------------------------------
rule fig_transition_species:
    input: nodes=rules.ancestral_niche.output.nodes
    output: f"{TRFIG}/{{sid}}/tree.png"
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/fig_transition_species.R --config config/config.yaml "
        f"--nodes {{input.nodes}} --species-id {{wildcards.sid}} "
        f"--out {TRFIG}/{{wildcards.sid}}/tree"

rule fig_transition_popgen:
    input: div=rules.popgen_species.output.div
    output: f"{TRFIG}/{{sid}}/popgen.png"
    params: wd=f"{TR}/work/{{sid}}/popgen"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_transition_popgen.py --config config/config.yaml "
        f"--popgen-dir {{params.wd}} --species-id {{wildcards.sid}} "
        f"--out {TRFIG}/{{wildcards.sid}}/popgen"

# ---- gather: meta-analysis + report ----------------------------------------
def _all_verdicts(wildcards):
    return [f"{TR}/work/{s}/verdict/transition_verdict.tsv" for s in _species_ids(wildcards)]

def _all_species_figs(wildcards):
    out = []
    for s in _species_ids(wildcards):
        out += [f"{TRFIG}/{s}/tree.png", f"{TRFIG}/{s}/popgen.png"]
    return out

rule transition_meta:
    input: _all_verdicts
    output:
        dirn=f"{TR}/meta_directionality.tsv",
        summ=f"{TR}/meta_summary.json",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/transition_meta.py --config config/config.yaml "
        f"--verdict-glob '{TR}/work/*/verdict/transition_verdict.tsv' "
        f"--out-prefix {TR}/meta"

rule fig_transition_meta:
    input:
        dirn=rules.transition_meta.output.dirn,
        summ=rules.transition_meta.output.summ,
    output: f"{TRFIG}/meta/transition_meta.png"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/fig_transition_meta.py --config config/config.yaml "
        f"--directionality {{input.dirn}} --all-calls {TR}/meta_all_calls.tsv "
        f"--summary {{input.summ}} --out {TRFIG}/meta/transition_meta"

rule transition_report:
    input:
        meta=rules.transition_meta.output.summ,
        metafig=rules.fig_transition_meta.output,
        sppfigs=_all_species_figs,
    output:
        html=f"{RESULTS}/report/transition_report.html",
        md=f"{RESULTS}/report/transition_report.md",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/make_transition_report.py --config config/config.yaml "
        f"--results {RESULTS} --out-prefix {RESULTS}/report/transition_report"
