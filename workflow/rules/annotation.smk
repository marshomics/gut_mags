# =============================================================================
# annotation.smk -- de-novo annotation (dbCAN / antiSMASH / AMRFinder) and
# parsing of all five annotation sources into one uniform long table per source.
#
# Prokka, eggNOG and KOfam already exist (provided by the user); only dbCAN,
# antiSMASH and AMRFinder are generated here. All heavy steps are CHUNKED so the
# DAG holds ~hundreds of jobs, not 581k. Each chunk rule loops over the genomes
# in its chunk file; a sentinel marks completion and the parser reads the
# per-genome outputs back.
#
# No rule here templates a file path. Every input location is read from
# annotation_paths.tsv (rule annotation_paths in the Snakefile), which
# resolve_annotation_paths.py builds from either a manifest of explicit paths or
# the {genome} templates in the config. annot_lookup() projects one column of
# that table for the genomes in a chunk, as "<genome>\t<path>" lines.
# =============================================================================

ANN = config["inputs"]["annotations"]
REF = config["references"]

checkpoint split_genomes:
    input:
        samples=rules.ingest.output.parquet,
        reps=rules.representatives.output.txt,
    output:
        directory(f"{RESULTS}/02_annot/chunks"),
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/split_list.py --samples {{input.samples}} "
        f"--representatives {{input.reps}} --chunk-size {CHUNK} "
        f"--out-dir {{output}}"


def _chunk_ids(prefix):
    ck = checkpoints.split_genomes.get().output[0]
    return sorted(glob_wildcards(os.path.join(ck, prefix + "_{i}.txt")).i)


# ---- de-novo annotation (chunked shell loops) -------------------------------
rule run_dbcan:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        paths=ANNOT_PATHS,
    output:
        sentinel=f"{RESULTS}/02_annot/dbcan/chunk_{{i}}.done",
    params:
        db=REF["dbcan_db"],
        lookup=lambda wc, input: annot_lookup(input.chunk, ["prokka_faa", "dbcan_overview"],
                                              input.paths),
    threads: config["resources"]["default_threads"]
    conda: "../envs/annotation.yaml"
    shell:
        r"""
        set -euo pipefail
        {params.lookup} | while IFS=$'\t' read -r g faa ov; do
          outdir=$(dirname "$ov")
          mkdir -p "$outdir"
          if [ ! -s "$ov" ]; then
            run_dbcan "$faa" protein --db_dir {params.db} \
              --out_dir "$outdir" --dia_cpu {threads} --hmm_cpu {threads} || true
          fi
        done
        touch {output.sentinel}
        """

rule run_amrfinder:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        paths=ANNOT_PATHS,
    output:
        sentinel=f"{RESULTS}/02_annot/amrfinder/chunk_{{i}}.done",
    params:
        db=REF["amrfinder_db"],
        lookup=lambda wc, input: annot_lookup(input.chunk, ["prokka_faa", "amrfinder_tsv"],
                                              input.paths),
    threads: config["resources"]["default_threads"]
    conda: "../envs/annotation.yaml"
    shell:
        r"""
        set -euo pipefail
        {params.lookup} | while IFS=$'\t' read -r g faa tsv; do
          mkdir -p "$(dirname "$tsv")"
          if [ ! -s "$tsv" ]; then
            amrfinder -p "$faa" --database {params.db} --threads {threads} \
              -o "$tsv" || true
          fi
        done
        touch {output.sentinel}
        """

rule run_antismash:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/reps_{{i}}.txt",
        paths=ANNOT_PATHS,
    output:
        sentinel=f"{RESULTS}/02_annot/antismash/chunk_{{i}}.done",
    params:
        db=REF["antismash_db"],
        lookup=lambda wc, input: annot_lookup(input.chunk, ["assembly", "antismash_json"],
                                              input.paths),
    threads: config["resources"]["default_threads"]
    conda: "../envs/annotation.yaml"
    shell:
        r"""
        set -euo pipefail
        {params.lookup} | while IFS=$'\t' read -r g fna js; do
          outdir=$(dirname "$js")
          mkdir -p "$outdir"
          if [ ! -s "$js" ]; then
            antismash "$fna" --output-dir "$outdir" --databases {params.db} \
              --cpus {threads} --genefinding-tool prodigal \
              --output-basename "$g" || true
          fi
        done
        touch {output.sentinel}
        """


# ---- parse each source, per chunk -------------------------------------------
rule parse_kofam:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        paths=ANNOT_PATHS,
    output: f"{RESULTS}/02_annot/parsed/kofam/chunk_{{i}}.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/parse_annotations.py --config config/config.yaml --type kofam "
        f"--genome-list {{input.chunk}} --annotation-paths {{input.paths}} "
        f"--out {{output}}"

rule parse_eggnog:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        paths=ANNOT_PATHS,
    output: f"{RESULTS}/02_annot/parsed/eggnog/chunk_{{i}}.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/parse_annotations.py --config config/config.yaml --type eggnog "
        f"--genome-list {{input.chunk}} --annotation-paths {{input.paths}} "
        f"--out {{output}}"

rule parse_dbcan:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        sentinel=f"{RESULTS}/02_annot/dbcan/chunk_{{i}}.done",
        paths=ANNOT_PATHS,
    output: f"{RESULTS}/02_annot/parsed/dbcan/chunk_{{i}}.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/parse_annotations.py --config config/config.yaml --type dbcan "
        f"--genome-list {{input.chunk}} --annotation-paths {{input.paths}} "
        f"--out {{output}}"

rule parse_amrfinder:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/all_{{i}}.txt",
        sentinel=f"{RESULTS}/02_annot/amrfinder/chunk_{{i}}.done",
        paths=ANNOT_PATHS,
    output: f"{RESULTS}/02_annot/parsed/amrfinder/chunk_{{i}}.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/parse_annotations.py --config config/config.yaml --type amrfinder "
        f"--genome-list {{input.chunk}} --annotation-paths {{input.paths}} "
        f"--out {{output}}"

rule parse_antismash:
    input:
        chunk=f"{RESULTS}/02_annot/chunks/reps_{{i}}.txt",
        sentinel=f"{RESULTS}/02_annot/antismash/chunk_{{i}}.done",
        paths=ANNOT_PATHS,
    output: f"{RESULTS}/02_annot/parsed/antismash/chunk_{{i}}.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/parse_annotations.py --config config/config.yaml --type antismash "
        f"--genome-list {{input.chunk}} --annotation-paths {{input.paths}} "
        f"--out {{output}}"


# ---- aggregate parsed chunks into one long table per source -----------------
def _agg(prefix, parsed):
    def f(wildcards):
        ids = _chunk_ids(prefix)
        return [f"{RESULTS}/02_annot/parsed/{parsed}/chunk_{i}.parquet" for i in ids]
    return f

rule aggregate_kofam:
    input: _agg("all", "kofam")
    output: f"{RESULTS}/02_annot/ko.parquet"
    conda: "../envs/python.yaml"
    shell: f"{PY}/concat_parquet.py --inputs {{input}} --out {{output}}"

rule aggregate_eggnog:
    input: _agg("all", "eggnog")
    output: f"{RESULTS}/02_annot/eggnog.parquet"
    conda: "../envs/python.yaml"
    shell: f"{PY}/concat_parquet.py --inputs {{input}} --out {{output}}"

rule aggregate_dbcan:
    input: _agg("all", "dbcan")
    output: f"{RESULTS}/02_annot/dbcan.parquet"
    conda: "../envs/python.yaml"
    shell: f"{PY}/concat_parquet.py --inputs {{input}} --out {{output}}"

rule aggregate_amrfinder:
    input: _agg("all", "amrfinder")
    output: f"{RESULTS}/02_annot/amrfinder.parquet"
    conda: "../envs/python.yaml"
    shell: f"{PY}/concat_parquet.py --inputs {{input}} --out {{output}}"

rule aggregate_antismash:
    input: _agg("reps", "antismash")
    output: f"{RESULTS}/02_annot/antismash.parquet"
    conda: "../envs/python.yaml"
    shell: f"{PY}/concat_parquet.py --inputs {{input}} --out {{output}}"
