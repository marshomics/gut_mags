# =============================================================================
# comparative.smk -- the analytical core.
#   scaffold -> differential prep -> {phyloglm (chunked), CMH, resampling}
#   -> consensus -> permutation null; then combined signatures feed phylogenetic
#   signal, PGLS, ancestral/convergence; plus ordination + variation partition.
# Every per-(layer,contrast) branch uses the SAME prepared inputs.
# =============================================================================

FCHUNK = config["resources"]["feature_chunk_size"]

# ---- GTDB scaffold ----------------------------------------------------------
rule scaffold:
    input:
        species=rules.species_units.output.table,
    output:
        tree=f"{RESULTS}/04_phylo/species_tree.nwk",
        tipmap=f"{RESULTS}/04_phylo/tip_map.tsv",
        cov=f"{RESULTS}/04_phylo/scaffold_coverage.json",
    conda: "../envs/phylo.yaml"
    shell:
        f"{PY}/build_species_scaffold.py --config config/config.yaml "
        f"--species-table {{input.species}} --out-dir {RESULTS}/04_phylo"

# ---- differential prep (per layer x contrast) -------------------------------
rule differential_prep:
    input:
        samples=rules.ingest.output.parquet,
        species=rules.species_units.output.table,
        prevalence=f"{RESULTS}/03_profiles/prevalence_{{layer}}.parquet",
    output:
        analysis=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_analysis_species.tsv",
        presence=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_presence.parquet",
        feats=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_tested_features.txt",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/differential_prep.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} "
        f"--prevalence {{input.prevalence}} --layer {{wildcards.layer}} "
        f"--contrast {{wildcards.contrast}} "
        f"--out-prefix {RESULTS}/05_diff/{{wildcards.layer}}/{{wildcards.contrast}}"

# ---- METHOD A: phyloglm (feature-chunked) -----------------------------------
checkpoint split_features:
    input:
        feats=rules.differential_prep.output.feats,
    output:
        directory(f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_fchunks"),
    run:
        os.makedirs(output[0], exist_ok=True)
        feats = [l.strip() for l in open(input.feats) if l.strip()]
        for k in range(0, max(len(feats), 1), FCHUNK):
            with open(os.path.join(output[0], f"fc_{k // FCHUNK:04d}.txt"), "w") as fh:
                fh.write("\n".join(feats[k:k + FCHUNK]) + "\n")

rule phyloglm_chunk:
    input:
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
        analysis=rules.differential_prep.output.analysis,
        presence=rules.differential_prep.output.presence,
        fchunk=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_fchunks/fc_{{k}}.txt",
    output:
        f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_phyloglm/fc_{{k}}.tsv",
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/phyloglm_enrichment.R --tree {{input.tree}} --tip-map {{input.tipmap}} "
        f"--analysis {{input.analysis}} --presence {{input.presence}} "
        f"--features {{input.fchunk}} "
        f"--btol {config['stats']['phyloglm']['btol']} "
        f"--boot {config['stats']['phyloglm']['boot']} --out {{output}}"

def _phyloglm_shards(wildcards):
    ck = checkpoints.split_features.get(layer=wildcards.layer,
                                        contrast=wildcards.contrast).output[0]
    ks = sorted(glob_wildcards(os.path.join(ck, "fc_{k}.txt")).k)
    return [f"{RESULTS}/05_diff/{wildcards.layer}/{wildcards.contrast}_phyloglm/fc_{k}.tsv"
            for k in ks]

rule phyloglm_combine:
    input: _phyloglm_shards
    output: f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_phyloglm.tsv"
    conda: "../envs/python.yaml"
    shell: f"{PY}/combine_tsv.py --inputs {{input}} --out {{output}}"

# ---- METHOD B: clade-stratified CMH -----------------------------------------
rule cmh:
    input:
        analysis=rules.differential_prep.output.analysis,
        presence=rules.differential_prep.output.presence,
    output: f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_cmh.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/cmh_stratified.py --config config/config.yaml "
        f"--analysis {{input.analysis}} --presence {{input.presence}} --out {{output}}"

# ---- METHOD C: balanced resampling ------------------------------------------
rule resampling:
    input:
        analysis=rules.differential_prep.output.analysis,
        presence=rules.differential_prep.output.presence,
        species=rules.species_units.output.table,
    output: f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_resampling.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/balanced_resampling.py --config config/config.yaml "
        f"--analysis {{input.analysis}} --presence {{input.presence}} "
        f"--species-table {{input.species}} --contrast {{wildcards.contrast}} "
        f"--out {{output}}"

# ---- consensus + permutation null -------------------------------------------
# Scoary2 enters as a 4th method for the contrasts in scoary.contrasts (its parse
# output is produced in scoary.smk). consensus_signatures.py treats a method as
# required only if its input is present, so other contrasts use the trio alone.
def _scoary_input(wildcards):
    sc = config.get("scoary", {})
    if sc.get("enabled") and wildcards.contrast in sc.get("contrasts", []) \
            and wildcards.layer in sc.get("layers", []):
        return f"{RESULTS}/scoary/species/{wildcards.layer}/{wildcards.contrast}/scoary.tsv"
    return []

rule consensus:
    input:
        pg=rules.phyloglm_combine.output,
        cmh=rules.cmh.output,
        rs=rules.resampling.output,
        scoary=_scoary_input,
    output:
        tsv=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.tsv",
        json=f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_signatures.json",
    params:
        scoary_arg=lambda wc, input: f"--scoary {input.scoary}" if input.scoary else "",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/consensus_signatures.py --config config/config.yaml "
        f"--phyloglm {{input.pg}} --cmh {{input.cmh}} --resampling {{input.rs}} {{params.scoary_arg}} "
        f"--layer {{wildcards.layer}} --contrast {{wildcards.contrast}} "
        f"--out-prefix {RESULTS}/05_diff/{{wildcards.layer}}/{{wildcards.contrast}}_signatures"

rule permutation_null:
    input:
        analysis=rules.differential_prep.output.analysis,
        presence=rules.differential_prep.output.presence,
    output: f"{RESULTS}/05_diff/{{layer}}/{{contrast}}_permnull.json"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/permutation_null.py --config config/config.yaml "
        f"--analysis {{input.analysis}} --presence {{input.presence}} "
        f"--contrast {{wildcards.contrast}} --out {{output}}"

# ---- combine signatures + presence across layers (per contrast) -------------
rule combine_signatures:
    input:
        sig=expand(f"{RESULTS}/05_diff/{{layer}}/{{{{contrast}}}}_signatures.tsv",
                   layer=DIFF_LAYERS),
    output: f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/combine_signatures.py --inputs {{input.sig}} "
        f"--layers {' '.join(DIFF_LAYERS)} --out {{output}}"

rule combine_presence:
    input:
        prev=expand(f"{RESULTS}/03_profiles/prevalence_{{layer}}.parquet",
                    layer=DIFF_LAYERS),
    output: f"{RESULTS}/05_diff/combined/presence_all.parquet"
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/concat_parquet.py --inputs {{input.prev}} --out {{output}}"

# ---- phylogenetic comparative on the combined signatures --------------------
rule phylo_signal:
    input:
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
        traits=rules.species_traits.output,
        sig=f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv",
        presence=rules.combine_presence.output,
    output: directory(f"{RESULTS}/06_comparative/{{contrast}}_signal")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/phylo_signal.R --config config/config.yaml --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --traits {{input.traits}} "
        f"--signatures {{input.sig}} --presence {{input.presence}} "
        f"--out-dir {{output}}"

rule pgls:
    input:
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
        traits=rules.species_traits.output,
    output: f"{RESULTS}/06_comparative/pgls_results.tsv"
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/pgls_continuous.R --config config/config.yaml --tree {{input.tree}} "
        f"--tip-map {{input.tipmap}} --traits {{input.traits}} --out {{output}}"

rule ancestral:
    input:
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
        sig=f"{RESULTS}/05_diff/combined/{{contrast}}_signatures_all.tsv",
        presence=rules.combine_presence.output,
    output: directory(f"{RESULTS}/06_comparative/{{contrast}}_ancestral")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/ancestral_convergence.R --config config/config.yaml "
        f"--tree {{input.tree}} --tip-map {{input.tipmap}} "
        f"--signatures {{input.sig}} --presence {{input.presence}} "
        f"--out-dir {{output}}"

# ---- ordination / variation partitioning ------------------------------------
rule ordination:
    input:
        prevalence=f"{RESULTS}/03_profiles/prevalence_{ORD_LAYER}.parquet",
        traits=rules.species_traits.output,
        tree=rules.scaffold.output.tree,
        tipmap=rules.scaffold.output.tipmap,
    output: directory(f"{RESULTS}/07_ordination")
    conda: "../envs/r.yaml"
    shell:
        f"{RS}/ordination_varpart.R --config config/config.yaml "
        f"--prevalence {{input.prevalence}} --species-traits {{input.traits}} "
        f"--tree {{input.tree}} --tip-map {{input.tipmap}} --out-dir {{output}}"
