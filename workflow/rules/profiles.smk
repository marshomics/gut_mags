# =============================================================================
# profiles.smk -- species-level prevalence profiles, KEGG module completeness,
# and the continuous species trait table. This is where genome-level annotations
# become species-level units (strains-per-species confounder removed).
# =============================================================================

rule species_profiles:
    input:
        samples=rules.ingest.output.parquet,
        kofam=f"{RESULTS}/02_annot/ko.parquet",
        eggnog=f"{RESULTS}/02_annot/eggnog.parquet",
        dbcan=f"{RESULTS}/02_annot/dbcan.parquet",
        antismash=f"{RESULTS}/02_annot/antismash.parquet",
        amrfinder=f"{RESULTS}/02_annot/amrfinder.parquet",
    output:
        ko=f"{RESULTS}/03_profiles/prevalence_ko.parquet",
        pfam=f"{RESULTS}/03_profiles/prevalence_pfam.parquet",
        cog=f"{RESULTS}/03_profiles/prevalence_cog.parquet",
        ec=f"{RESULTS}/03_profiles/prevalence_ec.parquet",
        cazyme=f"{RESULTS}/03_profiles/prevalence_cazyme.parquet",
        bgc=f"{RESULTS}/03_profiles/prevalence_bgc.parquet",
        amr=f"{RESULTS}/03_profiles/prevalence_amr.parquet",
        ko_conc=f"{RESULTS}/03_profiles/ko_source_concordance.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/species_prevalence_profiles.py --config config/config.yaml "
        f"--samples {{input.samples}} --kofam {{input.kofam}} "
        f"--eggnog {{input.eggnog}} --dbcan {{input.dbcan}} "
        f"--antismash {{input.antismash}} --amrfinder {{input.amrfinder}} "
        f"--out-dir {RESULTS}/03_profiles"

rule module_completeness:
    input:
        ko=rules.species_profiles.output.ko,
        moddef=config["references"]["kegg_module_def"],
    output:
        f"{RESULTS}/03_profiles/prevalence_module.parquet",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/kegg_module_completeness.py --config config/config.yaml "
        f"--ko-prevalence {{input.ko}} --module-def {{input.moddef}} "
        f"--out {{output}}"

rule species_traits:
    input:
        samples=rules.ingest.output.parquet,
        species=rules.species_units.output.table,
        profiles=rules.species_profiles.output.ko,   # dependency anchor
    output:
        f"{RESULTS}/03_profiles/species_traits.tsv",
    conda: "../envs/python.yaml"
    shell:
        f"{PY}/species_trait_table.py --config config/config.yaml "
        f"--samples {{input.samples}} --species-table {{input.species}} "
        f"--profiles-dir {RESULTS}/03_profiles --out {{output}}"
