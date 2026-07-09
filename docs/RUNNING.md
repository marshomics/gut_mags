# Running the pipeline

Ordered from "nothing configured" to "all stages complete". Each stage names the
inputs it actually needs, so you can start the taxonomy stage today without
waiting for the annotation stack.

## 1. Inputs

### What you already have

| File | Rows | Used for | Config key |
|---|---|---|---|
| `total_metadata_qc_bbmap_animals_extra.txt` | 581,395 genomes | master sample sheet: drives every stage | `inputs.metadata` |
| `input/aa.presence.tsv` | 342,759 | amino-acid biosynthesis → auxotrophy layer | `community.inputs.auxotrophy` |
| `input/carbon.presence.tsv` | 342,759 | carbon-source utilisation layer | `community.inputs.carbon` |
| `input/predictions.tsv` | 342,759 | wide KEGG module presence → `modulep` layer | `community.inputs.modules` |
| `input/predictions_reliable.tsv` | 342,759 | predicted phenotypes → `trait` layer | `community.inputs.traits` |

The four `input/` tables key on the genome id in their first column
(`FREE000051`, `animal000005`, `GUT000001`, …), which matches `metadata.genome`
directly. `community_ingest.py` joins on that column whatever it is named. All
342,759 ids resolve.

Two files in `input/` are **not** read by the pipeline: `modules.long.tsv` (long
format, superseded by the wide `predictions.tsv`) and `genome_map.tsv` (G-id ↔
Prokka path). The `community.inputs.genome_map` config key is dead — no script
reads it. Leave it or delete it; it changes nothing.

### What you must point the config at

**Per-genome annotation.** Two ways to declare where the files are; pick either,
or mix them.

*Manifest.* Set `inputs.annotation_manifest` to a tab-separated table, one row
per genome, giving the full path to each file. Use this when the layout is
irregular — different directories per study, original basenames, symlink farms.
See `resources/annotation_manifest.example.tsv`.

```
genome        prokka_faa                              eggnog                                  kofam
GUT000001     /data/gut/prokka/GUT000001.faa          /data/gut/eggnog/GUT000001.emapper.a…   /data/gut/kofam/GUT000001.tsv
FREE000051    /scratch/env/prokka/PROKKA_0815.faa     /scratch/env/eggnog/SRR999.emapper.a…   /scratch/env/kofam/SRR999.tsv
```

Recognised columns, matched case-insensitively (aliases in brackets): `genome`
[`genome_id`, `id`, `file`], `assembly` [`fna`, `fasta`], `prokka_gff` [`gff`],
`prokka_faa` [`faa`, `proteins`], `prokka_ffn` [`ffn`, `cds`], `eggnog`
[`emapper`], `kofam`, and optionally `dbcan_overview`, `antismash_json`,
`amrfinder_tsv`. Unrecognised columns are ignored with a warning. A manifest may
be **partial**: any kind it omits falls back to the template, so you can list the
awkward files and template the rest.

*Templates.* Leave `annotation_manifest: null` and fill the `{genome}` paths in
`inputs.annotations`: `assemblies_fasta`, `prokka_gff`, `prokka_faa`,
`prokka_ffn`, `eggnog_annotations`, `kofam_annotations`.

Either way, `resolve_annotation_paths.py` writes
`results/00_ingest/annotation_paths.tsv` once and every rule reads only that. Two
consequences worth knowing:

- Every QC-passed genome must resolve to a path for each kind in
  `inputs.required_annotations` (default `prokka_faa`, `eggnog`, `kofam`). Add
  `assembly` and `prokka_gff` before running `transition_all` or the de-novo
  annotation rules; add `prokka_ffn` for the HyPhy selection analysis. A genome
  with no path is an error, not a skipped genome, because a silently dropped
  genome biases every prevalence estimate.
- `inputs.validate_paths` (`all` / `sample` / `none`, default `sample`) stats the
  declared files before the first job is submitted.

Mash and Panaroo name each sample after its file's basename, and the tree's tip
labels and the niche map key on the genome id. With a manifest those need not
coincide, so those rules first build a symlink farm of `<genome>.fna` /
`<genome>.gff` (`scripts/sh/symlink_farm.sh`). You do not have to do anything;
it just means arbitrary basenames are safe.

Reference data:

- `references.gtdb_data_dir`, `bac120_tree`, `ar53_tree`, `bac120_taxonomy`,
  `ar53_taxonomy` — GTDB release trees and taxonomy. **The taxonomy stage needs
  only these plus the metadata.**
- `references.gtdbtk_db` — only if you graft unplaced species.
- `references.dbcan_db`, `antismash_db`, `amrfinder_db` — only for the de-novo
  annotation rules.

Fetched once, over the internet:

```bash
python resources/fetch_kegg_modules.py  --out resources/kegg_modules.tsv
python resources/fetch_kegg_genesets.py --out-dir resources/genesets
```

The first is required by the KEGG module completeness rule; the second by the
enrichment stage's KO → pathway gene sets.

### What is missing

`western_nonwestern` is not a column in the metadata. Until you add it (values
matching `population.values`, human genomes only), `population.enabled` must stay
`false`, and `snakemake western_all` produces the redundancy report without the
Western vs non-Western contrast. Nothing else is blocked. When you add the
column, set `population.enabled: true` and the `western_vs_nonwestern` contrast is
injected automatically into the differential, Scoary2 and enrichment stacks — no
other edit.

## 2. Environment and scheduler

```bash
conda install -n base -c conda-forge -c bioconda snakemake'>=8' mamba
```

Per-rule conda environments live in `envs/` and are built on first use by
`--use-conda`. Nothing else needs installing by hand.

Three cluster profiles ship with the pipeline. Every command below is written
against SLURM; substitute the profile you use.

| Profile | Scheduler | Needs |
|---|---|---|
| `config/slurm` | SLURM | `snakemake-executor-plugin-slurm` |
| `config/sge` | SGE / UGE / OGS | `pip install snakemake-executor-plugin-sge` |
| `config/sge-generic` | SGE, via the generic executor | `pip install snakemake-executor-plugin-cluster-generic`; only `qsub`, `qstat`, `qacct`, `qdel` |

```bash
snakemake --workflow-profile config/sge --use-conda -j 200
```

Prefer `config/sge` for array-job submission; prefer `config/sge-generic` if you
would rather not depend on a third-party executor plugin. Validate either against
the live cluster before submitting anything:

```bash
python scripts/python/check_sge_profile.py config/sge
```

That checks every rule's `runtime` against the `h_rt` ceiling of the queue it is
routed to, that the parallel environment and queues exist, and that the memory
complex being requested is actually *consumable*.

It works from any submit host, which takes some doing. `qconf -sc`, `qconf -spl`
and `qconf -sq` are permitted only from an **admin host** (`qconf -sh` lists
them; `node519` is one, `morty` is not). Elsewhere they answer

```
denied: host "morty.eb.local" is not an admin host
```

which is a refusal to answer, not the answer "no". So the checker draws on three
sources and says which one each fact came from:

- **qconf**, when it answers. Its values are cross-checked against
  `config/sge/cluster.yaml` and any disagreement is reported as `DRIFT`.
- **`config/sge/cluster.yaml`**, the facts recorded earlier on an admin host,
  with provenance. Regenerate with `bash config/sge/record_cluster_facts.sh >
  config/sge/cluster.yaml` when the cluster changes.
- **`qsub -w v`**, a dry scheduling run against an empty cluster. It submits
  nothing, works anywhere, and rejects a request whose `h_rt` exceeds the queue
  ceiling or whose queue, PE or complex is unknown. It cannot see whether a
  complex is *consumable*, which is the one thing that still needs qconf.

Anything none of the three can settle is printed as `unverifiable from this
host`, with the command to run on an admin host. It does not fail the check —
an unanswered question is not a failed answer.

Both profiles are tuned to this cluster:

| | |
|---|---|
| `standard.q` | `h_rt` 24 h, 865 slots — everything under ~20 h |
| `long.q` | `h_rt` 672 h, 2112 slots — `scaffold`, `iqtree_species`, `clade_panaroo`, `selection_hyphy` |
| parallel env | `parallel` (there is no `smp`), used as `-pe parallel <threads>` |
| memory | `h_vmem`, requestable **and consumable** |

Three things about memory on this cluster, all of which bite silently.

`mem_free` is requestable but **not consumable**, so requesting it reserves
nothing: it is checked once at dispatch and never tracked, and jobs oversubscribe
the node until the OOM killer picks a winner. Only `h_vmem` reserves. Set
`SGE_MEM_FREE_TOO=1` to co-request `mem_free` as a dispatch filter if you want
both.

`h_vmem` is consumed **per slot**. Under `-pe parallel 8`, `h_vmem=40G` reserves
320 GB. `mem_mb` in the profiles is the job total and is divided by the thread
count before submission.

`h_vmem` caps **address space, not resident memory**. A tool that maps far more
than it touches — diamond, some JVM and glibc-arena patterns — is killed with a
small RSS. If a job dies on memory while `qacct` shows modest `maxvmem`, raise
the request before doubting the tool.

Threads are stated explicitly in every profile under `set-threads`, because they
decide two things at once: the `-pe parallel N` slot count, and the divisor that
turns a job's total `mem_mb` into its per-slot `h_vmem`. Four rules (`scaffold`,
`phyloglm_chunk`, `species_profiles`, `ordination`) declare no `threads:` and so
get 1 — their parallelism comes from chunking, not cores. Reserving 8 CPUs for
them would idle 7.

The heaviest single job is `scaffold` (pruning the GTDB bac120 tree): single
threaded, so its 320 GB lands on **one slot** as `h_vmem=320000M`. Every node has
at least 504 GB, so it fits, but it only starts once a node has 320 GB free —
which is why it is routed to `long.q` despite needing only 24 h.

## 3. Configure

Fill the `CHANGE_ME` paths:

```bash
grep -n CHANGE_ME config/config.yaml
```

Metadata, GTDB reference data, GTDB-Tk DB, the community `input/` tables, and the
annotation databases. Per-genome annotations are either
`inputs.annotation_manifest` or the `{genome}` templates in `inputs.annotations`,
as above.

Then set your scheduler settings: partition and account in
`config/slurm/config.yaml`, or the queue names in `config/sge/config.yaml` /
`config/sge-generic/config.yaml`.

Everything scientific is in `config/config.yaml` with the reasoning beside it.
Nothing scientific is hard-coded in a script.

## 4. Verify before you burn cluster time

```bash
bash tests/run_tests.sh          # compiles, validates config, runs the statistics
                                 # and end-to-end synthetic-data tests
snakemake -n -p                  # dry-run the DAG for the default target
snakemake -n -p functional_all   # dry-run the heavy stage
```

`run_tests.sh` takes a few minutes and needs no cluster. It will catch a
malformed config, a rule pointing at a script that does not exist, and an R
package attached but absent from `envs/*.yaml`.

## 5. Stage 1 — taxonomy by niche (default target)

Needs: metadata + GTDB trees. Nothing else. Runs today.

```bash
snakemake --workflow-profile config/slurm --use-conda -j 500
open results/report/taxonomy_report.html
```

Produces per-rank niche specificity with a permutation null, indicator taxa
(IndVal.g), taxon enrichment, the curveball cross-niche overlap null, β-diversity
turnover/nestedness, phylogenetic community structure (PD/NRI/NTI), rarefied
novelty, and the host-resolved (mouse) analysis.

## 6. Stage 2 — annotation and species profiles

The expensive step, and where the compute decision lives. `functional_annotation.scope`
currently says `dbcan: all_genomes` and `amrfinder: all_genomes` — that is
~342k–581k dbCAN and AMRFinder jobs. CAZyme prevalence genuinely wants all
genomes, but if you need to bound the cost, set `dbcan: representatives` and note
the change in the methods. antiSMASH is already restricted to representatives.

```bash
snakemake results/03_profiles/prevalence_ko.parquet \
  --workflow-profile config/slurm --use-conda -j 500
```

Or just run stage 3, which pulls this in.

## 7. Stage 3 — gene/function signatures

```bash
snakemake functional_all --workflow-profile config/slurm --use-conda -j 500
open results/report/report.html
```

This runs the four-method consensus (phyloglm, clade-stratified CMH, balanced
bootstrap, Scoary2) across every layer and contrast, plus the permutation null,
ordination/PERMANOVA/varpart, and the comparative stack.

Scoary2 and enrichment are wired into `functional_all`, but each also has a
standalone target if you want to iterate on one:

```bash
snakemake scoary_all     --workflow-profile config/slurm --use-conda -j 500
snakemake enrichment_all --workflow-profile config/slurm --use-conda -j 500
```

Read `results/05_diff/<layer>/<contrast>_signatures.json` before the figures. It
reports `n_tier1_all_methods`, `n_tier2_partial` and `n_untestable_per_method`.
If CMH is untestable for most features, niche is nested within genus at the
stratifying rank and you should say so in the manuscript, or coarsen
`stats.cmh.stratify_rank`.

## 8. Stage 4 — within-species niche transitions

Independent of stage 3. Needs assemblies and Prokka output. Per-species Panaroo →
Gubbins → IQ-TREE, so cost scales with the number of qualifying species (those
with ≥ 10 near-complete strains in each of ≥ 2 niches).

```bash
snakemake transition_all --workflow-profile config/slurm --use-conda -j 500
open results/report/transition_report.html
```

Check `results/*/directionality.tsv` for `call_status`. `unresolved` is a real
answer, not a failure: it means the two populations do not differ in nestedness,
private alleles or π once sample size is equalised.

## 9. Stage 5 — communities, redundancy, synthesis

```bash
snakemake community_all   --workflow-profile config/slurm --use-conda -j 200
snakemake redundancy_all  --workflow-profile config/slurm --use-conda -j 200
snakemake western_all     --workflow-profile config/slurm --use-conda -j 200  # needs the column
snakemake synthesis_all   --workflow-profile config/slurm --use-conda -j 200
open results/report/synthesis_report.html
```

`community_all` needs only the metadata and the four `input/` tables, so it can
run alongside stage 1. `synthesis_all` depends on stages 1, 3 and the enrichment
stack, and builds the three human-specific catalogues and the master figure.

## 10. Sensitivity

After the main run, not before:

```bash
bash scripts/sh/run_sensitivity.sh "--workflow-profile config/slurm --use-conda -j 300"
```

Sweeps the presence threshold (0.1 / 0.5 / 0.9), HQ-only genomes, drop-mouse, and
phylogeny-off, writing each to its own `results_*` directory for comparison.

## Dependency summary

```
metadata + GTDB ──> all (taxonomy)          ──┐
                                              ├──> synthesis_all
metadata + input/ ─> community_all            │
                                              │
assemblies+prokka ─> annotation ─> functional_all ─> enrichment_all ─┘
                  └> transition_all                └> scoary_all

metadata + profiles ─> redundancy_all ─> western_all (needs western_nonwestern)
```

## Order of operations, if you want the short version

1. `bash tests/run_tests.sh`, then (SGE) `python scripts/python/check_sge_profile.py config/sge`
2. fill the `CHANGE_ME` paths, and the queue (SGE) or partition/account (SLURM)
3. `python resources/fetch_kegg_modules.py --out resources/kegg_modules.tsv`
4. `snakemake -n -p`
5. `snakemake ... -j 500` (taxonomy) and `snakemake community_all ...` in parallel
6. `snakemake functional_all ...` once annotation paths resolve
7. `snakemake transition_all ...`
8. add `western_nonwestern`, flip `population.enabled`, `snakemake western_all ...`
9. `snakemake synthesis_all ...`
10. `bash scripts/sh/run_sensitivity.sh ...`
