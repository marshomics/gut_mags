# Human gut niche signature pipeline

A Snakemake workflow to identify, defensibly, the genes and functions that
distinguish human-gut microbial genomes from animal-gut and free-living genomes.
The guiding question: what makes the human gut *human* at the gene/function
level, once you remove every reason a difference could be an artefact rather
than biology.

The dataset is 581,395 metagenome-assembled and isolate genomes spanning three
niches (human gut, animal gut, free-living), with GTDB taxonomy, CheckM-style
quality, and host taxonomy for animal genomes. The same metadata table drives
the whole workflow as the master sample sheet.

## The one design decision that everything follows from

Closely related genomes share gene content because they share an ancestor, not
because the human gut selected for it. The most-sequenced species are human-gut
commensals with thousands of near-identical genomes each; the animal niche is
~80% mouse; free-living holds 46k species against human's 6k. A naive
"compare gene frequencies between niches" would mostly measure these sampling
and ancestry artefacts.

So the pipeline does three things consistently, everywhere:

1. **The species is the unit of analysis.** Genomes are collapsed to species-
   level functional profiles before any comparison, so a 9,606-genome species
   and a singleton each count once.
2. **Phylogeny is controlled, not ignored.** Every differential call must pass a
   phylogenetic-covariance model *and* a clade-matched stratified test *and* a
   balance-corrected bootstrap before it is reported.
3. **Quality and genome size are covariates in every model**, never silent.

A finding is only called a "human signature" when all three confounder controls
agree in direction and significance. The full mapping of confounder to control
is in [`docs/DESIGN_RATIONALE.md`](docs/DESIGN_RATIONALE.md); read that first if
you are reviewing the science.

## Staging: taxonomy first, function later

The current stage is the species/taxonomy-by-niche analysis. `snakemake` (default
target) builds the **taxonomy report** and only its dependencies (the metadata
and the GTDB scaffold); it does not need the functional annotations, so it runs
now. The gene/function stack is built later with `snakemake functional_all`. The
taxonomy stage answers, with the same confounder discipline applied to taxonomy:
which taxa characterise each niche, how phylogenetically restricted each niche
is, how much overlap and novelty there is, and how much of the animal signal is
mouse.

A third, separate stage (`snakemake transition_all`) handles the species that
span more than one niche: for each one with enough near-complete strains per
niche, it tests whether a niche is a recent acquisition. Per species it
dereplicates strains within each niche, builds a Panaroo core genome, masks
recombination with Gubbins, infers an outgroup-rooted IQ-TREE, polarizes the
niche history by ancestral-state reconstruction, and looks for the founder
signature (reduced diversity, negative Tajima's D, nestedness, gene gains) in the
candidate derived niche, then asks across all such species whether acquisitions
go in a consistent direction. See `docs/METHODS.md` (Within-species niche
transitions) and `resources/fastsimcoal2_recipe.md` for the optional coalescent
split/migration model.

## What it produces

- Species-level functional profiles for KO, KEGG modules, Pfam, COG, CAZymes,
  BGCs and AMR genes.
- Consensus human-signature features per functional layer, for two contrasts
  (human vs everything, human vs animal gut), each with effect size, CI, FDR,
  per-method support and an empirical permutation FDR.
- Phylogenetic-signal, PGLS and ancestral-state/convergence results that say
  whether each signature was inherited once or acquired repeatedly.
- Ordination, PERMANOVA and variation partitioning that quantify how much
  functional variation niche explains after phylogeny is removed.
- Publication figures, each as PNG (400 dpi) and SVG with editable text.
- An HTML/markdown report and a sensitivity comparison across four alternative
  analytical choices.

## Prerequisites

- Snakemake >= 8 and conda/mamba. Per-rule environments are in `envs/`.
- A cluster profile: SLURM, or SGE/UGE via `config/sge` or `config/sge-generic`.
- The genomes' existing annotations: Prokka, eggNOG-mapper, KofamScan outputs.
- Assemblies (for de novo dbCAN / antiSMASH / AMRFinder, which the workflow runs).
- GTDB release reference trees + taxonomy (bac120, ar53), and GTDB-Tk DB if you
  graft unplaced species.

## Install

The target cluster has **no internet**, and no package may come from the
`defaults` or `anaconda` channels. Both constraints are handled; the full
procedure is [`docs/INSTALL.md`](docs/INSTALL.md).

Every file in `envs/` declares `conda-forge`, `bioconda`, `nodefaults` in that
order; `config/condarc` pins strict channel priority and empties
`default_channels`; and the build script passes `--override-channels` so a stray
`~/.condarc` cannot reintroduce them. `check_env_channels.py` refuses to build if
any of that is violated, and `verify_envs.sh` afterwards reads `conda list
--explicit` to check the URL each package *actually* came from.

Offline deployment works because Snakemake stores each environment in
`<conda-prefix>/<hash>`, where the hash is the MD5 of the environment file's
content — not of the machine or the path. So environments built elsewhere are
already named what the cluster looks for. Only the absolute paths inside them
need care: build at a path that exists on both machines, or use `--pack`
(conda-pack) and let `conda-unpack` rewrite them.

On a Linux x86_64 machine with internet:

```bash
mamba create -p ./bootstrap -c conda-forge --override-channels snakemake=8.20 mamba
conda activate ./bootstrap
export CONDARC=$PWD/config/condarc

bash scripts/sh/build_envs.sh /shared/hgn/conda-envs          # add --pack if paths differ
rsync -a /shared/hgn/conda-envs/ cluster:/shared/hgn/conda-envs/
rsync -a deploy/ cluster:/path/to/repo/deploy/
```

Then on the cluster:

```bash
bash scripts/sh/deploy_envs.sh /shared/hgn/conda-envs
```

which verifies the environments and finishes with `snakemake --sdm conda
--conda-create-envs-only` on the air-gapped host — so a missing or
hash-mismatched environment fails in seconds rather than three hours into the
first annotation job.

The build machine must be Linux x86_64: `qhost` reports `lx-amd64` for every
node, and conda environments do not cross platforms.

## Configure

Everything scientific lives in [`config/config.yaml`](config/config.yaml), with
the rationale beside each parameter. Set the `CHANGE_ME` paths (metadata, GTDB
data, reference DBs) and declare where the per-genome annotations live — either
`inputs.annotation_manifest` (a TSV of explicit paths) or the `{genome}`
templates in `inputs.annotations`. Nothing is hard-coded in the scripts.

## Run

Step-by-step, with the input inventory and what each stage actually needs:
[`docs/RUNNING.md`](docs/RUNNING.md). Verify first with `bash tests/run_tests.sh`;
the correctness audit is in [`docs/AUDIT.md`](docs/AUDIT.md).

```bash
# environments installed and KEGG tables fetched already: see docs/INSTALL.md
export CONDA_ENVS=/shared/hgn/conda-envs

# check the statistics and the wiring before spending cluster time
bash tests/run_tests.sh

# sanity-check the DAG without running anything
snakemake -n -p

# CURRENT STAGE: taxonomy by niche (default target; no functional data needed)
snakemake --workflow-profile config/slurm --sdm conda --conda-prefix $CONDA_ENVS -j 500
open results/report/taxonomy_report.html

# WITHIN-SPECIES niche transitions (separate, per-species stage)
snakemake transition_all --workflow-profile config/slurm --sdm conda --conda-prefix $CONDA_ENVS -j 500
open results/report/transition_report.html

# LATER: gene/function stage
snakemake functional_all --workflow-profile config/slurm --sdm conda --conda-prefix $CONDA_ENVS -j 500
open results/report/report.html

# stability checks (after the main run)
bash scripts/sh/run_sensitivity.sh "--workflow-profile config/slurm --sdm conda --conda-prefix $CONDA_ENVS -j 300"
```

The author runs this; the pipeline is not executed during its construction.
Cluster profiles for SLURM (`config/slurm`) and SGE/UGE (`config/sge`, or
`config/sge-generic` which needs only qsub/qstat/qacct) are in `config/`; tune
per-rule resources there. `bash tests/run_tests.sh` checks that every per-rule
key in a profile names a real rule, because Snakemake ignores a typo silently;
`python scripts/python/check_sge_profile.py config/sge` checks the queue routing,
the parallel environment and the memory complex against the live cluster.

## Layout

```
config/        config.yaml (all parameters) + condarc + slurm/ sge/ sge-generic/ profiles
envs/          per-rule conda environments (conda-forge + bioconda + nodefaults)
workflow/      Snakefile + rules/{annotation,profiles,comparative,figures,report}.smk
scripts/python annotation parsing, profiles, differential methods, figures, report
scripts/R      phyloglm, signal, PGLS, ancestral, ordination/varpart, tree figures
scripts/sh     env build/deploy/verify, sensitivity sweep, symlink farm
resources/     KEGG module fetcher + reference-data notes
docs/          INSTALL, RUNNING, DESIGN_RATIONALE, METHODS, OUTPUTS, AUDIT
```

## Scope and honesty about limits

Sample source (study, extraction, assembler) is partly confounded with niche and
cannot be fully removed without per-sample provenance; the species-level,
phylogenetically controlled design mitigates it and the limitation is stated in
`docs/METHODS.md`. Phylogenetic tests run on the GTDB-placeable species
(coverage is reported per niche). GTDB placeholder/de-novo species clusters are
valid units but are not cross-database named species, so cross-niche overlap is
reported for named and placeholder clusters separately.
