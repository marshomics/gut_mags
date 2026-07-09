# Design rationale: confounders and controls

This document is the defensibility argument. For each threat to validity it
states the problem, the control, where the control lives in the code, and what a
reviewer should check. The data motivating each point come from the provided
metadata (581,395 genomes; human 342,875 / animal 142,183 / free 96,337;
59,622 species).

## The analysis unit: species, not genome

Genome counts per species range from 1 to 9,606 (median 1; 65% singletons). The
top species by genome count are all human-gut commensals (*Bacteroides
uniformis*, *Alistipes putredinis*, *Escherichia coli*, *Phocaeicola vulgatus*,
*Agathobacter rectalis*). Counting genomes would let a handful of well-sequenced
human commensals dominate every comparison.

Control: collapse genomes to species before any comparison. A species' functional
profile is the prevalence of each feature across its conspecific genomes;
"present" means prevalence above a threshold (default 0.5, swept 0.1/0.5/0.9).
Code: `species_prevalence_profiles.py`, `define_species_units.py`. Prevalence
across many genomes also buffers single-genome incompleteness.

Check: `results/00_ingest/species_units_report.json` (singleton count, max genomes
per species) and the prevalence threshold sweep in the sensitivity comparison.

## Confounder 1 — strains per species

Covered by the species-unit decision above. Additionally, the balanced bootstrap
draws one genome-equivalent per species per iteration, so within-species
oversampling cannot leak into the resampling CIs.
Code: `balanced_resampling.py`.

## Confounder 2 — species per niche (human 6k vs free 46k)

Unequal group sizes distort any test that compares counts, and inflate
PERMANOVA through unequal dispersion.

Controls:
- Prevalences are proportions within a niche, not counts.
- Diversity is compared at matched sampling effort: sample-size rarefaction and
  coverage-based Hill numbers (q=0,1,2). Code: `rarefaction_diversity.py`.
- The balanced bootstrap draws an equal number of species per niche every
  iteration (N = the smaller group). Code: `balanced_resampling.py`.
- The differential prevalence filter is applied per niche (max across niches),
  so removing rare features does not favour the larger niche. Code:
  `differential_prep.py`.
- PERMANOVA is paired with `betadisper` so a significant niche term is not a
  dispersion artefact. Code: `ordination_varpart.R`.

## Confounder 3 — host imbalance within the animal niche (~80% mouse)

The animal niche is 113k mouse genomes out of 142k. "Animal gut" could quietly
mean "mouse gut".

Controls:
- Animal species enter the balanced bootstrap with weights inversely
  proportional to their dominant host's frequency, so mouse-derived species
  cannot dominate the animal contribution (`cap_dominant_host_fraction`).
  Code: `balanced_resampling.py`.
- A sensitivity run drops mouse entirely (`config/sensitivity/nomouse.yaml`) and
  the signature overlap is reported. Code: `make_sensitivity_configs.py`,
  `compare_sensitivity.py`.
- Host breadth per species is recorded for interpretation. Code:
  `define_species_units.py`.

## A fourth method: Scoary2 pan-GWAS

Scoary2 is added as an independent fourth association method. Its pairwise-
comparisons (Maddison) algorithm corrects for clonal and phylogenetic structure
by a different route than phyloglm (continuous covariance) or CMH (discrete
clade matching), so agreement across all four is strong evidence. It is run at
the species level on the functional presence matrices (species as the unit, so
strain bias is already removed) with the GTDB species tree, and a hit must pass
the corrected Fisher q, the label-switching permutation p (which flags signals
that merely track lineages), the count of supporting pairwise comparisons and
the best-pairwise p. The consensus treats a method as required only if it ran for
that contrast, so Scoary2 strengthens the calls where it applies without forcing
all contrasts through the heavier trio. A clade-stratified run (Panaroo per
genus) gives ortholog-family resolution where orthology is meaningful. The input
is deliberately the functional matrices, not a global 581k-genome de-novo
pangenome, because ortholog families do not span ~60k species and such a matrix
would re-detect taxonomy rather than niche.

## Confounder 4 — phylogenetic non-independence (the central threat)

Phylogenetic signal of niche membership is expected to be high (niche is
clustered on the tree), so an enriched gene may simply mark a human-associated
clade, not a human-selected function.

Controls, applied as three independent methods whose agreement is required:
- **Method A, phyloglm** (`phyloglm_enrichment.R`): phylogenetic logistic
  regression, `feature ~ niche + covariates`, with the GTDB species tree's
  covariance. Models ancestry as a continuous correlation structure.
- **Method B, CMH** (`cmh_stratified.py`): Cochran-Mantel-Haenszel stratified by
  genus, so human species are compared only with non-human species in the same
  genus. Matches ancestry discretely.
- **Method C, balanced bootstrap** (`balanced_resampling.py`): species-level,
  niche- and host-balanced resampling.

A consensus signature requires every method **applicable to that feature** to
agree in direction and pass its threshold (`consensus_signatures.py`). The
distinction between applicable and significant is the crux, and it cuts against
the confounder control rather than with it. CMH is undefined for a feature when
no genus contains species from both niches — which is to say, when the feature
lives in a clade that is itself niche-specific. Requiring CMH unconditionally
would therefore delete the very signatures the study is about, and would do it
silently, reporting zero rather than an error. So a call needs at least two
applicable methods, must include phyloglm (the anchor), and must not be
contradicted by any method that could be run. Calls with all four methods
applicable are tier1; calls where one was untestable are tier2, counted and
labelled separately so the manuscript can quote both numbers and say which method
was missing. `stats.consensus.strict_all_methods` restores the unconditional
rule for a sensitivity analysis.

The phylogenetic signal of niche
and of each top feature is measured explicitly (Fritz-Purvis D, Pagel's lambda,
Blomberg's K) so the threat is quantified, not assumed away (`phylo_signal.R`).
Ancestral-state reconstruction then tests whether a signature was acquired once
or convergently across the tree (`ancestral_convergence.R`).

A sensitivity run turns phylogenetic control off and reports how many more calls
a naive analysis would have made (`config/sensitivity/nophylo.yaml`). The size of
that jump is the quantitative case for the controls.

## Confounder 5 — genome quality / completeness

Annotation counts scale with completeness; an apparent gene gain can be a
completeness difference.

Controls:
- Completeness and contamination are covariates in every model (phyloglm, PGLS).
  Code: `config.yaml stats.covariates`, the R model formulas.
- Species-level prevalence across multiple genomes buffers single-genome
  incompleteness.
- A high-quality-only run (completeness > 90, contamination < 5) is part of the
  sensitivity sweep. Code: `config/sensitivity/hqonly.yaml`.
- QC thresholds are re-applied by this pipeline, not inherited; recorded in
  `results/00_ingest/samples_qc_report.json`.

## Confounder 6 — genome size

Larger genomes encode more functions, so they look "enriched" for many things at
once.

Controls: log10 genome size is a covariate in every model and a separate block
in variation partitioning, so its contribution is partialled out rather than
attributed to niche. Code: model formulas; `ordination_varpart.R varpart`.

## Confounder 7 — compositionality and multiple testing

Controls:
- Presence/absence at the species level sidesteps much of the compositional
  problem; the count-based ordination uses a CLR transform option.
- Benjamini-Hochberg FDR is applied per method, then consensus across three
  methods is required, then a phylogeny-aware permutation null gives an empirical
  FDR (`permutation_null.py`). Effect sizes and CIs are always reported, not just
  p-values.

## Annotation-source confounders

KO is taken from KofamScan (curated HMM thresholds) with eggNOG as a cross-check;
CAZymes from a dedicated dbCAN run (two-tool consensus) with eggNOG as a
cross-check. Per-genome concordance between the two sources is written out
(`*_source_concordance.tsv`) so the choice of primary source is auditable.

## Positive and negative controls

- Positive: host-glycan/mucin CAZyme families (GH33, GH20, GH29, GH95, GH101)
  and a sialidase KO are expected to be human-gut enriched; the report flags them
  if they are not recovered (`consensus_signatures.py`, configurable).
- Negative: the phylogeny-aware label permutation estimates how many signatures
  arise by chance under the same pipeline.

## Within-species niche-transition analysis (`transition_all`)

The "is this niche a recent acquisition" question has its own failure modes, each
controlled:

- Unequal strains per niche distorts diversity, the SFS and Tajima's D. Only
  species with enough near-complete strains in each niche are tested, and every
  population statistic is computed on equal-n subsamples with bootstrapping.
  Code: `transition_select.py`, `popgen_sfs.py`.
- Clonal/epidemic oversampling fakes low diversity or a clade. Strains are
  dereplicated within each niche (Mash); near-identical strains across niches are
  kept because they are the recent-transfer signal. Code: `dereplicate_strains.py`.
- Recombination distorts trees, branch lengths and diversity. The core alignment
  is masked with Gubbins before the tree, the SFS and the diversity statistics.
  Code: `gubbins_species` rule.
- Directionality needs a root. Trees are rooted on a congeneric outgroup; the
  ancestral niche is the reconstructed root state, never assumed. Code:
  `transition_select.py` (outgroup), `ancestral_niche.R`.
- No single statistic is trusted. A recent-acquisition call requires phylogenetic
  structure (Slatkin-Maddison) plus agreement among transition depth, diversity
  reduction, Tajima's D, nestedness and gene gain. Code: `transition_verdict.py`.
- One species could be a fluke. The headline is the cross-species binomial test
  for a consistent transition direction. Code: `transition_meta.py`.

## Functional enrichment of the signatures (`08_enrichment`)

Turning a signature feature list into "the human gut is enriched for pathway X"
is where enrichment analyses are most often done indefensibly. The controls:

- Background is the tested, annotatable feature set for that layer and contrast,
  not the whole database. Testing against the database inflates every category;
  testing against what was actually examined does not. Code: `ora_enrichment.py`.
- Direction is separated (enriched vs depleted) so opposing signals cannot
  cancel into a false null.
- Two methods must agree: a hypergeometric over-representation test and a
  threshold-free preranked GSEA (fgsea) on the signed association statistic. ORA
  depends on the significance cutoff; GSEA does not; requiring both guards
  against artefacts of either. Code: `ora_enrichment.py`, `gsea_enrichment.R`,
  `enrichment_combine.py`.
- Multiple testing is BH-controlled within each system and direction; set sizes
  are bounded; leading-edge features are reported so a call is traceable to
  specific genes.
- Membership comes from canonical sources (KEGG module definitions, KEGG REST
  KO-to-pathway) or deterministic rules (CAZyme class, COG group, BGC class);
  systems with no available mapping are skipped, not guessed.

## What is not fully controlled

Sample provenance (study, DNA extraction, assembler) is partly confounded with
niche and cannot be removed without per-sample metadata. The species-level,
phylogenetically controlled, multi-method design reduces its influence; it is
stated as a limitation in `docs/METHODS.md` rather than hidden.
