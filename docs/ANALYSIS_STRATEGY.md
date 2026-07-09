# Analysis strategy: from data to the six questions

This maps the manuscript questions to the evidence the pipeline produces. The
questions are not answered one analysis each. They share a backbone — the
species-level, phylogenetically controlled signatures — and the "why/what drives"
questions add mechanistic analyses on top. The synthesis stage
(`snakemake synthesis_all`) integrates everything into catalogues, a master
figure and a report.

Two groups: the descriptive questions (which species/genes/functions are
human-specific) and the mechanistic ones (what drives adaptation and speciation,
and why the human gut). The descriptive answers are inputs to the mechanistic
ones.

## Which species are human-specific?

Three independent lines, integrated into one ranked catalogue with a confidence
tier (`synthesis/human_specific_species.tsv`):
  * niche specialists — species occupying only the human niche (occupancy on the
    QC-passed, species-collapsed data);
  * indicator species — IndVal.g, which combines exclusivity and fidelity, with a
    permutation p and FDR;
  * the permutation specificity test — a species' niche restriction tested against
    a label-shuffling null.
A species is called human-specific with high confidence when all three agree.
Novelty (undescribed lineages) and strain counts are carried along for context.
This is defensible because each line controls a different failure mode (raw
counts, fidelity, sampling), and the species, not the genome, is the unit
throughout.

## Which genes are human-specific?

Two resolutions. At the functional-feature level, the four-method consensus
(phylogenetic logistic regression, clade-stratified CMH, balanced bootstrap,
Scoary2) calls features enriched or depleted in human, with every confounder
(strain imbalance, species-per-niche, host imbalance, phylogeny, quality, genome
size) controlled. At the true ortholog-family level, clade-stratified Scoary2 and
the within-species accessory-genome analysis identify gene families gained with
the human niche. Both feed `synthesis/human_specific_features.tsv`. Cross-database
"genes" across all 60k species are deliberately not attempted, because ortholog
families do not span that range; the functional features are the consistent
cross-genome currency, with ortholog families resolved within clades.

## Which functions are human-specific?

Enrichment lifts the gene-level signatures to functional categories: KEGG
pathways and modules, CAZyme classes, COG groups, BGC classes, tested by both
over-representation (proper tested-feature background) and preranked GSEA, with
agreement required. The curated ecological-pressure sets (below) add hypothesis-
driven function groups. These populate `synthesis/human_specific_functions.tsv`.

## What drives niche adaptation?

Adaptation mode decomposition (`adaptation_mode.py`) classifies each human
signature by HOW it arises: gene acquisition (present in human species, largely
absent elsewhere — gain/HGT), gene loss (present elsewhere, reduced in human),
or sequence-level change (the gene is shared but evolves differently). Gain/loss
come from the species-level prevalence contrast; sequence-level change comes from
the HyPhy selection scan (BUSTED for episodic positive selection, RELAX for
intensified selection on human-associated branches) run on the human-associated
ortholog families from the clade pangenomes. Genome architecture (size, GC,
coding density) is tested as a parallel driver by PGLS. The output states the
relative contribution of acquisition, loss, selection and architecture — that
breakdown is the answer, not any single mechanism.

## What drives functional adaptation and/or speciation?

Functional adaptation: the adaptation-mode breakdown above, plus the variation
partitioning that already attributes functional variance to niche after
phylogeny, genome size and quality. Speciation: the diversification-rate analysis
(`diversification_rate.R`) computes the per-tip DR statistic on the GTDB scaffold
and tests whether human/host-association is associated with faster or slower
diversification than free-living, with a phylogeny-aware permutation. Linking the
two, the synthesis reports whether the niches that diversify faster are also those
with more gene gain, i.e. whether functional acquisition accompanies
diversification.

## Why have these species specifically adapted to the human gut?

The "why" is made testable rather than narrated. Curated gene sets representing
known human-gut selective pressures — host-glycan and mucin foraging, bile
tolerance, oxidative-stress handling for a near-anaerobic gut, short-chain fatty
acid fermentation, and dietary-fiber and starch degradation — are tested for
human enrichment (`ecological_pressure_test.py`) using the same defensible
over-representation framework. The pressures whose gene sets are human-enriched,
and which also carry signature genes under positive selection, are the candidate
reasons these species occupy the human gut. The sets are curated and editable;
they encode established hypotheses so the test is explicit and falsifiable.

## Integration

`synthesis_catalogues.py`, `fig_synthesis_master.py` and the synthesis report
combine the above into the three catalogues, a single master figure (counts of
human-specific species/genes/functions, the adaptation-mode breakdown, the top
human-gut pressures, and diversification by niche), and a document that answers
each question with its evidence and links to the supporting tables and figures.
