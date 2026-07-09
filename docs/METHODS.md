# Methods (manuscript draft)

Edit freely; numbers in braces are filled from the run reports. This text
matches what the pipeline actually does.

## Genome set and quality control

We analysed {n_passed} genomes spanning three niches (human gut, animal gut,
free-living) drawn from {n_input} metagenome-assembled and isolate genomes.
Genomes were retained at completeness >= 50% and contamination <= 5% (CheckM),
quality score (completeness - 5x contamination) >= 50, N50 >= 5 kb and <= 1000
contigs. Taxonomy follows GTDB ({gtdb_release}). A high-quality subset
(completeness > 90%, contamination < 5%) was defined for sensitivity analysis.

## Species units

Genomes were grouped into species by their GTDB species assignment, including
GTDB placeholder/de-novo clusters. The species was the unit of all comparative
analyses. For each species and each functional feature we computed prevalence as
the fraction of conspecific genomes carrying the feature; a feature was present
in a species when prevalence was at least 0.5 (varied 0.1-0.9 in sensitivity
analysis). One representative genome per species (highest quality score, then
N50, then size) was used for computationally heavy de novo annotation and for
mapping species to the phylogeny.

A species was scored as occupying a niche if at least one of its genomes came
from that niche. Single-niche species (specialists) were the basis of the
primary contrast; multi-niche species (generalists) were analysed separately.
Niche breadth was summarised with the standardised Levins index.

## Taxonomic structure by niche

All taxonomic comparisons used the species as the unit (each species counted
once); genome-weighted views were reported alongside only to show strain
sampling bias. Richness was compared across niches by sample-size rarefaction
and by Chao1 and ACE asymptotic estimators with sample-coverage, at matched
effort. Niche breadth of each taxon (phylum to genus) was summarised by the
standardised Levins index, and a taxon's niche restriction was tested against a
null that shuffled species' niche labels (preserving niche sizes), with
Benjamini-Hochberg correction within rank. Indicator taxa for each niche and
niche combination were identified with group-size-corrected IndVal (indicspecies
multipatt, permutation test, FDR within rank). Over- and under-representation of
each taxon in each niche was tested with Fisher's exact test on species counts
(FDR within rank), complemented by a rank-level G-test of independence with
standardised residuals.

Cross-niche species overlap was compared to a null produced by randomising the
species-by-niche occupancy matrix preserving both margins (curveball), reported
as a standardised effect size and empirical p per niche pair and the triple
intersection. Between-niche dissimilarity was partitioned into turnover and
nestedness components (Baselga, Sorensen family). Phylogenetic community
structure per niche was computed on the GTDB scaffold as Faith's PD with its
standardised effect size, and as NRI and NTI (from ses.mpd and ses.mntd) against
the configured null model; phylogenetic beta diversity between niches used
UniFrac. Taxonomic novelty (GTDB placeholder lineages) was quantified per niche
and rank, rarefied to a common effort so it did not reflect sequencing depth.
The animal niche was resolved by host species: per-host richness was rarefied to
a common genome count, and animal richness, novelty and composition were
recomputed with and without the dominant host (mouse) to show its contribution.

## Within-species niche transitions

For each species present in more than one niche with at least ten near-complete
genomes (completeness >= 90%, contamination <= 5%) in each of at least two
niches, we tested whether a niche was a recent acquisition. Strains were
dereplicated within each niche with Mash (collapsing pairs above ~99.99% ANI,
keeping the highest-quality representative) so clonal oversampling could not bias
diversity; near-identical strains spanning niches were retained as they indicate
recent transfer. A congeneric species present in the dataset (the nearest sister
on the GTDB scaffold) was added as an outgroup for rooting.

The core genome was built with Panaroo and aligned, recombination was removed
with Gubbins, and a maximum-likelihood tree was inferred from the recombination-
filtered sites with IQ-TREE (ModelFinder, ultrafast bootstrap and SH-aLRT,
ascertainment-bias correction), rooted on the outgroup. Niche was reconstructed
as a discrete trait on the rooted ingroup tree under an Mk model (equal-rates vs
all-rates-different by AIC) with stochastic character mapping, giving the
ancestral (root) niche, directed transition counts and the depth of each
transition (shallow transitions indicate recent acquisition). A Slatkin-Maddison
permutation test confirmed that niche was phylogenetically structured before
directionality was interpreted.

Within each niche population, on the recombination-masked SNPs and after
subsampling to equal strain number with bootstrapping, we computed nucleotide
diversity, Watterson's theta, the folded site-frequency spectrum, Tajima's D and
the singleton fraction; a recently founded niche is expected to show reduced
diversity and a negative Tajima's D from post-founder expansion. Between niches
we computed Hudson's Fst, dxy and net divergence da, and inferred the direction
of derivation from three statistics measured at equal sample size: nestedness of
segregating variation (the fraction of a population's segregating sites at which
its full allele set is contained in the other), the number of private alleles,
and nucleotide diversity. A derived population carries a subset of the source's
variation, so it is the more nested, the less diverse, and the one with fewer
private alleles. Each statistic is antisymmetric under exchanging the two
populations; confidence intervals came from bootstrap subsampling and
significance from permuting the niche labels across the pooled genomes, since a
bootstrap within fixed populations reproduces whatever asymmetry those samples
contain and cannot test the null. A pair with no significant asymmetry was
reported as unresolved rather than assigned a direction.

Accessory genes differentially present between niches (Fisher's exact test,
FDR-corrected) identified gene content gained with the new niche. Because the
power of that test depends on sample size and the two niche populations of a
species are rarely the same size, the per-niche counts of differentially present
genes were taken from draws in which both populations were subsampled to the
smaller; with equal row margins the exact two-sided probability is twice the
smaller hypergeometric tail, and a gene was counted when it was significant in at
least half of the draws. An optional coalescent model
(fastsimcoal2 on the joint folded SFS; see resources) estimates split time,
effective sizes and directional migration.

A niche was called a recent acquisition for a species only when niche was
phylogenetically structured and the evidence agreed across lines: shallow
transitions, reduced diversity, negative Tajima's D, nestedness, and gene gains.
Finally, across all tested species, we tested whether supported acquisitions ran
in a consistent direction with a binomial sign test, since repeated independent
acquisitions in one direction are far stronger evidence than any single species.

## Functional annotation

Existing per-genome annotations (Prokka gene calls, eggNOG-mapper v2, KofamScan)
were ingested. KEGG orthologs were taken from KofamScan (significant assignments)
with eggNOG KO as a concordance cross-check. Pfam, COG categories and EC numbers
were taken from eggNOG. CAZymes were annotated de novo with run_dbcan (families
supported by at least two of HMMER, dbCAN-sub and DIAMOND), with eggNOG CAZy as a
cross-check. Biosynthetic gene clusters were predicted with antiSMASH on species
representatives, and antimicrobial resistance genes with AMRFinderPlus. KEGG
module completeness was computed from species KO sets by evaluating each module's
KEGG definition (stepwise block completeness); a module was present at >= 0.66
completeness. Per-genome source concordance (Jaccard) was recorded for KO and
CAZyme.

## Diversity

Taxonomic composition was reported species-weighted and genome-weighted to expose
strain-sampling bias. Species richness was compared across niches by sample-size
rarefaction and by Hill numbers (q = 0, 1, 2) with the Chao coverage estimator,
both at matched sampling effort.

## Differential function

For each functional layer and each contrast (human vs all others; human vs
animal gut), features reaching 5% prevalence in at least one group were tested
with three methods sharing identical inputs:

1. Phylogenetic logistic regression (phylolm::phyloglm, logistic MPLE) of feature
   presence on niche with completeness, log10 genome size and GC as covariates
   and the GTDB species-tree covariance.
2. Cochran-Mantel-Haenszel test stratified by genus, giving a clade-matched
   common odds ratio.
3. Balanced bootstrap (1000 iterations) drawing equal species per niche and one
   genome-equivalent per species, with animal species weighted by inverse host
   frequency; the prevalence difference and Haldane-corrected odds ratio were
   recorded.

Effect sizes were expressed as log2 odds ratios, clipped to +/-10 so that
complete separation contributes a bounded, signed estimate rather than an
infinite one or a discarded feature. p-values were FDR-corrected
(Benjamini-Hochberg) within each method.

Not every method can be evaluated for every feature. The Cochran-Mantel-Haenszel
test is undefined when no clade at the stratifying rank contains species from
both groups, Scoary2 discards invariant features, and the phylogenetic regression
occasionally fails to converge; in each case the method is inapplicable rather
than negative, and treating the two as equivalent would suppress precisely those
features restricted to niche-specific clades. A feature was therefore called a
consensus signature when every method applicable to it met its threshold
(FDR < 0.05 and |log2 OR| >= 1 for the regression, the stratified test and the
pan-GWAS; a bootstrap confidence interval excluding zero with sign consistency
>= 0.95 for the balanced resampling), when at least two methods were applicable,
when the phylogenetic regression was among them, and when all supporting methods
agreed in direction. Features supported by all four methods and features for
which one method was inapplicable are reported and counted separately, and the
inapplicable method is named for each feature.

A phylogeny-aware permutation null (niche labels shuffled within phylum, 1000
permutations) provided an empirical false-discovery estimate; the number of
species that the within-clade restriction actually leaves free to move is
reported alongside it, because a null that cannot move labels cannot calibrate
anything. Host-glycan CAZymes and a sialidase KO served as positive controls.

## Pan-GWAS (Scoary2)

As an independent fourth association method we ran Scoary2 (Roder et al. 2024),
whose pairwise-comparisons algorithm corrects for clonal and phylogenetic
structure. Scoary2 was applied at the species level on the same functional
presence matrices (one species per row; features = KO, KEGG modules, Pfam, COG,
CAZyme, BGC and AMR presence), with the GTDB species tree supplying the pairwise
correction and niche as the binary trait. We tested human, animal and free each
against the rest and host-associated (human+animal) against free-living. A
feature was a Scoary2 hit when the Benjamini-Hochberg-corrected Fisher q, the
label-switching permutation p (which flags merely lineage-specific signals), the
number of supporting pairwise comparisons and the best-pairwise p all passed.
For the contrasts shared with the differential trio, Scoary2 entered the
consensus as a fourth method, so a reported signature is supported by a
phylogenetic regression, a clade-stratified test, a balanced bootstrap and a
pan-GWAS, and the concordance among methods is reported.

We also ran Scoary2 clade-stratified: within each well-sampled genus we built a
Panaroo pangenome from species representatives and tested niche association on
real ortholog families, then summarised how many genera show niche-associated
families. This complements the functional-category analysis with gene-family
resolution where orthology is meaningful.

## Functional enrichment of niche signatures

The niche-signature features were tested for over-representation of higher-level
functional categories: KEGG pathways and modules for KOs, CAZyme classes,
COG super-groups, biosynthetic-cluster classes, and (where mapping files are
supplied) GO terms, Pfam clans, CAZyme substrates and AMR drug classes. Module
membership was taken from the KEGG module definitions; KO-to-pathway from the
KEGG REST API; CAZyme class, COG group and BGC class are determined directly from
the feature identifier.

Two methods were applied and required to agree. Over-representation analysis used
the hypergeometric test with the features actually TESTED for that layer and
contrast (restricted to those annotatable in the system) as the background, run
separately for features enriched and depleted in the niche so opposing signals
could not cancel; p-values were Benjamini-Hochberg corrected within each system
and direction, over every category passing the size filter including those with
no overlapping signature feature, since those categories were tested hypotheses
and omitting them would shrink the correction denominator. Preranked gene-set
enrichment (fgsea) used the signed association
statistic (sign of the consensus log odds ratio times -log10 of the phylogenetic-
regression q) across all tested features, which uses the full ranking without a
significance cutoff. A category was reported as enriched when significant by both
methods; single-method results are retained at a lower tier. Results are
summarised as a per-contrast enrichment dot plot and a category-by-contrast
heatmap of signed enrichment.

## Functional redundancy

Functional redundancy of the human gut was measured with the species as the unit
so over-sequenced species do not inflate it. For each functional layer we
recorded the occupancy of every function (how many species carry it) and the
fraction that are core versus rare; built rarefied accumulation curves of unique
functions against species sampled (functions saturating while species keep
accumulating is the redundancy signature); computed the Ricotta functional-
redundancy index (Gini-Simpson diversity D minus Rao's quadratic entropy Q on
functional Jaccard distances, FR = D - Q, relative FR = 1 - Q/D, on a species
subsample with bootstraps); and recorded the phylogenetic breadth of each
function's carriers as the number of distinct families. The same measures were
computed within the Western and non-Western populations when available.

## Western versus non-Western guts

To ask whether the functional repertoire is conserved despite species turnover,
taxonomic and functional turnover were placed side by side. Species composition
dissimilarity between populations (Sorensen, partitioned into turnover and
nestedness; rarefied to equal species number) was compared with the dissimilarity
of the function sets per layer. For every function present in both populations we
computed the Jaccard overlap of its carrier species across populations; low
overlap means the same function is supplied by different species (functional
redundancy expressed as taxonomic substitution). The functions that genuinely
differ were identified by adding a within-human western_vs_nonwestern contrast
that flows through the same four-method consensus and the enrichment analysis, so
the divergent functions are called with the full confounder control. These
analyses activate when the western_nonwestern lifestyle label is present in the
metadata.

## Community and metabolic interaction

Per-genome metabolic and phenotype predictions (amino-acid biosynthesis, carbon-
source utilisation, KEGG module presence, and predicted traits) were analysed at
two levels, with both confounders controlled consistently. Each genome is a
strain, so capabilities were collapsed to species-level prevalence (present when
carried by at least half of a species' strains); amino-acid biosynthesis was
inverted so the feature is auxotrophy (the species cannot make the amino acid).
At the microbial level, auxotrophy, carbon use, traits and modules were treated
as feature layers and tested for human-gut enrichment through the same
phylogenetic-consensus stack as the other functional layers.

At the community level, each niche's species pool was rarefied to the same number
of species (bootstrapped) before any pool metric, so the species-per-niche
imbalance cannot drive the result. We computed auxotrophy dependency (mean
auxotrophies per species and the per-amino-acid prototroph-to-auxotroph ratio,
which is the pool's capacity to provision auxotrophs), byproduct cross-feeding
(for each exchange metabolite the potential producer-to-consumer links, with
acetate producers taken from a specific module and other organic-acid producers
from a fermenter-trait proxy), metabolic division of labour (collective KEGG-
module coverage relative to the mean per species), and community trait composition
(the fraction of the pool that is anaerobic, spore-forming, bile-susceptible,
motile and so on). Niches were compared on all of these.

## Multivariate structure

Species functional repertoires (Jaccard distance on presence/absence, balanced-
subsampled per niche) were ordinated by principal coordinates analysis.
Differences among niches were tested by PERMANOVA (adonis2, marginal terms, with
completeness, genome size and GC as covariates) and checked for dispersion
artefacts with betadisper. Variation partitioning apportioned functional variance
among niche, phylogeny (leading principal coordinates of patristic distance),
genome size and quality; the niche-unique fraction estimates differentiation not
attributable to ancestry.

## Phylogenetic comparative analysis

Species were placed on the GTDB reference trees (bac120, ar53) via their species
representatives; bacteria and archaea were analysed separately, with bacteria as
the primary analysis. Phylogenetic signal of niche membership and of top
signature features used the Fritz-Purvis D statistic; continuous functional
traits used Pagel's lambda and Blomberg's K. Continuous functional load
(CAZyme/BGC/AMR richness, genome size, CDS count) was modelled by phylogenetic
generalised least squares. For the top human-enriched signature features,
stochastic character mapping (100 simulations, all-rates-different model)
estimated the number of independent gains; a feature gained at least three times
independently was scored convergent.

## Reproducibility

The workflow is implemented in Snakemake with per-rule conda environments. All
parameters are in a single config file; all stochastic steps use seeds derived
deterministically from a global seed. Provenance (seed, timestamps, parameters)
is written beside each major output. Sensitivity analyses re-ran the differential
stage under alternative species-presence thresholds, a high-quality-only genome
set, exclusion of the dominant animal host, and without phylogenetic correction.

## Limitations

Sample provenance (study, extraction protocol, assembler) is partly confounded
with niche and was not directly modelled; the species-level, phylogenetically
controlled, multi-method design limits but does not eliminate its influence.
Phylogenetic tests were restricted to species placeable on the GTDB reference
trees (coverage reported per niche). GTDB placeholder species are valid clusters
but not cross-database named species; cross-niche overlap was reported separately
for named and placeholder clusters.
