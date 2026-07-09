# Audit: what each test measures, why it is defensible, and how it was verified

This document exists so a reviewer can check the pipeline without running it. For
every analysis it states the hypothesis, the statistic, the confounder controls,
the known limitations, and the verification status. It also records the bugs
found during the audit and how each was fixed, because "we looked and found
nothing" is a weaker claim than "we looked, found these, and fixed them".

Reproduce the verification with:

```bash
bash tests/run_tests.sh
```

which compiles every Python script, validates the config, checks that every
script referenced by a Snakemake rule exists, runs `tests/test_statistics.py`
(hand-derived numerical checks) and `tests/test_pipeline_smoke.py` (end-to-end
runs on synthetic data with planted signals).

## Verification status, by component

Three levels are used throughout the tables below.

**Executed** means the code was run on synthetic data in which the answer is
known by construction, and the test asserts that the answer comes back.
**Hand-checked** means the statistic was computed by hand or reduced to an
analytic limit and compared against the implementation. **Reviewed** means the
code was read line by line against the method's definition but not executed here,
because the environment has no R, no Snakemake and none of the external tools;
those run on your cluster.

Everything in `scripts/python/` outside the figure scripts is Executed or
Hand-checked. Everything in `scripts/R/` is Reviewed. The external tools
(GTDB-Tk, Mash, Panaroo, Gubbins, IQ-TREE, Scoary2, HyPhy, dbCAN, antiSMASH,
AMRFinder, fgsea) are Reviewed at the level of their invocation and output
parsing; their internals are the responsibility of their authors.

## Stage 1 — species units and quality control

| What | Statistic | Defensibility | Verified |
|---|---|---|---|
| Genome QC | completeness ≥ 90 %, contamination ≤ 5 %, contig and N50 floors | thresholds are the MIMAG "high-quality draft" convention, applied before anything else so quality never enters a model as signal | Executed |
| Species unit | GTDB species; strains collapse to a species-level prevalence vector | the unit of ecological and evolutionary inference is the species, not the assembly. Strains-per-species is the largest imbalance in the metadata (some species have thousands of genomes, most have one) and collapsing removes it at the source rather than adjusting for it later | Executed |
| Niche assignment | a species is a specialist when ≥ `species.specialist_threshold` of its genomes come from one niche | species with genuinely mixed occupancy are held out of the specialist contrasts and analysed separately in the transition module, instead of being forced into one niche | Executed |
| Population label | `western_nonwestern` metadata column | the whole population module stays inert until the column exists, so the pipeline cannot silently analyse an absent variable | Executed |
| Annotation paths | one resolved table, from a manifest of explicit paths or from `{genome}` templates | a genome with no path is a hard error, never a skipped genome: silently dropping one biases every prevalence estimate. Declared paths are stat-ed before the first job is submitted. Mash and Panaroo see a symlink farm of `<genome>.fna` / `<genome>.gff`, so an arbitrary basename cannot become a tree tip label | Executed |

## Stage 2 — taxonomy and ecology

| What | Statistic | Defensibility | Verified |
|---|---|---|---|
| Niche specificity of a taxon | standardised Levins' niche breadth `B_std`, 0 = one niche, 1 = uniform | breadth is computed on the species-per-niche matrix, not the genome matrix, so a heavily sequenced species cannot make its genus look niche-specific | Hand-checked (`test_levins`) |
| Is that specificity more than expected? | permutation of niche labels across species within rank; one-sided empirical p `(1 + #{B_null ≤ B_obs}) / (N + 1)` | the add-one is required: without it a p of exactly 0 is reportable, which is not a probability | Executed |
| Taxon enrichment per niche | Fisher's exact test per taxon × niche, Haldane-corrected log2 OR, BH across all cells | BH runs over every taxon × niche cell tested, not the significant subset | Executed |
| Cross-niche overlap | observed shared species vs a curveball (fixed-fixed) null | the curveball preserves BOTH margins — species per niche and niches per species — so overlap is compared against communities with the same richness and the same generalist/specialist mix. A fixed-degree null preserves only one margin and would report overlap that is a consequence of richness | Hand-checked (`test_curveball` asserts both margins invariant) |
| Novelty | fraction of species with placeholder GTDB names, rarefied to equal sampling effort | novelty scales with the number of genomes examined; without rarefaction the best-sampled niche always looks the most novel | Executed |
| Richness / diversity | Chao1 (bias-corrected), ACE, Hill numbers q = 0,1,2, coverage-based rarefaction | Hill numbers make richness, Shannon and Simpson one parametric family, so they cannot be cherry-picked | Hand-checked (`test_richness`) |
| β-diversity partition | Baselga's Sørensen decomposition into turnover and nestedness | separates "different species" from "fewer species", which is the difference between niche specialisation and sampling depth | Hand-checked (`test_baselga`) |
| Indicator taxa | IndVal.g with restricted permutations (`permute::how`) | the `.g` correction removes the group-size bias in the original IndVal, which matters because the niches have very different species counts | Reviewed |
| Phylogenetic community structure | Faith's PD, NRI, NTI against a taxa-shuffle null on the GTDB scaffold | states whether a niche's species are more closely related than chance, which is the taxonomic prerequisite for any claim about niche-specific clades | Reviewed |

## Stage 3 — differential function (the consensus)

Four methods, each controlling a different confounder, must agree.

| Method | What it controls | Statistic |
|---|---|---|
| `phyloglm` (phylolm) | shared ancestry, as a continuous covariance | phylogenetic logistic regression of presence on niche + completeness + log10 genome size + GC |
| CMH | shared ancestry, discretely: human species are only ever compared with non-human species in the same genus | Mantel–Haenszel common odds ratio and CMH χ² |
| Balanced resampling | species-per-niche AND strains-per-species imbalance, plus host imbalance in the animal niche | 1000 draws of equal species per niche, one genome per species, dominant host capped; bootstrap CI of the prevalence difference and sign consistency |
| Scoary2 | population structure, via the pairwise-comparisons (Maddison) correction | Fisher q, label-switching empirical p, supporting pairwise comparisons |

The consensus rule and its justification are in the module docstring of
`scripts/python/consensus_signatures.py`. The point worth restating here:
**a method that cannot be applied to a feature is not a method that disagrees.**
CMH is undefined when no genus contains both groups, Scoary drops invariant
features, phyloglm can fail to converge. Requiring all four unconditionally
would silently produce zero signatures precisely where the biology is
strongest — features confined to niche-specific clades — and would leave no
record of having done so. So: every *applicable* method must agree, at least
`min_applicable_methods` (2) must be applicable, and the anchor method
(`phyloglm`) must be applicable and significant. Calls where all four applied
are `tier1_consensus`; calls where one was untestable are
`tier2_consensus_partial`, counted separately, with `methods_untestable` naming
the method per feature. Set `stats.consensus.strict_all_methods: true` to report
only tier1.

Effect sizes: complete separation gives an infinite odds ratio. Rather than drop
those features (which would discard the strongest signal) or report infinity,
every method's log2 OR is clipped to ±10 with the direction preserved, and the
CMH point estimate falls back to a Haldane-corrected pooled OR when the
Mantel–Haenszel estimate is 0 or ∞. The *p*-value is unaffected.

Null model: `permutation_null.py` shuffles niche labels within clade, so the
clade–niche association is preserved and the null is conservative. It now also
reports how many species were actually permutable; if a stratum holds only one
group its labels are invariant, and if that is true of almost every stratum the
null is degenerate and its empirical FDR means nothing. That case is detected and
warned about rather than reported as a number.

Verified: Executed. `tests/test_pipeline_smoke.py` plants one feature present in
all six human species and absent from all twelve others, plus a background
feature present everywhere, and asserts that (a) the planted feature is called
and labelled `human_enriched`, (b) the background feature is not called and its
direction is `ns` rather than a spurious depletion, (c) CMH is correctly reported
as untestable in a design where every species is niche-specific, (d) an
applicable-but-null CMH vetoes the call, and (e) an untestable anchor blocks
every call. CMH orientation (which niche the odds ratio points at) is checked
against a hand-computed pooled OR of 16 in `test_cmh_orientation`.

## Stage 4 — within-species niche transition

The question is whether one niche is a recent acquisition for a species that
occupies more than one. No single statistic answers it; five lines of evidence
are combined in `transition_verdict.py`, behind a gate.

| Line | Statistic | Note |
|---|---|---|
| gate | Slatkin–Maddison test that niche is phylogenetically structured | if niche is randomly scattered across the tree there is no transition to date |
| e1 | mean depth of niche transitions from stochastic mapping (Mk ER vs ARD by AIC) | shallow transitions = recent |
| e2 | π lower in the derived niche than the source | founder effect |
| e3 | Tajima's D negative in the derived niche | post-founder expansion |
| e4 | the derived niche's variation is nested inside the source's | directionality |
| e5 | accessory genes gained in the derived niche | colonisation cargo |

Tier: strong = gate + e4 + ≥ 3 of the rest; moderate = gate + ≥ 2; else weak.

Sequence handling: Mash dereplication is applied *within niche only*, so a clonal
expansion cannot inflate one niche's diversity while cross-niche differences are
preserved. Gubbins masks recombinant tracts before the tree is built, because
recombination violates the tree model that everything downstream assumes.
Rooting uses the nearest sister species in the dataset; without a root there is
no ancestral state and therefore no direction.

Population genetics: π, Watterson's θ, the folded SFS, Tajima's D, and Hudson's
Fst/dxy/da, all computed on bootstrap subsamples of equal size, because S and the
private-allele count both grow with n and would otherwise track sample size
rather than biology.

Directionality (`demography_directionality.py`) was rewritten during this audit;
see the bug list below. It now measures three antisymmetric statistics at equal
sample size — allele-set containment (nestedness), private-allele count, and π —
takes CIs from a bootstrap and *significance from a label permutation*, and
returns `unresolved` rather than a coin flip when nothing separates the two
populations.

Verified: Executed. The smoke test builds a source population of 20 haplotypes,
founds a derived population from 3 of them, expands it to 16 with drift and 3
novel mutations, and asserts that π falls, that S falls at equal n, that the
derived population is called correctly with the permutation p significant, that
all three lines agree, and that it carries fewer private alleles. It then
shuffles the labels and asserts that no directional call is made. Nestedness,
private alleles and π are additionally checked against a hand-worked four-site
example in `test_nestedness`; Tajima's D against Tajima's own worked n = 4, S = 2
case (D = 0.59142) in `test_tajima`; Hudson's Fst against its two analytic limits
in `test_hudson`.

Gene gain/loss (`accessory_differentiation.py`) was also rewritten. Fisher's
exact test on the full data is valid, but its *power* depends on sample size, so
with unequal strain counts more genes clear significance in the larger
population, and "genes gained in X vs Y" — evidence line e5 — is biased toward
whichever niche has more strains. The counts now come from balanced draws with
both populations subsampled to the smaller; with equal row margins the exact
two-sided p is exactly twice the smaller hypergeometric tail, which is verified
against `scipy.stats.fisher_exact` for every table with n = 3…10 in
`test_balanced_fisher`. The unbalanced full-data counts are retained in the
output under `unbalanced_*` for reference.

## Stage 5 — functional enrichment

Over-representation analysis and preranked GSEA must agree before a category is
called (`confidence == "both"`).

ORA background is the set of features actually **tested** for that layer and
contrast that carry at least one annotation in the system — not the whole KEGG or
CAZy database. Using the database as background inflates every enrichment by the
annotation rate. Directions are tested separately (up, down, all): pooling them
lets opposing signals cancel. BH is applied over every category passing the size
filter, **including categories with no overlap** (p = 1), because those were
hypotheses too; the common practice of dropping them before correction shrinks
the denominator and inflates significance. Only the k ≥ 1 rows are written out.

Preranked GSEA (fgsea, multilevel) is run on the consensus log2 OR, so it uses
the full ranking rather than a thresholded list and can see coordinated shifts
that no single feature would show.

Curated ecological-pressure sets (`resources/genesets/human_gut_pressures.tsv`,
43 entries across host glycan/mucin foraging, bile tolerance, oxidative stress,
SCFA fermentation, dietary fibre, starch/glycogen) are tested the same way. They
are a hypothesis, not a result: they were assembled by hand and are meant to be
edited. Provenance for each entry is in the file.

Verified: Executed. The smoke test constructs 20 tested features, 10 GH and 10
GT, makes all 8 human-enriched signatures GH, and asserts the background is the
20 tested features, the fold enrichment is exactly 2, and the p equals the
hand-computed hypergeometric 45/125970. Orientation (that ORA is testing
over- and not under-representation) is asserted separately in
`test_ora_orientation`.

## Stage 6 — synthesis, redundancy, populations

| What | Statistic | Defensibility | Verified |
|---|---|---|---|
| Human-specific species / genes / functions | tiered catalogues, each row carrying the evidence that produced it | a species is "human-specific" for a reason that is printed next to it, not by assertion | Executed |
| Adaptation mode | acquisition vs loss, from prevalence in the focal niche versus the comparator | distinguishes gene gain from gene loss, which have different evolutionary interpretations and are otherwise conflated by a two-sided test | Executed |
| Diversification rate | Jetz's DR statistic per niche, within-phylum permutation null | asks whether the human gut holds more species because lineages there diversify faster, or because more lineages entered | Reviewed |
| Functional redundancy | Ricotta's FR = D − Q (Gini–Simpson minus Rao's quadratic entropy), plus function-accumulation curves | FR is the part of species diversity not explained by functional dissimilarity, which is exactly "how many species do the same job" | Hand-checked (`test_ricotta`: relFR = 1 for identical repertoires, 0 for disjoint) |
| Western vs non-Western | taxonomic Sørensen versus functional Sørensen, per layer, plus carrier-Jaccard for shared functions | the question is whether different species supply the same functions. High taxonomic turnover with low functional turnover and low carrier overlap is exactly the carrier-substitution signature | Executed |

The redundancy denominator was a bug and is now fixed: every quantity for a layer
is computed on the same set of species, namely the species annotated in that
layer, and both `n_species_pool` and `n_species` are written to the output so the
denominator is visible.

## Stage 7 — interacting communities

Built from the files in `input/`: amino-acid biosynthesis, carbon utilisation,
KEGG module completeness, predicted traits.

| What | Statistic | Confounder control |
|---|---|---|
| Auxotrophy | fraction of amino acids a species cannot synthesise | prevalence per species, then rarefaction of species per niche |
| Cross-feeding | producer → consumer byproduct links (acetate, lactate, propionate, succinate, ethanol) | link density computed on rarefied species pools of equal size, since edge counts grow superlinearly with pool size |
| Module complementarity | collective module completeness of a rarefied pool, versus mean per-species completeness | the gap between the two is the community-level function that no single species carries |
| Predicted traits | prevalence per niche | rarefied |

Every community statistic is computed on pools of equal species number, drawn
`community.rarefaction.bootstrap` times. Without that, the human niche — which
has the most species — would win every count.

The byproduct table (`resources/genesets/metabolite_exchange.tsv`) is honest
about its own weakness: acetate production is grounded in KEGG module M00579;
D-/L-lactate, propionate, succinate, ethanol, pyruvate and fumarate use the
"Glucose fermenter" predicted trait as a proxy. Those rows are flagged in the
file and are meant to be replaced with pathway-level evidence before the
cross-feeding results are published as more than exploratory.

Verified: Executed for ingest and the auxotrophy gradient (human species are
asserted to be more auxotrophic than free-living ones, from planted data). The
species-ID join against the metadata was checked on the real input files: 0 of
342,759 genome identifiers failed to match.

## Bugs found and fixed during this audit

Crashes:

1. `ancestral_niche.R` called `parsimony()` and `phyDat()` without attaching
   `phangorn`; `indicator_species.R` called `how()` without attaching `permute`.
   Those two packages, plus `jsonlite` and `yaml` — the last attached by
   *fourteen* of the sixteen R scripts to read the config — were missing from
   `envs/r.yaml`, so every R rule would have failed on first contact with a clean
   conda environment. `tests/run_tests.sh` now cross-checks every `library()` call
   in `scripts/R/` against `envs/*.yaml` so this class of error cannot recur.
2. `ancestral_niche.R` dropped tip names when building `phyDat` from a vector, so
   the permutation compared states against the wrong tips. Fixed by constructing
   a named single-column matrix.
3. `kegg_module_completeness.py` emitted an extra `name` column, so
   `combine_presence` crashed concatenating layers with different schemas. Names
   are now written to a sidecar `_names.tsv` and the parquet schema is uniform;
   `concat_parquet.py` additionally takes the column intersection.
4. `make_report.py` and `make_taxonomy_report.py` formatted `'NA'` with `{:,}`,
   which raises `ValueError`.
5. Three missing-column guards in `enrichment_aggregate.py`,
   `synthesis_catalogues.py` and `fig_synthesis_master.py` for the case where a
   layer × contrast produced no enrichment.

Silent statistical errors — the ones that matter, because they produce a number
rather than a traceback:

6. `overlap_nullmodel.py` read `null_model: fixed_fixed` from the config but only
   branched on the string `curveball`, so it ran the weaker fixed-degree null and
   reported it as fixed-fixed. Both spellings now select the curveball, and the
   null actually used is written to the output.
7. `community_interactions.py` dropped species with no annotated feature when
   subsetting a rarefied pool, shrinking the denominator and inflating every
   per-species rate. Now reindexed with zero fill, so the denominator is the pool
   size.
8. `functional_redundancy.py` used different species sets for occupancy, the
   accumulation curve and Rao's Q within the same layer. All three now use the
   species annotated in that layer, and the counts are reported.
9. `consensus_signatures.py` treated an undefined method as a failed method,
   which made the four-method consensus unsatisfiable for exactly the features
   confined to niche-specific clades. Applicability and significance are now
   distinct; see Stage 3.
10. `cmh_stratified.py` returned `p = NA` both when the test failed and when no
    stratum was informative, and gave up on the point estimate under complete
    separation. It now reports `status` explicitly and falls back to a
    Haldane-corrected pooled OR for the effect size, keeping the valid p.
11. `scoary_parse.py` applied `.fillna(1)` to `empirical_p`. Scoary2 writes that
    column only when run with permutations, so a run without them made every gene
    non-significant and silently muted the fourth consensus method. Criteria are
    now applied only when the column exists, the applied set is written to
    `*.criteria.txt`, and a run with no structure-aware criterion is flagged.
12. `demography_directionality.py` measured nestedness as "is Y's *minor* allele
    present in X". Which allele is minor depends on frequency, not ancestry, so
    the statistic was not a nestedness measure. Replaced with allele-set
    containment. It also used a bootstrap CI as a significance test: resampling
    within fixed populations reproduces whatever asymmetry those samples contain,
    so a CI excluding zero says nothing about the null. Significance now comes
    from a label permutation. Both errors are caught by the shuffled-label test,
    which the old code failed at p = 0.02.
13. `accessory_differentiation.py` compared gene-gain counts between populations
    of unequal size using a test whose power depends on size. Counts are now taken
    from balanced draws; see Stage 4.
14. `ora_enrichment.py` excluded zero-overlap categories before BH, shrinking the
    multiple-testing denominator.
15. `permutation_null.py` reported an empirical FDR that could exceed 1 and gave
    no permutation p; and it could not tell a conservative null from a degenerate
    one. It now caps the FDR, reports `p_empirical` with the add-one correction,
    and flags the degenerate case.
16. `check_sge_profile.py` — the cluster-profile validator written to catch this
    class of error — committed it. `qconf -sq` is permitted only from an admin
    host, and from anywhere else answers `denied: host "..." is not an admin
    host`. The validator read that stderr and reported `queue 'standard.q' not
    found`, turning a refusal to answer into the answer "no". It now separates
    VERIFIED from ASSUMED (recorded in `config/sge/cluster.yaml`, with
    provenance) from UNVERIFIABLE, never reports absence from a denied query,
    cross-checks recorded facts against qconf when qconf answers, and falls back
    to `qsub -w v` — a dry scheduling run that works from any submit host.
17. The cluster profiles requested `cpus_per_task: 8` for four rules that declare
    no `threads:` directive and therefore run single-threaded, reserving seven
    idle cores each; and `check_sge_profile.py` verified its `qsub -w v` requests
    against an assumed 8 threads rather than the number Snakemake would actually
    pass, so the request it validated was not the request that would be
    submitted. Threads are now stated per rule in `set-threads` in all three
    profiles, the per-slot memory divisor follows from that number, and
    `tests/run_tests.sh` fails if the profiles disagree with each other or if
    `cpus_per_task` contradicts `set-threads`.

Test-harness errors — recorded because they were mistaken for code bugs:

16. The first `test_hudson` reused the same two sequences in both populations,
    creating artificial self-pairs; at n = 2 the unbiased within-population π and
    the cross-population dxy are both consistent but high-variance, so `da` can
    legitimately be negative. The test, not the code, was wrong.
17. The first community test put the genome column last; the real input files put
    it first. The script was correct.

## Known limitations

These are stated so they can be addressed or disclosed, not hidden.

The animal niche is roughly 80 % mouse. Host is balanced by aggregation and by a
cap on the dominant host in resampling, and a host-resolved analysis is run
separately, but "animal gut" in this dataset largely means "mouse gut" and no
statistical device changes that.

Ortholog gene families do not span 60,000 species, so Scoary2 at the species
level operates on functional annotation matrices (KO, Pfam, CAZyme, …), not on
gene families. True ortholog resolution comes from the per-genus Panaroo
pangenomes, which necessarily cover a smaller slice of the data.

Directionality rests on within-species polymorphism. Species with few strains in
one niche give `unresolved`, and that is the correct answer, not a failure.

fastsimcoal2 and HyPhy are documented, optional, manual paths. Shipping untested
demographic models would be worse than shipping none.

The cross-feeding byproduct assignments are proxies for six of seven metabolites,
as described in Stage 7.

Nothing in `scripts/R/` and no Snakemake DAG was executed in the audit
environment. Their correctness rests on line-by-line review, on the static checks
in `tests/run_tests.sh` (every rule's script exists; every attached R package is
in `envs/r.yaml`), and on their first real run on your cluster.
