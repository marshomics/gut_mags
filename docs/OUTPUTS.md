# Outputs catalogue

All paths are under `results/` (or `results_<tag>/` for sensitivity runs).

## 00_ingest
- `samples.parquet` / `samples.tsv` — QC-passed per-genome table with derived
  covariates (completeness, log10 genome size, GC, HQ flag, placeholder flag).
- `samples_qc_failures.tsv` — excluded genomes with the failing reason.
- `samples_qc_report.json` — counts before/after each filter, per niche/domain.
- `species_table.tsv` — one row per species: genomes per niche, occupancy,
  specialist/generalist, primary niche, contrast niche, taxonomy, host breadth.
- `occupancy_matrix.tsv` — species x niche binary occupancy.
- `strains_per_species.tsv` — genome count per species (confounder figure input).
- `species_units_report.json` — specialist/generalist and strain-skew summary.
- `representatives.tsv` / `.txt` — species representative genome ids.
- `annotation_paths.tsv` — genome x annotation kind, holding the resolved path to
  every per-genome file (assembly, prokka gff/faa/ffn, eggnog, kofam, and the
  dbCAN/antiSMASH/AMRFinder outputs the workflow writes). Built once from either
  `inputs.annotation_manifest` or the `{genome}` templates; every downstream rule
  reads this and nothing else, so the two ways of declaring paths become
  indistinguishable past this point.

## 01_ecology
- `composition/composition_<rank>.tsv` — taxonomic composition per niche,
  species-weighted and genome-weighted.
- `div_rarefaction_curves.tsv`, `div_hill_numbers.tsv` — diversity at matched
  effort with bootstrap CIs.
- `breadth_species_breadth.tsv`, `breadth_upset_sets.tsv`,
  `breadth_summary.json` — niche breadth and occupancy overlap (named vs
  placeholder species).

## taxonomy/  (species-by-niche stage — current focus)
- `01_ecology/div_richness_estimators.tsv` — observed richness, Chao1 (+CI), ACE,
  sample coverage per niche.
- `02_specificity/taxon_specificity_<rank>.tsv` — per taxon: species per niche,
  Levins breadth, specificity, specialist flag, permutation p/q for niche
  restriction. `specialist_summary.tsv` — FDR-significant niche-specific taxa per
  rank and niche.
- `03_enrichment/enrichment_<rank>.tsv` — per taxon x niche log2 OR, p, q,
  direction. `gtest_<rank>.tsv` — rank-level G-test + standardised residuals.
- `03_enrichment_indicator/indicator_<rank>.tsv` — IndVal.g indicator taxa per
  niche/combination with permutation p and FDR q.
- `04_overlap/overlap_nullmodel.tsv` — observed vs null cross-niche overlap (SES,
  empirical p). `overlap_observed.tsv` — raw shared-species counts.
- `05_novelty/novelty_by_niche.tsv` — novel (placeholder) species fraction,
  observed and rarefied with CI. `novelty_by_rank.tsv` — placeholder-taxa
  fraction per rank x niche.
- `06_phylo_community/phylo_community.tsv` — per niche Faith's PD + SES, NRI, NTI
  with p. `phylo_beta.tsv` — pairwise UniFrac.
- `07_beta/beta_pairwise.tsv`, `beta_multi.tsv` — Sorensen turnover/nestedness.
- `08_host/host_host_summary.tsv` — per-host genomes/species + rarefied richness;
  `host_animal_with_without_mouse.tsv`; `host_host_phylum_composition.tsv`;
  `host_mouse_dependence_by_family.tsv`.
- figures under `figures/taxonomy/`: fig1 dataset overview, fig2 diversity/
  composition, fig3 occupancy, fig_specificity, fig_taxa_signatures,
  fig_phylo_community, fig_beta_overlap, fig_novelty, fig_host_resolved,
  fig_tree_niche. Report: `report/taxonomy_report.html` / `.md`.

## transition/  (within-species niche transitions — `snakemake transition_all`)
- `selection/manifest.tsv` — qualifying species, per-niche strain counts, chosen
  outgroup. `selection/species/<id>/candidate_genomes.tsv` — genome list + roles.
- `work/<id>/derep/` — dereplicated genome list + within-niche clusters.
- `work/<id>/panaroo/` — core_gene_alignment.aln, gene_presence_absence.Rtab.
- `work/<id>/gubbins/` — recombination-filtered SNP alignment.
- `work/<id>/tree/<id>.treefile` — IQ-TREE (model-selected, supported).
- `work/<id>/ancestral/ancestral_summary.tsv` — root (ancestral) niche, Mk model,
  directed transitions, transition depths, Slatkin-Maddison p, monophyly;
  `transitions.tsv`; `ancestral_nodes.rds` (figure).
- `work/<id>/popgen/popgen_diversity.tsv` — per-niche pi, theta, Tajima's D,
  singleton fraction (equal-n, bootstrapped); `sfs_<niche>.tsv`.
- `work/<id>/demography/directionality.tsv` — source/derived call from nestedness,
  private/shared sites; `joint_sfs_<X>_<Y>.tsv` (for optional fastsimcoal2).
- `work/<id>/accessory/` — genes differentially present between niches (gains).
- `work/<id>/verdict/transition_verdict.tsv` — per-derived-niche evidence flags +
  confidence tier (strong/moderate/weak).
- `meta_directionality.tsv`, `meta_summary.json`, `meta_all_calls.tsv` —
  cross-species direction counts + binomial test.
- figures under `figures/transition/<id>/` (tree, popgen) and
  `figures/transition/meta/`; report `report/transition_report.html` / `.md`.

## 02_annot
- `chunks/` — genome id chunk files for array-style annotation.
- `dbcan/`, `antismash/`, `amrfinder/` — de novo annotation outputs + sentinels.
- `parsed/<type>/chunk_*.parquet` — per-chunk uniform long tables.
- `ko.parquet`, `eggnog.parquet`, `dbcan.parquet`, `antismash.parquet`,
  `amrfinder.parquet` — aggregated long tables (genome, layer, feature, count).

## 12_community/  (community & metabolic interaction — `community_all`)
- species-level layers written into 03_profiles: `prevalence_auxotrophy.parquet`
  (present = auxotroph), `prevalence_carbon.parquet`, `prevalence_trait.parquet`,
  `prevalence_modulep.parquet` — these also enter the differential consensus, so
  `05_diff/*/{contrast}_signatures.tsv` gains human-specific auxotrophies, carbon
  uses, traits and modules.
- `community_summary.tsv` — per niche (rarefied): mean auxotrophies/species,
  carbon breadth, cross-feeding potential + active exchanges, module collective
  coverage and complementarity, with CIs.
- `auxotrophy_by_aa.tsv` — per niche x amino acid: auxotroph fraction and
  prototroph:auxotroph ratio. `crossfeeding_by_metabolite.tsv` — producers,
  consumers, potential per metabolite. `trait_composition.tsv` — trait fractions.
- figures `figures/community/community.{png,svg}` (+ `_heatmaps`); report
  `report/community_report.html`. Curated `resources/genesets/metabolite_exchange.tsv`.

## 10_redundancy/  (functional redundancy of the human gut — `redundancy_all`)
- `redundancy_summary.tsv` — per layer (and population): n species/features, median
  occupancy, % core, % rare, Ricotta D/Q/FR/relFR.
- `occupancy_<layer>.tsv` — per function: carriers, core/intermediate/rare class.
- `accumulation_<layer>.tsv` — rarefied unique-functions-vs-species curve.
- `spread_<layer>.tsv` — per function: carrier families/genera (phylogenetic breadth).
- figure `figures/redundancy/redundancy.{png,svg}`.

## 11_population/  (Western vs non-Western — `western_all`; needs the label)
- `population_turnover.tsv` — taxonomic vs functional Sorensen (+ rarefied),
  % shared species/functions, median carrier-substitution Jaccard, per layer.
- `carrier_substitution_<layer>.tsv` — per shared function, carrier-species Jaccard
  between populations.
- the divergent functions: `05_diff/*/western_vs_nonwestern_signatures.tsv` and
  enrichment `figures/enrichment/dotplot_western_vs_nonwestern.png` (4-method).
- figure `figures/population/population.{png,svg}`.
- combined `report/redundancy_report.html`.

## 09_synthesis/  (capstone answering the six questions — `synthesis_all`)
- `catalogues/human_specific_species.tsv` — human specialists + indicator-taxon
  evidence + tier; `human_specific_genes.tsv` — signature features + adaptation
  mode; `human_specific_functions.tsv` — enriched categories + pressures;
  `catalogue_counts.json`.
- `adaptation/adaptation_mode_features.tsv` + `_summary.tsv` + `_drivers.json` —
  each signature classified as acquisition / loss / differential retention, with
  genome-architecture (PGLS) and selection counts.
- `diversification/diversification_dr.tsv`, `_summary.tsv`, `_test.json` — per-tip
  DR by niche + phylogeny-aware test.
- `ecological/ecological_pressure_<contrast>.tsv` — curated human-gut pressure
  enrichment (host-glycan/mucin, bile, oxidative stress, SCFA, fiber, starch).
- `selection/clades/<genus>/...` HyPhy per family; `selection_all.tsv`,
  `selection_summary.tsv` — BUSTED positive-selection and RELAX intensification.
- figure `figures/synthesis/master.{png,svg}`; report `report/synthesis_report.html`.
- strategy: `docs/ANALYSIS_STRATEGY.md` maps the six questions to this evidence.

## 08_enrichment/  (functional enrichment of niche signatures)
- `genesets/genesets_<layer>.tsv` — feature -> category memberships per system.
- `<layer>/<contrast>_ora.tsv` — hypergeometric over-representation (per direction,
  fold enrichment, q, overlapping features) on the tested-feature background.
- `<layer>/<contrast>_gsea.tsv` — preranked GSEA (NES, padj, size, leading edge).
- `<layer>/<contrast>_combined.tsv` — ORA + GSEA merged with an agreement flag.
- `enrichment_all.tsv`, `enrichment_top_by_contrast.tsv` — aggregated; the top
  table is the both-method-significant categories per contrast.
- figures under `figures/enrichment/` (per-contrast dot plot, category-by-contrast
  themes heatmap). Also summarised in `report/report.html`.
- one-time input: `resources/genesets/kegg_pathway.tsv` from
  `resources/fetch_kegg_genesets.py`; optional GO/Pfam-clan/CAZyme-substrate/AMR
  map files dropped into the same directory enable those systems.

## scoary/  (Scoary2 pan-GWAS — 4th method; part of `functional_all`)
- `species/<layer>/<contrast>/input/` — genes.tsv, traits.tsv, tree.nwk built from
  the species-level functional matrices + GTDB species tree.
- `species/<layer>/<contrast>/out/traits/<contrast>/results.tsv` — Scoary2 output.
- `species/<layer>/<contrast>/scoary.tsv` — parsed uniform table (log2 OR, fisher
  q, empirical p, supporting pairwise comparisons, best-pairwise p, scoary_sig);
  merged into the consensus where the contrast is shared with the trio.
- `clade/selection/manifest.tsv` + `clades/<genus>/species_genomes.tsv` — clades
  chosen for clade-stratified runs.
- `clade/clades/<genus>/panaroo/gene_presence_absence.Rtab` — per-genus pangenome.
- `clade/clades/<genus>/scoary/<contrast>.tsv` — per-genus ortholog-family results.
- `clade/clade_scoary_all.tsv`, `clade/clade_scoary_summary.tsv` — aggregated.
- figures under `figures/scoary/` (per-contrast concordance, clade summary).
  Scoary2 results also surface in the functional `report/report.html`.

## 03_profiles
- `prevalence_<layer>.parquet` — species x feature prevalence, presence call,
  mean copy number (layers: ko, pfam, cog, ec, cazyme, bgc, amr, module).
- `ko_source_concordance.tsv`, `cazyme_source_concordance.tsv` — per-genome
  Jaccard between primary and cross-check annotation sources.
- `species_traits.tsv` — per-species continuous traits (richness/load per layer
  + genome traits) for comparative methods.

## 04_phylo
- `species_tree.nwk` (bacteria), `species_tree_archaea.nwk` — pruned GTDB
  scaffolds, tips = representative accessions.
- `tip_map.tsv` — tip accession to species mapping.
- `scaffold_coverage.json` — placed vs unplaced species per domain and niche.

## 05_diff  (per `<layer>/<contrast>`)
- `<contrast>_analysis_species.tsv` — species, group, covariates, genus/family.
- `<contrast>_presence.parquet` — tested species x feature presence.
- `<contrast>_tested_features.txt` — features passing the prevalence filter.
- `<contrast>_phyloglm.tsv` — Method A results (log OR, SE, p, status).
- `<contrast>_cmh.tsv` — Method B results (MH log OR, CI, p, strata used).
- `<contrast>_resampling.tsv` — Method C results (prevalence diff, CI, sign).
- `<contrast>_signatures.tsv` / `.json` — merged per-feature calls, consensus
  flag, tier, direction, effect size; summary with positive-control recovery.
- `<contrast>_permnull.json` — observed vs null signature count, empirical FDR.
- `combined/<contrast>_signatures_all.tsv` — signatures across all layers.
- `combined/presence_all.parquet` — presence across all layers (tree/ancestral).

## 06_comparative
- `<contrast>_signal/signal_binary.tsv`, `signal_continuous.tsv` — phylogenetic
  signal of niche and features (D, lambda, K).
- `pgls_results.tsv` — PGLS coefficients for continuous functional load traits.
- `<contrast>_ancestral/convergence_summary.tsv` — independent gains per feature,
  convergence flag; `simmaps/<feature>.rds` for plotting.

## 07_ordination
- `pcoa_coords.tsv`, `pcoa_varexp.tsv` — ordination coordinates and axis
  variance.
- `permanova.tsv`, `betadisper.tsv`, `betadisper_permutest.txt` — PERMANOVA and
  dispersion check.
- `varpart.tsv`, `varpart.rds` — variation partitioning fractions.

## figures  (each as .png at 400 dpi and .svg with editable text)
- `fig1_dataset_overview` — confounders made visible (strain skew, niche
  imbalance, quality, host dominance).
- `fig2_diversity_composition` — rarefaction, Hill numbers, species- vs
  genome-weighted composition.
- `fig3_occupancy` (+ `_breadth`) — UpSet of niche occupancy; niche breadth.
- `volcano/<layer>_<contrast>` — volcano + phyloglm-vs-CMH concordance.
- `heatmap/<contrast>_signature_heatmap` — top consensus features x niche.
- `fig_function_landscape` — CAZyme/BGC/AMR load per niche with PGLS annotation.
- `fig_ordination_varpart` — functional PCoA + variation partitioning.
- `tree/<contrast>_tree_annotated` — bacterial scaffold with niche + signature
  rings.
- `tree/<contrast>_convergence` — stochastic map of the most convergent feature.

## report
- `report.html` / `report.md` — headline numbers, confounder recap, signature
  tables, permutation FDR, figure gallery.
- `sensitivity_comparison.tsv` — consensus-call overlap (Jaccard, gained/lost)
  across the four sensitivity scenarios.
