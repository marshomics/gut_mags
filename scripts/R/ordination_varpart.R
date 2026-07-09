#!/usr/bin/env Rscript
# ordination_varpart.R
# ---------------------------------------------------------------------------
# Multivariate view of functional repertoires, with the niche effect separated
# from phylogeny, genome size and quality.
#
#   1. Species x feature presence -> Jaccard distance -> PCoA (for the figure).
#   2. PERMANOVA (adonis2, marginal) of niche WITH covariates, so the niche
#      term is the variance it explains after completeness, size and GC.
#   3. betadisper + permutest: groups of unequal size/spread can inflate
#      PERMANOVA; this checks whether a significant niche term could be a
#      dispersion artefact rather than a location difference.
#   4. Variation partitioning (varpart) of functional distance among four
#      blocks: niche | phylogeny (patristic PCoA axes) | genome size | quality.
#      The phylogeny block is the key control: it reports how much functional
#      variation niche explains that phylogeny does NOT.
#
# Species are balanced-subsampled to max_species_per_niche so the analysis is
# not dominated by the larger niche and the distance matrix stays tractable.
#
# Outputs: pcoa_coords.tsv, permanova.tsv, betadisper.tsv, varpart.tsv

suppressPackageStartupMessages({
  library(optparse); library(arrow); library(vegan); library(ape); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--prevalence"),     # prevalence_<layer>.parquet
  make_option("--species-traits"), # species_traits.tsv (has per-species covariates)
  make_option("--tree"),
  make_option("--tip-map"),
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
ocfg <- cfg$ordination
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

sp  <- read.delim(opt$`species-traits`, stringsAsFactors = FALSE)
sp$niche <- if (isTRUE(ocfg$use_specialists_only)) sp$specialist_niche else sp$niche_primary
sp <- sp[!is.na(sp$niche) & sp$niche != "", ]

# balanced subsample per niche
cap <- ocfg$max_species_per_niche
sel <- do.call(rbind, lapply(split(sp, sp$niche), function(d) {
  if (nrow(d) > cap) d[sample(nrow(d), cap), ] else d
}))
species_keep <- sel$species

prev <- as.data.frame(read_parquet(opt$prevalence))
prev <- prev[prev$species %in% species_keep & prev$present == 1, ]
# species x feature binary matrix
feats <- sort(unique(prev$feature))
M <- matrix(0L, nrow = length(species_keep), ncol = length(feats),
            dimnames = list(species_keep, feats))
M[cbind(match(prev$species, rownames(M)), match(prev$feature, colnames(M)))] <- 1L
# drop empty species / ultra-rare features
M <- M[rowSums(M) > 0, , drop = FALSE]
M <- M[, colSums(M) >= 2, drop = FALSE]
meta <- sel[match(rownames(M), sel$species), ]

# --- distance + PCoA --------------------------------------------------------
D <- vegdist(M, method = ocfg$presence_distance)
pco <- cmdscale(D, k = 4, eig = TRUE)
eig <- pco$eig; varexp <- round(100 * eig[1:4] / sum(eig[eig > 0]), 2)
coords <- data.frame(species = rownames(M),
                     PCo1 = pco$points[, 1], PCo2 = pco$points[, 2],
                     PCo3 = pco$points[, 3], PCo4 = pco$points[, 4],
                     niche = meta$niche)
attr(coords, "varexp") <- varexp
write.table(coords, file.path(opt$`out-dir`, "pcoa_coords.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(paste(c("PCo1", "PCo2", "PCo3", "PCo4"), varexp, sep = "\t"),
           file.path(opt$`out-dir`, "pcoa_varexp.tsv"))

# --- PERMANOVA (marginal, covariate-adjusted) -------------------------------
meta$log10_genome_size <- log10(pmax(meta$genome_size_mean, 1))
ad <- adonis2(D ~ niche + completeness_mean + log10_genome_size + gc_mean,
              data = meta, by = "margin",
              permutations = ocfg$permanova_permutations, na.action = na.omit)
write.table(as.data.frame(ad), file.path(opt$`out-dir`, "permanova.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# --- betadisper -------------------------------------------------------------
if (isTRUE(ocfg$betadisper)) {
  bd <- betadisper(D, meta$niche)
  pt <- permutest(bd, permutations = 999)
  out_bd <- data.frame(group = names(bd$group.distances),
                       mean_dist_to_centroid = as.numeric(bd$group.distances))
  write.table(out_bd, file.path(opt$`out-dir`, "betadisper.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(capture.output(print(pt)),
             file.path(opt$`out-dir`, "betadisper_permutest.txt"))
}

# --- variation partitioning -------------------------------------------------
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
tree <- read.tree(opt$tree)
meta$tip <- tipmap$tip_label[match(meta$species, tipmap$species)]
on_tree <- !is.na(meta$tip) & meta$tip %in% tree$tip.label
if (sum(on_tree) > 50) {
  sub_tree <- keep.tip(tree, intersect(tree$tip.label, meta$tip[on_tree]))
  coph <- cophenetic(sub_tree)
  pp <- cmdscale(as.dist(coph), k = min(ocfg$phylo_pcoa_axes, nrow(coph) - 1))
  phylo_ax <- pp[meta$tip[on_tree], , drop = FALSE]
  Msub <- M[on_tree, , drop = FALSE]
  metasub <- meta[on_tree, ]
  X_niche <- model.matrix(~ niche, metasub)[, -1, drop = FALSE]
  X_phylo <- phylo_ax
  X_size  <- metasub[, "log10_genome_size", drop = FALSE]
  X_qual  <- metasub[, c("completeness_mean"), drop = FALSE]
  vp <- varpart(vegdist(Msub, method = ocfg$presence_distance),
                X_niche, X_phylo, X_size, X_qual)
  frac <- vp$part$indfract
  write.table(data.frame(fraction = rownames(frac), frac),
              file.path(opt$`out-dir`, "varpart.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(vp, file.path(opt$`out-dir`, "varpart.rds"))
} else {
  writeLines("Too few species on tree for variation partitioning.",
             file.path(opt$`out-dir`, "varpart.tsv"))
}
cat("ordination_varpart.R done\n")
