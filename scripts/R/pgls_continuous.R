#!/usr/bin/env Rscript
# pgls_continuous.R
# ---------------------------------------------------------------------------
# Phylogenetic generalised least squares for continuous functional load traits
# (CAZyme / BGC / AMR richness, genome size, CDS count) vs niche, controlling
# for genome quality and GC and for shared ancestry. Answers questions like
# "do human-gut genomes carry more CAZymes than free-living relatives, beyond
# what their phylogeny and genome size predict?" without the pseudo-replication
# that an ordinary regression over 581k genomes would suffer.
#
#   trait ~ niche + completeness + log10(genome size) + GC ,  Pagel's lambda
#   correlation estimated jointly (phylolm, model = "lambda").
#
# Output: pgls_results.tsv  (trait, term, estimate, se, t, p, lambda, n)

suppressPackageStartupMessages({
  library(optparse); library(ape); library(phylolm); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--tree"),
  make_option("--tip-map"),
  make_option("--traits"),
  make_option("--out")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
focal <- cfg$inputs$focal_niche

tree <- read.tree(opt$tree)
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
tr <- read.delim(opt$traits, stringsAsFactors = FALSE)
tr$tip <- tipmap$tip_label[match(tr$species, tipmap$species)]
tr <- tr[!is.na(tr$tip) & tr$tip %in% tree$tip.label, ]
tr <- tr[!is.na(tr$niche_primary) & tr$niche_primary != "", ]
tree <- keep.tip(tree, intersect(tree$tip.label, tr$tip))
tr <- tr[match(tree$tip.label, tr$tip), ]

tr$niche <- relevel(factor(tr$niche_primary), ref = focal)
tr$log10_genome_size <- log10(pmax(tr$genome_size_mean, 1))

want <- cfg$comparative$pgls_traits
present_traits <- intersect(c(want, paste0(sub("_count", "", want), "_richness")),
                            names(tr))
# map common names to columns produced by species_trait_table.py
alias <- c(cazyme_count = "cazyme_richness", bgc_count = "bgc_richness",
           amr_count = "amr_richness", genome_size = "genome_size_mean",
           cds_number = "cds_number_mean")
targets <- unique(unname(ifelse(want %in% names(alias), alias[want], want)))
targets <- intersect(targets, names(tr))

rows <- list()
for (tt in targets) {
  y <- tr[[tt]]
  if (length(unique(y[is.finite(y)])) < 5) next
  d <- data.frame(y = y, niche = tr$niche,
                  completeness = scale(tr$completeness_mean),
                  log10_genome_size = scale(tr$log10_genome_size),
                  gc = scale(tr$gc_mean))
  rownames(d) <- tree$tip.label
  m <- tryCatch(phylolm(y ~ niche + completeness + log10_genome_size + gc,
                        data = d, phy = tree, model = "lambda"),
                error = function(e) NULL)
  if (is.null(m)) next
  co <- summary(m)$coefficients
  lam <- if (!is.null(m$optpar)) m$optpar else NA
  for (term in rownames(co)) {
    rows[[paste(tt, term)]] <- data.frame(
      trait = tt, term = term,
      estimate = co[term, 1], se = co[term, 2],
      t = co[term, 3], p = co[term, 4],
      lambda = lam, n = length(y))
  }
}
res <- do.call(rbind, rows)
write.table(res, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("pgls_continuous.R done: %d trait-terms\n", nrow(res)))
