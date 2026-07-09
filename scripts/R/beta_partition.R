#!/usr/bin/env Rscript
# beta_partition.R
# ---------------------------------------------------------------------------
# Decomposes between-niche taxonomic dissimilarity into turnover (species
# replacement) and nestedness (species loss/subset) following Baselga. This
# distinguishes two very different biological stories that a single overlap
# number hides: niches that share few species because they hold DIFFERENT
# lineages (turnover) versus one niche being a depauperate SUBSET of another
# (nestedness). For "what makes the human gut human", high turnover against
# free-living and animal niches is the stronger claim.
#
# Uses species-level incidence (each species once) so strain sampling does not
# inflate any niche. Sorensen family: beta.sim = turnover, beta.sne = nestedness,
# beta.sor = total.
#
# Outputs: beta_pairwise.tsv (turnover/nestedness/total per niche pair),
#          beta_multi.tsv (multiple-site components across all niches)

suppressPackageStartupMessages({
  library(optparse); library(betapart); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--species-table"),
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
niches <- unlist(cfg$inputs$niche_levels)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)
ncols <- paste0("n_", niches)
# niche x species incidence (species occupies niche if it has >=1 genome there)
inc <- t(sapply(niches, function(n) as.integer(sp[[paste0("n_", n)]] > 0)))
colnames(inc) <- sp$species
rownames(inc) <- niches
inc <- inc[, colSums(inc) > 0, drop = FALSE]

bp <- beta.pair(inc, index.family = "sorensen")
pairs <- t(combn(niches, 2))
rows <- lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs[i, 1]; b <- pairs[i, 2]
  data.frame(niche_a = a, niche_b = b,
             turnover_sim = as.matrix(bp$beta.sim)[a, b],
             nestedness_sne = as.matrix(bp$beta.sne)[a, b],
             total_sor = as.matrix(bp$beta.sor)[a, b])
})
write.table(do.call(rbind, rows), file.path(opt$`out-dir`, "beta_pairwise.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bm <- beta.multi(inc, index.family = "sorensen")
write.table(data.frame(component = names(unlist(bm)), value = unlist(bm)),
            file.path(opt$`out-dir`, "beta_multi.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("beta_partition.R done\n")
