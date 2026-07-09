#!/usr/bin/env Rscript
# diversification_rate.R
# ---------------------------------------------------------------------------
# Does host-association change speciation rate? The per-tip DR statistic (Jetz
# et al. 2012, the inverse equal-splits measure) is computed on the GTDB species
# scaffold and compared across niches. DR needs only a tree with branch lengths
# (no fossil calibration), so it is computable here; branch lengths are in
# substitutions/site, making this a RELATIVE comparison of diversification among
# niches, which is what the question asks.
#
# Significance is assessed two ways: a Kruskal-Wallis test across niches, and a
# phylogeny-aware permutation that shuffles niche labels WITHIN phylum (so the
# null preserves the phylogenetic clustering of niche) and recomputes the
# human-minus-free median DR difference.
#
# Outputs: diversification_dr.tsv, diversification_summary.tsv, diversification_test.json

suppressPackageStartupMessages({ library(optparse); library(ape); library(yaml); library(jsonlite) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"), make_option("--species-table"), make_option("--tree"),
  make_option("--tip-map"), make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
dv <- cfg$synthesis$diversification
niches <- unlist(cfg$inputs$niche_levels)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(opt$tree)

# --- DR statistic per tip (inverse equal-splits) ---------------------------
ntip <- length(tree$tip.label)
root <- ntip + 1
parent <- setNames(tree$edge[, 1], tree$edge[, 2])          # child -> parent
elen <- setNames(tree$edge.length, tree$edge[, 2])          # child -> branch length
dr <- numeric(ntip)
for (t in seq_len(ntip)) {
  node <- t; j <- 1; es <- 0
  while (node != root) {
    es <- es + elen[as.character(node)] * (0.5 ^ (j - 1))
    node <- parent[as.character(node)]; j <- j + 1
    if (is.na(node)) break
  }
  dr[t] <- if (es > 0) 1 / es else NA
}
names(dr) <- tree$tip.label

# --- map tips -> species -> niche, phylum ----------------------------------
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)
sp$niche <- ifelse(!is.na(sp$specialist_niche) & sp$specialist_niche != "",
                   sp$specialist_niche, sp$niche_primary)
m <- data.frame(tip = names(dr), DR = as.numeric(dr), stringsAsFactors = FALSE)
m$species <- tipmap$species[match(m$tip, tipmap$tip_label)]
m$niche <- sp$niche[match(m$species, sp$species)]
m$phylum <- sp$gtdb_phylum[match(m$species, sp$species)]
m <- m[!is.na(m$niche) & m$niche %in% niches & is.finite(m$DR), ]

write.table(m, file.path(opt$`out-dir`, "diversification_dr.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
summ <- aggregate(DR ~ niche, m, function(x) c(n = length(x), median = median(x), mean = mean(x)))
summ <- do.call(data.frame, summ)
write.table(summ, file.path(opt$`out-dir`, "diversification_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- tests ------------------------------------------------------------------
kw <- tryCatch(kruskal.test(DR ~ as.factor(niche), data = m)$p.value, error = function(e) NA)
obs_diff <- NA; perm_p <- NA
if (all(c("human", "free") %in% m$niche)) {
  med <- tapply(m$DR, m$niche, median)
  obs_diff <- med["human"] - med["free"]
  B <- dv$permutations
  null <- numeric(B)
  for (b in seq_len(B)) {
    permn <- m$niche
    for (ph in unique(m$phylum)) {
      idx <- which(m$phylum == ph)
      if (length(idx) > 1) permn[idx] <- sample(m$niche[idx])
    }
    md <- tapply(m$DR, permn, median)
    null[b] <- (if ("human" %in% names(md)) md["human"] else NA) -
               (if ("free" %in% names(md)) md["free"] else NA)
  }
  null <- null[is.finite(null)]
  perm_p <- (1 + sum(abs(null) >= abs(obs_diff))) / (length(null) + 1)
}
write_json(list(kruskal_p = kw, human_minus_free_median_DR = unname(obs_diff),
                permutation_p = perm_p, null_rank = dv$null_rank,
                permutations = dv$permutations),
           file.path(opt$`out-dir`, "diversification_test.json"),
           auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("diversification_rate.R done: KW p=%.3g, human-free DR diff=%.4f, perm p=%.3g\n",
            kw, obs_diff, perm_p))
