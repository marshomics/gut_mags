#!/usr/bin/env Rscript
# gsea_enrichment.R
# ---------------------------------------------------------------------------
# Threshold-free preranked GSEA (fgsea) of functional categories on the signed
# niche-association statistic across ALL tested features. Unlike the
# over-representation test, this uses the entire ranking rather than a hit
# cutoff, so it detects coordinated shifts of whole categories even when few
# individual features clear significance.
#
# Rank metric (config rank_metric): sign(consensus log2 OR) * -log10(phyloglm q),
# so categories enriched among positively-associated features get a positive NES.
# fgsea is run per category system; NES, p, BH-adjusted p, set size and the
# leading-edge features are reported.
#
# Output: enrichment_gsea_<layer>_<contrast>.tsv

suppressPackageStartupMessages({
  library(optparse); library(fgsea); library(data.table); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"), make_option("--signatures"), make_option("--genesets"),
  make_option("--layer"), make_option("--contrast"), make_option("--out")
)))

cfg <- read_yaml(opt$config)
en <- cfg$enrichment
set.seed(cfg$seed)

sig <- read.delim(opt$signatures, stringsAsFactors = FALSE)
sig$feature <- as.character(sig$feature)
# rank metric: signed effect x evidence
or <- suppressWarnings(as.numeric(sig$consensus_log2or))
q  <- suppressWarnings(as.numeric(sig$pg_q))
q[is.na(q)] <- suppressWarnings(as.numeric(sig$cmh_q))[is.na(q)]
q[is.na(q)] <- 0.5
stat <- sign(or) * -log10(pmax(q, 1e-300))
stat[is.na(stat)] <- 0
ranks <- tapply(stat, sig$feature, function(x) x[which.max(abs(x))])
ranks <- sort(ranks, decreasing = TRUE)

gs <- read.delim(opt$genesets, stringsAsFactors = FALSE)
gs$feature <- as.character(gs$feature)
gs <- gs[gs$feature %in% names(ranks), ]
name_of <- setNames(gs$category_name, gs$category_id)

res_all <- list()
for (system in unique(gs$system)) {
  sub <- gs[gs$system == system, ]
  pathways <- split(sub$feature, sub$category_id)
  pathways <- pathways[lengths(pathways) >= en$min_set_size &
                       lengths(pathways) <= en$max_set_size]
  if (length(pathways) < 1) next
  fg <- tryCatch(
    fgsea(pathways = pathways, stats = ranks, eps = 0.0,
          minSize = en$min_set_size, maxSize = en$max_set_size,
          nPermSimple = en$gsea_permutations),
    error = function(e) NULL)
  if (is.null(fg) || nrow(fg) == 0) next
  fg$system <- system
  fg$category_name <- name_of[fg$pathway]
  fg$leadingEdge <- vapply(fg$leadingEdge, function(x) paste(head(x, 25), collapse = ","), "")
  res_all[[system]] <- as.data.frame(fg)
}

if (length(res_all)) {
  out <- do.call(rbind, res_all)
  out <- out[, c("system", "pathway", "category_name", "NES", "ES", "pval",
                 "padj", "size", "leadingEdge")]
  colnames(out)[colnames(out) == "pathway"] <- "category_id"
} else {
  out <- data.frame(system = character(), category_id = character(),
                    category_name = character(), NES = numeric(), ES = numeric(),
                    pval = numeric(), padj = numeric(), size = integer(),
                    leadingEdge = character())
}
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s/%s GSEA: %d categories, %d padj<%.2f\n", opt$layer, opt$contrast,
            nrow(out), sum(out$padj < en$fdr_alpha, na.rm = TRUE), en$fdr_alpha))
