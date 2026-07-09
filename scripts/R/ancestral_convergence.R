#!/usr/bin/env Rscript
# ancestral_convergence.R
# ---------------------------------------------------------------------------
# Where on the tree did each human-signature function arise? Stochastic
# character mapping (simmap) reconstructs gains and losses of each top signature
# feature over the bacterial scaffold. A feature gained ONCE in a human-
# associated clade is consistent with shared ancestry; a feature gained
# INDEPENDENTLY many times in human-associated lineages is convergent and is
# much stronger evidence of selection by the human gut environment. This turns
# a correlation ("present in human species") into an evolutionary statement
# ("repeatedly acquired in the human gut").
#
# For each of the top signature features:
#   * make.simmap (phytools) with the configured transition model, N reps;
#   * mean number of 0->1 transitions (independent gains) and 1->0 (losses);
#   * convergence flag if mean gains >= convergence_min_independent_gains.
# One representative simmap per feature is saved for the figure module.
#
# Outputs: convergence_summary.tsv, simmaps/<feature>.rds

suppressPackageStartupMessages({
  library(optparse); library(ape); library(phytools); library(arrow); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--tree"),
  make_option("--tip-map"),
  make_option("--signatures"),    # combined consensus signatures (with layer column)
  make_option("--presence"),      # combined presence parquet (all layers)
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
cc <- cfg$comparative
dir.create(file.path(opt$`out-dir`, "simmaps"), showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(opt$tree)
if (!is.ultrametric(tree)) tree <- force.ultrametric(tree, method = "extend")
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)

sig <- read.delim(opt$signatures, stringsAsFactors = FALSE)
sig <- sig[sig$consensus_signature %in% c(TRUE, "True", "TRUE") &
           sig$direction == "human_enriched", ]
sig <- sig[order(-abs(sig$consensus_log2or)), ]
feats <- head(unique(sig$feature), cc$max_features_ancestral)

pres <- as.data.frame(read_parquet(opt$presence))
pres <- pres[pres$feature %in% feats, ]
pres$tip <- tipmap$tip_label[match(pres$species, tipmap$species)]

rows <- list()
for (f in feats) {
  x <- setNames(rep(0L, length(tree$tip.label)), tree$tip.label)
  hit <- pres$tip[pres$feature == f & pres$present == 1]
  x[intersect(hit, names(x))] <- 1L
  if (sum(x) < 3 || sum(x) > length(x) - 3) next
  states <- setNames(as.character(x), names(x))
  sm <- tryCatch(make.simmap(tree, states, model = cc$ancestral_transition_model,
                             nsim = cc$ancestral_reps, message = FALSE),
                 error = function(e) NULL)
  if (is.null(sm)) next
  ct <- summary(sm)$count          # columns include "0,1" (gains) and "1,0" (losses)
  gains_col <- grep("0,1", colnames(ct))
  loss_col  <- grep("1,0", colnames(ct))
  mean_gains <- if (length(gains_col)) mean(ct[, gains_col[1]]) else NA
  mean_loss  <- if (length(loss_col)) mean(ct[, loss_col[1]]) else NA
  rows[[f]] <- data.frame(feature = f, n_species_present = sum(x),
                          mean_independent_gains = mean_gains,
                          mean_losses = mean_loss,
                          convergent = !is.na(mean_gains) &&
                            mean_gains >= cc$convergence_min_independent_gains)
  saveRDS(sm[[1]], file.path(opt$`out-dir`, "simmaps", paste0(gsub("[^A-Za-z0-9]", "_", f), ".rds")))
}
res <- do.call(rbind, rows)
write.table(res, file.path(opt$`out-dir`, "convergence_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("ancestral_convergence.R done: %d features, %d convergent\n",
            nrow(res), sum(res$convergent, na.rm = TRUE)))
