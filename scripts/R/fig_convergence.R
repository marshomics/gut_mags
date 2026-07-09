#!/usr/bin/env Rscript
# fig_convergence.R
# ---------------------------------------------------------------------------
# Evolutionary figure: paints the stochastic map of the single most convergent
# human-enriched signature feature onto the bacterial tree. Branch colour is the
# reconstructed presence/absence state; clusters of independent green branches
# scattered across the tree are independent gains - the visual statement that
# the function was acquired repeatedly in human-associated lineages rather than
# inherited once. The independent-gain count from ancestral_convergence.R is
# annotated in the title.

suppressPackageStartupMessages({
  library(optparse); library(phytools); library(yaml)
})
source(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "theme_pub.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--convergence"),     # convergence_summary.tsv
  make_option("--simmap-dir"),      # dir with <feature>.rds
  make_option("--out")
)))

cfg <- read_yaml(opt$config)
conv <- read.delim(opt$convergence, stringsAsFactors = FALSE)
conv <- conv[is.finite(conv$mean_independent_gains), ]
conv <- conv[order(-conv$mean_independent_gains), ]
if (nrow(conv) == 0) { writeLines("no convergent features", paste0(opt$out, ".txt")); quit() }

feat <- conv$feature[1]
rds <- file.path(opt$`simmap-dir`, paste0(gsub("[^A-Za-z0-9]", "_", feat), ".rds"))
if (!file.exists(rds)) { writeLines("simmap missing", paste0(opt$out, ".txt")); quit() }
sm <- readRDS(rds)

cols <- setNames(c("grey85", "#0072B2"), c("0", "1"))
title <- sprintf("%s: ~%.1f independent gains",
                 feat, conv$mean_independent_gains[1])

save_base(
  {
    plot(sm, colors = cols, type = "fan", fsize = 0.0001, lwd = 0.3,
         ftype = "off")
    title(main = title, cex.main = 0.8)
    add.simmap.legend(colors = cols, x = -1, y = -1, prompt = FALSE,
                      vertical = TRUE, fsize = 0.7)
  },
  opt$out, cfg, width = 170, height = 170)
cat("fig_convergence.R done\n")
