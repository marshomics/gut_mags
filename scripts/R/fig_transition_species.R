#!/usr/bin/env Rscript
# fig_transition_species.R
# ---------------------------------------------------------------------------
# Per-species figure: the rooted, recombination-masked tree with tips coloured by
# niche and internal nodes shown as ancestral-state pies (marginal Mk
# probabilities). Read together with the root state this shows which niche is
# ancestral and where the transitions into the derived niche sit (shallow =
# recent). PNG + editable-text SVG.

suppressPackageStartupMessages({ library(optparse); library(ape); library(yaml) })
source(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "theme_pub.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"), make_option("--nodes"), make_option("--species-id"),
  make_option("--out")
)))

cfg <- read_yaml(opt$config)
pal <- niche_colors(cfg)
obj <- readRDS(opt$nodes)
tr <- obj$tree; states <- obj$states; ace <- obj$ace; niches <- obj$niches
tipcol <- pal[states[tr$tip.label]]
piecol <- pal[colnames(ace)]

save_base({
  par(mar = c(1, 1, 2, 1))
  plot(tr, show.tip.label = FALSE, edge.width = 0.4, no.margin = FALSE,
       main = sprintf("%s  (root: %s)", opt$`species-id`,
                      names(which.max(ace[1, ]))))
  tiplabels(pch = 19, col = tipcol, cex = 0.5)
  nodelabels(pie = ace, piecol = piecol, cex = 0.25)
  legend("bottomleft", legend = names(pal), pch = 19, col = unlist(pal),
         bty = "n", cex = 0.7)
}, opt$out, cfg, width = 140, height = 150)
cat("fig_transition_species.R done\n")
