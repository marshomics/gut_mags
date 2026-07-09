#!/usr/bin/env Rscript
# fig_tree_niche.R
# ---------------------------------------------------------------------------
# Taxonomic centrepiece: the bacterial GTDB scaffold with, as concentric rings,
# each species' occupancy of the human / animal / free niches and whether it is
# an undescribed (placeholder) species. Reading the rings shows whether niches
# map onto distinct clades (phylogenetic conservatism of niche) and where novel
# diversity sits - the visual companion to the phylo-community and specificity
# tests.

suppressPackageStartupMessages({
  library(optparse); library(ape); library(ggtree); library(ggplot2)
  library(ggnewscale); library(yaml)
})
source(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "theme_pub.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"), make_option("--tree"), make_option("--tip-map"),
  make_option("--species-table"), make_option("--out")
)))

cfg <- read_yaml(opt$config)
niches <- unlist(cfg$inputs$niche_levels)
tree <- read.tree(opt$tree)
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)
sp$tip <- tipmap$tip_label[match(sp$species, tipmap$species)]
sp <- sp[!is.na(sp$tip) & sp$tip %in% tree$tip.label, ]

# occupancy rings (one column per niche) + novelty
occ <- data.frame(row.names = tree$tip.label)
for (n in niches) {
  v <- rep(0L, length(tree$tip.label)); names(v) <- tree$tip.label
  hit <- sp$tip[sp[[paste0("n_", n)]] > 0]
  v[intersect(hit, names(v))] <- 1L
  occ[[n]] <- v[rownames(occ)]
}
nov <- rep(0L, length(tree$tip.label)); names(nov) <- tree$tip.label
nov[sp$tip[sp$species_is_placeholder %in% c(TRUE, "True", "TRUE")]] <- 1L

p <- ggtree(tree, layout = "circular", linewidth = 0.12)
# niche occupancy rings (binary -> niche colour vs white)
occ_long <- occ
p <- gheatmap(p, occ_long, width = 0.18, colnames = TRUE, colnames_angle = 90,
              font.size = 1.6, color = NA) +
  scale_fill_gradient(low = "white", high = "#444444", guide = "none")
p <- p + new_scale_fill()
p <- gheatmap(p, data.frame(novel = factor(nov[rownames(occ)]), row.names = rownames(occ)),
              width = 0.05, offset = 0.22, colnames = TRUE, colnames_angle = 90,
              font.size = 1.6, color = NA) +
  scale_fill_manual(values = c("0" = "white", "1" = "#AA3377"),
                    name = "undescribed", labels = c("named", "novel"))
p <- p + ggtitle("Bacterial scaffold: niche occupancy and novelty") +
  theme(plot.title = element_text(face = "bold"))

save_pub(p, opt$out, cfg, width = 180, height = 180)
cat("fig_tree_niche.R done\n")
