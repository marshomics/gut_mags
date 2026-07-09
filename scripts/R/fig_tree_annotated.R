#!/usr/bin/env Rscript
# fig_tree_annotated.R
# ---------------------------------------------------------------------------
# Centrepiece figure: the bacterial GTDB scaffold with niche and signature
# rings. A circular tree (tips unlabelled at this scale) carries an inner ring
# coloured by each species' primary niche and outer rings showing presence of
# the top human-enriched signature features. Read together, the figure shows
# whether the signature features track particular clades (shared ancestry) or
# are scattered across the tree (convergence) - the visual companion to the
# phylogenetic-signal and ancestral-state analyses.

suppressPackageStartupMessages({
  library(optparse); library(ape); library(ggtree); library(ggplot2)
  library(ggnewscale); library(arrow); library(yaml)
})
source(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "theme_pub.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"), make_option("--tree"), make_option("--tip-map"),
  make_option("--species-table"), make_option("--signatures"),
  make_option("--presence"), make_option("--out"), make_option("--n-features", default = 6)
)))

cfg <- read_yaml(opt$config)
tree <- read.tree(opt$tree)
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)

# tip -> niche
sp$tip <- tipmap$tip_label[match(sp$species, tipmap$species)]
niche_by_tip <- setNames(sp$niche_primary[match(tree$tip.label, sp$tip)], tree$tip.label)
dat <- data.frame(tip = tree$tip.label, niche = niche_by_tip)

# top human-enriched signature features
sig <- read.delim(opt$signatures, stringsAsFactors = FALSE)
sig <- sig[sig$consensus_signature %in% c(TRUE, "True", "TRUE") &
           sig$direction == "human_enriched", ]
sig <- sig[order(-abs(as.numeric(sig$consensus_log2or))), ]
feats <- head(unique(sig$feature), opt$`n-features`)

ring <- data.frame(row.names = tree$tip.label)
if (length(feats) > 0) {
  pres <- as.data.frame(read_parquet(opt$presence))
  pres <- pres[pres$feature %in% feats & pres$present == 1, ]
  pres$tip <- tipmap$tip_label[match(pres$species, tipmap$species)]
  for (f in feats) {
    v <- rep(0L, length(tree$tip.label)); names(v) <- tree$tip.label
    v[intersect(pres$tip[pres$feature == f], names(v))] <- 1L
    ring[[f]] <- v[rownames(ring)]
  }
}

p <- ggtree(tree, layout = "circular", linewidth = 0.15)
# inner niche ring
p <- gheatmap(p, data.frame(niche = dat$niche, row.names = dat$tip),
              width = 0.08, colnames = FALSE, color = NA) +
  scale_fill_manual(values = niche_colors(cfg), name = "niche", na.value = "grey90")

if (length(feats) > 0) {
  p <- p + new_scale_fill()
  p <- gheatmap(p, ring, width = 0.05 * length(feats), offset = 0.05,
                colnames = TRUE, colnames_angle = 90, colnames_offset_y = 0,
                font.size = 1.5, color = NA) +
    scale_fill_gradient(low = "white", high = "#222222", name = "feature\npresence")
}
p <- p + theme(legend.position = "right",
               plot.title = element_text(face = "bold")) +
  ggtitle("Bacterial scaffold: niche and human-signature features")

save_pub(p, opt$out, cfg, width = 180, height = 180)
cat("fig_tree_annotated.R done\n")
