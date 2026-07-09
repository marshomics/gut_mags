#!/usr/bin/env Rscript
# phylo_community.R
# ---------------------------------------------------------------------------
# Phylogenetic community structure per niche on the GTDB scaffold. Answers the
# taxonomic-level version of "what makes the human gut human": is the human gut
# a phylogenetically CLUSTERED slice of bacterial diversity (few deep lineages,
# densely sampled) or a DISPERSED one (many distant lineages)?
#
#   Faith's PD + SES(PD)        : phylogenetic richness vs richness-matched null
#   MPD  -> NRI = -ses.mpd.z    : clustering at deep tree scale
#   MNTD -> NTI = -ses.mntd.z   : clustering at the tips (recent radiation)
# Null model and runs from config. Presence-based (a species occupies a niche)
# so strain sampling does not bias the community matrix. Phylogenetic beta
# diversity (UniFrac) between niches is also reported.
#
# Outputs: phylo_community.tsv (per niche), phylo_beta.tsv (pairwise UniFrac)

suppressPackageStartupMessages({
  library(optparse); library(ape); library(picante); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--species-table"),
  make_option("--tree"),
  make_option("--tip-map"),
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
pc <- cfg$taxonomy$phylo_community
niches <- unlist(cfg$inputs$niche_levels)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(opt$tree)
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)
sp$tip <- tipmap$tip_label[match(sp$species, tipmap$species)]
sp <- sp[!is.na(sp$tip) & sp$tip %in% tree$tip.label, ]
tree <- keep.tip(tree, intersect(tree$tip.label, sp$tip))

# niche x tip community matrix (presence = species occupies niche)
comm <- sapply(niches, function(n) {
  v <- rep(0L, length(tree$tip.label)); names(v) <- tree$tip.label
  tips <- sp$tip[sp[[paste0("n_", n)]] > 0]
  v[intersect(tips, names(v))] <- 1L
  v
})
comm <- t(comm); rownames(comm) <- niches
cat(sprintf("Tree tips: %d; community matrix %d niches x %d tips\n",
            length(tree$tip.label), nrow(comm), ncol(comm)))

cphylo <- cophenetic(tree)
runs <- pc$runs; nm <- pc$null_model; aw <- isTRUE(pc$abundance_weighted)

pdres  <- ses.pd(comm, tree, null.model = nm, runs = runs, include.root = TRUE)
mpd_z  <- ses.mpd(comm, cphylo, null.model = nm, abundance.weighted = aw, runs = runs)
mntd_z <- ses.mntd(comm, cphylo, null.model = nm, abundance.weighted = aw, runs = runs)

out <- data.frame(
  niche = rownames(comm),
  richness = pdres$ntaxa,
  PD = pdres$pd.obs,
  PD_ses_z = pdres$pd.obs.z, PD_p = pdres$pd.obs.p,
  MPD = mpd_z$mpd.obs, NRI = -mpd_z$mpd.obs.z, MPD_p = mpd_z$mpd.obs.p,
  MNTD = mntd_z$mntd.obs, NTI = -mntd_z$mntd.obs.z, MNTD_p = mntd_z$mntd.obs.p)
write.table(out, file.path(opt$`out-dir`, "phylo_community.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# phylogenetic beta (UniFrac) between niches
uf <- as.matrix(unifrac(comm, tree))
ufdf <- data.frame(niche = rownames(uf), uf, check.names = FALSE)
write.table(ufdf, file.path(opt$`out-dir`, "phylo_beta.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("phylo_community.R done\n")
