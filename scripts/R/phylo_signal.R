#!/usr/bin/env Rscript
# phylo_signal.R
# ---------------------------------------------------------------------------
# Measures how phylogenetically clustered (a) niche membership and (b) each
# top signature feature are. This is the quantitative statement of the central
# caveat: if niche is strongly clustered on the tree, naive enrichment tests are
# confounded, which is exactly why the differential analysis is phylogenetically
# controlled. Reporting the signal makes the threat explicit and auditable.
#
#   binary traits  : Fritz & Purvis D (caper::phylo.d). D ~ 0 => Brownian
#                    clustering, D ~ 1 => random (overdispersed). Tests against
#                    both nulls.
#   continuous     : Pagel's lambda and Blomberg's K (phytools) for functional
#                    richness/load traits.
#
# Outputs: signal_binary.tsv, signal_continuous.tsv

suppressPackageStartupMessages({
  library(optparse); library(ape); library(caper); library(phytools); library(yaml)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--tree"),
  make_option("--tip-map"),
  make_option("--traits"),          # species_traits.tsv
  make_option("--signatures"),      # consensus signatures (any layer combined) for feature signal
  make_option("--presence"),        # combined presence parquet for signature features
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

tree <- read.tree(opt$tree)
tipmap <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)
traits <- read.delim(opt$traits, stringsAsFactors = FALSE)
traits$tip <- tipmap$tip_label[match(traits$species, tipmap$species)]
traits <- traits[!is.na(traits$tip) & traits$tip %in% tree$tip.label, ]
tree <- keep.tip(tree, intersect(tree$tip.label, traits$tip))
traits <- traits[match(tree$tip.label, traits$tip), ]

# --- binary: niche (human vs rest) ------------------------------------------
traits$is_human <- as.integer(traits$niche_primary == cfg$inputs$focal_niche)
df <- data.frame(tip = traits$tip, is_human = traits$is_human)
cd <- comparative.data(tree, df, names.col = "tip")
pd_niche <- phylo.d(cd, binvar = is_human, permut = 500)
binrows <- data.frame(trait = "niche_is_human", D = pd_niche$DEstimate,
                      p_random = pd_niche$Pval1, p_brownian = pd_niche$Pval0)

# --- binary: top signature features -----------------------------------------
if (!is.null(opt$signatures) && file.exists(opt$signatures)) {
  sig <- read.delim(opt$signatures, stringsAsFactors = FALSE)
  sig <- sig[sig$consensus_signature %in% c(TRUE, "True", "TRUE"), ]
  top <- head(sig$feature[order(-abs(sig$consensus_log2or))], 40)
  if (length(top) > 0 && file.exists(opt$presence)) {
    suppressPackageStartupMessages(library(arrow))
    pres <- as.data.frame(read_parquet(opt$presence))
    pres <- pres[pres$feature %in% top, ]
    pres$tip <- tipmap$tip_label[match(pres$species, tipmap$species)]
    for (f in top) {
      v <- rep(0L, length(tree$tip.label)); names(v) <- tree$tip.label
      hit <- pres$tip[pres$feature == f & pres$present == 1]
      v[intersect(hit, names(v))] <- 1L
      if (sum(v) < 3 || sum(v) > length(v) - 3) next
      d2 <- data.frame(tip = names(v), x = v)
      cd2 <- comparative.data(tree, d2, names.col = "tip")
      pdx <- tryCatch(phylo.d(cd2, binvar = x, permut = 200), error = function(e) NULL)
      if (!is.null(pdx))
        binrows <- rbind(binrows, data.frame(trait = paste0("feature:", f),
                          D = pdx$DEstimate, p_random = pdx$Pval1,
                          p_brownian = pdx$Pval0))
    }
  }
}
write.table(binrows, file.path(opt$`out-dir`, "signal_binary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- continuous traits ------------------------------------------------------
cont_traits <- grep("(_richness|_load)$", names(traits), value = TRUE)
crows <- list()
for (tr in cont_traits) {
  x <- traits[[tr]]; names(x) <- traits$tip
  x <- x[is.finite(x)]
  if (length(unique(x)) < 5) next
  lam <- tryCatch(phylosig(tree, x, method = "lambda", test = TRUE), error = function(e) NULL)
  kk  <- tryCatch(phylosig(tree, x, method = "K", test = TRUE, nsim = 500), error = function(e) NULL)
  crows[[tr]] <- data.frame(trait = tr,
    lambda = if (!is.null(lam)) lam$lambda else NA,
    lambda_p = if (!is.null(lam)) lam$P else NA,
    K = if (!is.null(kk)) kk$K else NA,
    K_p = if (!is.null(kk)) kk$P else NA)
}
write.table(do.call(rbind, crows), file.path(opt$`out-dir`, "signal_continuous.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("phylo_signal.R done\n")
