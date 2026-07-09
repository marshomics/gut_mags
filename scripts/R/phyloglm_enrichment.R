#!/usr/bin/env Rscript
# phyloglm_enrichment.R  --  METHOD A: phylogenetic logistic regression
# ---------------------------------------------------------------------------
# Tests, for each feature, whether presence is associated with the focal niche
# AFTER accounting for (i) shared ancestry among species (phylogenetic
# correlation structure from the GTDB species tree) and (ii) genome quality,
# size and GC as fixed covariates. This is the control for the single biggest
# threat to validity: closely related species share gene content by descent, so
# a naive test would mistake "trait of a human-associated clade" for "trait
# selected by the human gut".
#
#   present ~ group + completeness_mean + log10_genome_size_mean + gc_mean
#   with phylogenetic logistic regression (phylolm::phyloglm, logistic_MPLE)
#
# Run on a CHUNK of features for parallelism; combined + FDR-corrected later.
# Features with no presence variation or with complete separation are returned
# with NA estimates and a status flag rather than crashing the chunk.
#
# Output TSV: feature, estimate_log_or, se, z, p, n_species, n_present, status

suppressPackageStartupMessages({
  library(optparse); library(ape); library(phylolm); library(arrow)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--tree"),            # species tree, safe tip labels
  make_option("--tip-map"),         # TSV: tip_label <tab> species
  make_option("--analysis"),        # analysis_species.tsv (species, group, covariates)
  make_option("--presence"),        # presence parquet (species, feature, present)
  make_option("--features"),        # text file: features for THIS chunk
  make_option("--btol", type = "double", default = 20),
  make_option("--boot", type = "integer", default = 0),
  make_option("--out")
)))

tipmap  <- read.delim(opt$`tip-map`, stringsAsFactors = FALSE)         # tip_label, species
meta    <- read.delim(opt$analysis, stringsAsFactors = FALSE)
tree    <- read.tree(opt$tree)
chunk   <- readLines(opt$features)
chunk   <- chunk[nchar(chunk) > 0]

pres <- as.data.frame(read_parquet(opt$presence))
pres <- pres[pres$feature %in% chunk, ]

# align species <-> tip labels; keep only species that are on the tree
meta <- merge(meta, tipmap, by = "species")
meta <- meta[meta$tip_label %in% tree$tip.label, ]
tree <- keep.tip(tree, intersect(tree$tip.label, meta$tip_label))
rownames(meta) <- meta$tip_label
meta <- meta[tree$tip.label, ]                # order to match tree

# standardise covariates (helps phyloglm convergence; coefficients comparable)
for (cv in c("completeness_mean", "log10_genome_size_mean", "gc_mean")) {
  if (cv %in% names(meta)) meta[[cv]] <- as.numeric(scale(meta[[cv]]))
}

# wide presence matrix (species x feature), 0 for absent
pres$tip_label <- tipmap$tip_label[match(pres$species, tipmap$species)]
pres <- pres[!is.na(pres$tip_label), ]
mat <- matrix(0L, nrow = length(tree$tip.label), ncol = length(chunk),
              dimnames = list(tree$tip.label, chunk))
idx_r <- match(pres$tip_label, rownames(mat))
idx_c <- match(pres$feature, colnames(mat))
ok <- !is.na(idx_r) & !is.na(idx_c)
mat[cbind(idx_r[ok], idx_c[ok])] <- as.integer(pres$present[ok])

fit_one <- function(feat) {
  y <- mat[, feat]
  n_present <- sum(y == 1)
  res <- list(feature = feat, estimate_log_or = NA, se = NA, z = NA, p = NA,
              n_species = length(y), n_present = n_present, status = "ok")
  if (n_present < 3 || n_present > length(y) - 3) {
    res$status <- "low_variation"; return(as.data.frame(res))
  }
  d <- data.frame(y = y, group = meta$group,
                  completeness_mean = meta$completeness_mean,
                  log10_genome_size_mean = meta$log10_genome_size_mean,
                  gc_mean = meta$gc_mean)
  f <- y ~ group + completeness_mean + log10_genome_size_mean + gc_mean
  m <- tryCatch(
    phyloglm(f, data = d, phy = tree, method = "logistic_MPLE",
             btol = opt$btol, boot = opt$boot),
    error = function(e) NULL)
  if (is.null(m)) { res$status <- "fit_failed"; return(as.data.frame(res)) }
  co <- summary(m)$coefficients
  if (!"group" %in% rownames(co)) { res$status <- "no_group_term"; return(as.data.frame(res)) }
  res$estimate_log_or <- unname(co["group", "Estimate"])
  res$se <- unname(co["group", "StdErr"])
  res$z  <- unname(co["group", "z.value"])
  res$p  <- unname(co["group", "p.value"])
  as.data.frame(res)
}

out <- do.call(rbind, lapply(chunk, fit_one))
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("phyloglm: %d features, %d fitted ok\n",
            nrow(out), sum(out$status == "ok")))
