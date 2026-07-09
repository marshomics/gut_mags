#!/usr/bin/env Rscript
# indicator_species.R
# ---------------------------------------------------------------------------
# Indicator-taxa analysis: which taxa are diagnostic of each niche (and of niche
# combinations). Sampling units are species; the grouping is each species'
# primary niche; the "community matrix" is one-hot membership of each species in
# a taxon at the focal rank. IndVal.g (group-size corrected) is used because the
# niches contain very different species counts, so the uncorrected IndVal would
# favour the larger niche. Permutation p-values are FDR-corrected within rank.
#
# This is the formal version of "taxon X characterises the human gut": it
# combines specificity (how exclusive the taxon is to the niche) and fidelity
# (how consistently it occurs there) into one tested statistic.
#
# Output: indicator_<rank>.tsv  (taxon, niche/combination, IndVal.g stat, p, q)

suppressPackageStartupMessages({
  library(optparse); library(indicspecies); library(yaml)
  library(permute)   # how() permutation control used by multipatt
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--species-table"),
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
tcfg <- cfg$taxonomy
niches <- unlist(cfg$inputs$niche_levels)
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

sp <- read.delim(opt$`species-table`, stringsAsFactors = FALSE)
# assign each species one niche for grouping
if (tcfg$niche_assignment == "specialist") {
  sp <- sp[!is.na(sp$specialist_niche) & sp$specialist_niche != "", ]
  sp$grp <- sp$specialist_niche
} else {
  sp <- sp[!is.na(sp$niche_primary) & sp$niche_primary != "", ]
  sp$grp <- sp$niche_primary
}
sp <- sp[sp$grp %in% niches, ]

nperm <- tcfg$indicator$permutations
minsp <- tcfg$indicator$min_group_species

for (rank in tcfg$ranks) {
  if (rank == "gtdb_species") next
  lab <- sp[[rank]]
  keep <- !is.na(lab) & lab != ""
  labk <- lab[keep]; grpk <- sp$grp[keep]
  # one-hot species x taxon incidence
  taxa <- sort(unique(labk))
  taxa <- taxa[table(labk)[taxa] >= minsp]          # drop rare taxa
  if (length(taxa) < 2) next
  comm <- sapply(taxa, function(t) as.integer(labk == t))
  colnames(comm) <- taxa
  comm <- comm[, colSums(comm) > 0, drop = FALSE]

  mp <- tryCatch(
    multipatt(comm, grpk, func = tcfg$indicator$func,
              control = how(nperm = nperm)),
    error = function(e) NULL)
  if (is.null(mp)) next
  s <- mp$sign
  s$taxon <- rownames(s)
  # the indicated niche/combination = columns (one per niche) that are 1
  grp_cols <- setdiff(colnames(s), c("index", "stat", "p.value", "taxon"))
  s$niche_combination <- apply(s[, grp_cols, drop = FALSE], 1, function(r)
    paste(grp_cols[which(r == 1)], collapse = "+"))
  s$q.value <- p.adjust(s$p.value, method = "BH")
  out <- s[, c("taxon", "niche_combination", "stat", "p.value", "q.value")]
  out <- out[order(out$q.value, -out$stat), ]
  write.table(out, file.path(opt$`out-dir`, paste0("indicator_", sub("gtdb_", "", rank), ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s: %d taxa, %d indicators q<0.05\n", rank, nrow(out),
              sum(out$q.value < 0.05, na.rm = TRUE)))
}
