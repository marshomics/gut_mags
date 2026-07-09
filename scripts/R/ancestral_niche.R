#!/usr/bin/env Rscript
# ancestral_niche.R
# ---------------------------------------------------------------------------
# Polarize niche history on the rooted, recombination-masked within-species tree.
# Answers "which niche is ancestral and which is a recent acquisition" with:
#   * Mk model selection (ER vs ARD by AIC; phytools::fitMk).
#   * stochastic character mapping (make.simmap) -> posterior ancestral state at
#     the root, directed transition counts (e.g. human->animal vs animal->human),
#     and the DEPTH of each transition (transitions near the tips = recent).
#   * a Slatkin-Maddison-style permutation test of niche structure (observed
#     parsimony changes vs tip-label shuffles) so we know niche is non-random on
#     the tree before interpreting direction.
#   * per-niche monophyly / nestedness: a derived niche tends to be nested
#     (monophyletic or few shallow clusters) inside the ancestral niche.
#
# Branch lengths are in substitutions/site (no clock here); transition depths are
# RELATIVE recency within the species. Absolute split timing comes from the
# demographic model. The outgroup is used only to root, then dropped.
#
# Outputs: ancestral_summary.tsv (one row), transitions.tsv (directed counts +
# mean depth), ancestral_nodes.rds (for the figure).

suppressPackageStartupMessages({
  library(optparse); library(ape); library(phytools); library(yaml)
  library(phangorn)   # parsimony() and phyDat() for the Slatkin-Maddison test
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config"),
  make_option("--tree"),            # IQ-TREE rooted tree (Newick) incl. outgroup
  make_option("--niche-map"),       # dereplicated_genomes.tsv (genome,niche,role)
  make_option("--species-id"),
  make_option("--out-dir")
)))

cfg <- read_yaml(opt$config)
set.seed(cfg$seed)
ac <- cfg$transition$ancestral
dir.create(opt$`out-dir`, showWarnings = FALSE, recursive = TRUE)

tr <- read.tree(opt$tree)
nm <- read.delim(opt$`niche-map`, stringsAsFactors = FALSE)
og <- nm$genome[nm$role == "outgroup"]

# root on outgroup (if present), then drop it
if (length(og) >= 1 && all(og %in% tr$tip.label)) {
  tr <- root(tr, outgroup = og, resolve.root = TRUE)
  rooted_by <- "outgroup"
  tr <- drop.tip(tr, og)
} else {
  tr <- midpoint.root(tr)
  rooted_by <- "midpoint"
}
tr <- multi2di(tr)                      # binary tree for ASR
tr$edge.length[tr$edge.length <= 0] <- 1e-8

states <- setNames(nm$niche[match(tr$tip.label, nm$genome)], tr$tip.label)
states <- states[!is.na(states)]
tr <- keep.tip(tr, names(states))
niches <- sort(unique(states))
x <- setNames(factor(states, levels = niches), names(states))

# --- Mk model selection ----------------------------------------------------
fits <- lapply(ac$models, function(m) tryCatch(fitMk(tr, x, model = m), error = function(e) NULL))
names(fits) <- ac$models
fits <- fits[!sapply(fits, is.null)]
aic <- sapply(fits, AIC)
best <- names(which.min(aic))

# --- stochastic mapping ----------------------------------------------------
# Q is re-estimated by ML under the AIC-selected model (equivalent to fitMk's Q).
sm <- make.simmap(tr, x, Q = "empirical", model = best,
                  nsim = ac$simmap_reps, message = FALSE)
ss <- summary(sm)
root_node <- length(tr$tip.label) + 1
root_pp <- ss$ace[as.character(root_node), ]
root_state <- names(which.max(root_pp))

# directed transition counts (mean over simmaps)
ct <- ss$count
trans_cols <- grep(",", colnames(ct), value = TRUE)
trans_mean <- colMeans(ct[, trans_cols, drop = FALSE])

# transition depths: for one representative simmap, find branches where state
# changes and record their height above the tips (small = recent)
node_h <- max(nodeHeights(tr)) - nodeHeights(tr)   # height above tips per edge end
one <- sm[[1]]
depths <- list()
for (i in seq_len(nrow(one$edge))) {
  maps <- one$maps[[i]]
  if (length(maps) > 1) {                # a change occurred on this edge
    to_states <- names(maps)[-1]
    # height at the parent end of the edge (approx transition time)
    h <- node_h[i, 1]
    for (tos in to_states) depths[[length(depths) + 1]] <- c(to = tos, depth = h)
  }
}
dep_df <- if (length(depths))
  do.call(rbind, lapply(depths, function(z) data.frame(to = z["to"], depth = as.numeric(z["depth"]))))
  else data.frame(to = character(0), depth = numeric(0))
tree_h <- max(nodeHeights(tr))
mean_depth_into <- sapply(niches, function(n) {
  d <- dep_df$depth[dep_df$to == n]
  if (length(d)) mean(d) / tree_h else NA      # normalised 0(recent)..1(deep)
})

# --- Slatkin-Maddison permutation test of niche structure ------------------
# phyDat needs a taxa x site matrix with rownames = tip labels; build it
# explicitly so tip names survive the permutation (they were being dropped).
mk_phydat <- function(v, lev) {
  m <- matrix(as.character(v), ncol = 1, dimnames = list(names(v), NULL))
  phyDat(m, type = "USER", levels = lev)
}
states_chr <- setNames(as.character(x), names(x))
obs_changes <- parsimony(tr, mk_phydat(states_chr, niches), method = "fitch")
perm <- replicate(ac$simmap_reps, {
  xs <- setNames(sample(states_chr), names(states_chr))
  parsimony(tr, mk_phydat(xs, niches), method = "fitch")
})
# one-sided: how often does a random labelling need as FEW changes as observed?
sm_p <- (1 + sum(perm <= obs_changes)) / (length(perm) + 1)

# --- per-niche monophyly ---------------------------------------------------
mono <- sapply(niches, function(n) {
  tips <- names(states)[states == n]
  if (length(tips) < 2) return(NA)
  is.monophyletic(tr, tips)
})

# --- write -----------------------------------------------------------------
summ <- data.frame(
  species_id = opt$`species-id`, rooted_by = rooted_by,
  n_tips = length(tr$tip.label), niches = paste(niches, collapse = ","),
  mk_model = best, mk_aic = unname(aic[best]),
  root_state = root_state, root_pp = max(root_pp),
  total_transitions = sum(trans_mean),
  sm_parsimony = obs_changes, sm_p = sm_p,
  stringsAsFactors = FALSE)
for (n in niches) {
  summ[[paste0("mean_depth_into_", n)]] <- mean_depth_into[n]
  summ[[paste0("monophyletic_", n)]] <- mono[n]
}
write.table(summ, file.path(opt$`out-dir`, "ancestral_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(transition = names(trans_mean), mean_count = trans_mean),
            file.path(opt$`out-dir`, "transitions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(tree = tr, states = states, ace = ss$ace, simmap = sm[[1]],
             niches = niches), file.path(opt$`out-dir`, "ancestral_nodes.rds"))
cat(sprintf("[%s] root=%s (pp=%.2f), %s transitions, SM p=%.3f\n",
            opt$`species-id`, root_state, max(root_pp),
            round(sum(trans_mean), 1), sm_p))
