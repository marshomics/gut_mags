#!/usr/bin/env bash
# run_sensitivity.sh
# -----------------------------------------------------------------------------
# Re-runs the analysis under each setting in config.yaml's `sensitivity:` block,
# each into its own results_<tag>/ tree, then compares consensus calls against
# the main run. These turn "we found X" into "X holds when we vary the species-
# presence threshold, restrict to HQ genomes, drop the dominant host, and switch
# phylogenetic correction off". Run AFTER the main pipeline; cached annotation
# and profile outputs are reused where the override does not invalidate them.
#
# Usage:  bash scripts/sh/run_sensitivity.sh "<snakemake args>"
set -euo pipefail
SMK_ARGS="${1:---use-conda -j 8}"

# 1. generate one config file per scenario
python scripts/python/make_sensitivity_configs.py --config config/config.yaml

# 2. run each scenario to the combined-signatures stage (cheaper than full report)
for tag in prev01 prev09 hqonly nomouse nophylo; do
  echo "=== sensitivity: $tag ==="
  targets=""
  for c in $(python - <<'PY'
import yaml;print(" ".join(yaml.safe_load(open("config/config.yaml"))["stats"]["contrasts"]))
PY
); do
    for l in ko module pfam cog cazyme bgc amr; do
      targets="$targets results_${tag}/05_diff/${l}/${c}_signatures.tsv"
    done
  done
  snakemake $SMK_ARGS --configfile "config/sensitivity/${tag}.yaml" $targets
done

# 3. compare against the main run
python scripts/python/compare_sensitivity.py \
  --base results \
  --runs results_prev01 results_prev09 results_hqonly results_nomouse results_nophylo \
  --out results/report/sensitivity_comparison.tsv
echo "Sensitivity comparison -> results/report/sensitivity_comparison.tsv"
