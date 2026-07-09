#!/usr/bin/env bash
# build_envs.sh -- run this ON A MACHINE WITH INTERNET, not on the cluster.
#
#   bash scripts/sh/build_envs.sh /shared/hgn/conda-envs [--pack]
#
# Creates every per-rule conda environment plus the Snakemake launcher
# environment, downloads the reference data that needs the network, and leaves a
# tree that scripts/sh/deploy_envs.sh moves onto the air-gapped cluster.
#
# THE ONE THING THAT MAKES OFFLINE DEPLOYMENT WORK
#
# Snakemake stores each environment in <conda-prefix>/<hash>, where <hash> is the
# MD5 of the environment FILE'S CONTENT -- not of the prefix path, not of the
# machine. So the directory names produced here are exactly the names the cluster
# will look for, and an environment built anywhere lands in the right place.
#
# What does NOT travel is the absolute path baked into a conda environment:
# shebangs, pkg-config files, R's Makeconf. Two ways out, and this script
# supports both.
#
#   SAME PATH (default). Choose a PREFIX that you can create on both machines,
#   build here, rsync there. Nothing is rewritten because nothing moved. The
#   simplest and least fragile option; it needs one writable path in common.
#
#   RELOCATE (--pack). conda-pack each environment into a tarball; deploy_envs.sh
#   unpacks it at the cluster's prefix and runs conda-unpack, which rewrites the
#   embedded paths. Use this when the two paths cannot be made to match.
#
# The build machine must be Linux x86_64 (lx-amd64, as `qhost` reports for every
# node). An environment built on macOS or aarch64 will not run on the cluster,
# and conda-pack will not save you.
#
# Channels: conda-forge and bioconda only, with --override-channels so a stray
# ~/.condarc cannot slip `defaults` in. See config/condarc for why.
set -euo pipefail

PREFIX=${1:?usage: build_envs.sh <conda-prefix> [--pack]}
PACK=${2:-}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
OUT="$REPO/deploy"
export CONDARC="$REPO/config/condarc"

command -v mamba >/dev/null || { echo "mamba not found; install mambaforge/miniforge" >&2; exit 1; }
command -v snakemake >/dev/null || { echo "snakemake not on PATH (bootstrap: mamba create -p ./sm -c conda-forge --override-channels snakemake=8.20)" >&2; exit 1; }

case "$(uname -s)/$(uname -m)" in
  Linux/x86_64) ;;
  *) echo "build on Linux x86_64; the cluster is lx-amd64 and conda envs are not portable across platforms" >&2; exit 1 ;;
esac

mkdir -p "$PREFIX" "$OUT"

echo "== 1. reject any environment that could resolve from defaults/anaconda"
python3 "$REPO/scripts/python/check_env_channels.py" "$REPO/envs" || exit 1

echo "== 2. Snakemake launcher environment (for the cluster submit host)"
mamba env create --quiet --override-channels -c conda-forge -c bioconda \
  -f "$REPO/envs/snakemake.yaml" -p "$PREFIX/launcher" 2>&1 | tail -2 || \
  echo "   (already exists; skipping)"

echo "== 3. wheels for the pip-only packages (no internet on the cluster)"
mkdir -p "$OUT/wheels"
"$PREFIX/launcher/bin/pip" download --dest "$OUT/wheels" \
  snakemake-executor-plugin-sge 2>&1 | tail -1 || \
  echo "   WARNING: could not pre-download snakemake-executor-plugin-sge; use config/sge-generic"

echo "== 4. per-rule environments into $PREFIX"
# --conda-create-envs-only resolves and installs every `conda:` directive the DAG
# can reach, without running a single job. Ask for every target, or rules only
# reachable from transition_all / scoary_all / synthesis_all get missed and the
# first cluster job discovers it.
for target in all functional_all transition_all scoary_all enrichment_all \
              synthesis_all redundancy_all community_all; do
  echo "   -- $target"
  snakemake --directory "$REPO" "$target" \
    --sdm conda --conda-create-envs-only \
    --conda-prefix "$PREFIX" --conda-frontend mamba \
    --cores 1 --quiet 2>&1 | tail -2
done

echo "== 5. record which environment maps to which hash directory"
snakemake --directory "$REPO" --sdm conda --conda-prefix "$PREFIX" \
  --list-conda-envs > "$OUT/conda_env_manifest.txt" 2>/dev/null || true
cat "$OUT/conda_env_manifest.txt"

echo "== 6. prove nothing came from defaults"
bash "$REPO/scripts/sh/verify_envs.sh" "$PREFIX" --channels-only

echo "== 7. reference data that needs the network"
python3 "$REPO/resources/fetch_kegg_modules.py"  --out "$REPO/resources/kegg_modules.tsv"
python3 "$REPO/resources/fetch_kegg_genesets.py" --out-dir "$REPO/resources/genesets"

if [ "$PACK" = "--pack" ]; then
  echo "== 8. conda-pack each environment (paths will be rewritten on unpack)"
  mkdir -p "$OUT/packed"
  for env in "$PREFIX"/*/; do
    [ -x "$env/bin/python" ] || [ -x "$env/bin/R" ] || [ -d "$env/conda-meta" ] || continue
    name=$(basename "$env")
    echo "   -- $name"
    "$PREFIX/launcher/bin/conda-pack" -p "$env" -o "$OUT/packed/$name.tar.gz" --force
  done
  # the bookkeeping files Snakemake writes NEXT TO each env dir travel too
  find "$PREFIX" -maxdepth 1 -type f -print0 | tar --null -czf "$OUT/packed/_metadata.tar.gz" -T - 2>/dev/null || true
  ( cd "$OUT/packed" && sha256sum ./*.tar.gz > SHA256SUMS )
  echo
  echo "Ship $OUT/ to the cluster, then:"
  echo "  bash scripts/sh/deploy_envs.sh <cluster-conda-prefix> --packed deploy/packed"
else
  echo
  echo "Ship the tree to the cluster at the IDENTICAL path:"
  echo "  rsync -a --info=progress2 $PREFIX/ <cluster>:$PREFIX/"
  echo "  rsync -a $OUT/ <cluster>:$REPO/deploy/"
  echo "Then on the cluster:  bash scripts/sh/deploy_envs.sh $PREFIX"
fi
