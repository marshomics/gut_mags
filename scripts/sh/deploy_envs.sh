#!/usr/bin/env bash
# deploy_envs.sh -- run this ON THE CLUSTER, which has no internet.
#
#   bash scripts/sh/deploy_envs.sh /shared/hgn/conda-envs
#   bash scripts/sh/deploy_envs.sh /shared/hgn/conda-envs --packed deploy/packed
#
# Puts the environments built by scripts/sh/build_envs.sh where Snakemake will
# look for them, then proves they work before a single job is submitted.
#
# Snakemake looks for <conda-prefix>/<hash>, where <hash> is the MD5 of the
# environment file's content. That is why the directory names built elsewhere are
# already correct here: the hash does not depend on the machine or the path.
#
# Two modes, matching the two build modes:
#
#   default   the environments were rsync'd to this exact path. Nothing to
#             relocate; we only verify.
#   --packed  the environments arrive as conda-pack tarballs and must be
#             unpacked here and have their embedded absolute paths rewritten
#             (conda-unpack). Required when the build path and this path differ.
set -euo pipefail

PREFIX=${1:?usage: deploy_envs.sh <conda-prefix> [--packed <dir>]}
MODE=${2:-}
PACKED=${3:-}
REPO=$(cd "$(dirname "$0")/../.." && pwd)

if [ "$MODE" = "--packed" ]; then
  [ -d "${PACKED:-}" ] || { echo "usage: deploy_envs.sh <prefix> --packed <dir>" >&2; exit 1; }
  echo "== verifying checksums"
  ( cd "$PACKED" && sha256sum -c SHA256SUMS )

  echo "== unpacking into $PREFIX"
  mkdir -p "$PREFIX"
  # Snakemake's bookkeeping files sit next to the env directories, not inside
  # them; conda-pack does not carry them, so they are shipped separately.
  [ -f "$PACKED/_metadata.tar.gz" ] && tar -xzf "$PACKED/_metadata.tar.gz" -C /

  for tgz in "$PACKED"/*.tar.gz; do
    name=$(basename "$tgz" .tar.gz)
    [ "$name" = "_metadata" ] && continue
    dest="$PREFIX/$name"
    echo "   -- $name"
    rm -rf "$dest"; mkdir -p "$dest"
    tar -xzf "$tgz" -C "$dest"
    # rewrite the absolute paths baked into shebangs, pkg-config, R's Makeconf...
    if [ -x "$dest/bin/conda-unpack" ]; then
      "$dest/bin/conda-unpack"
    else
      echo "      WARNING: no conda-unpack in $name; paths may still point at the build machine" >&2
    fi
  done
else
  [ -d "$PREFIX" ] || { echo "$PREFIX does not exist. rsync it from the build machine, "\
"or use --packed if the build path differed." >&2; exit 1; }
  echo "== $PREFIX present; same-path deployment, nothing to relocate"
fi

echo
echo "== pip-only packages (no internet here)"
if [ -d "$REPO/deploy/wheels" ] && [ -x "$PREFIX/launcher/bin/pip" ]; then
  "$PREFIX/launcher/bin/pip" install --no-index --find-links "$REPO/deploy/wheels" \
    snakemake-executor-plugin-sge 2>&1 | tail -1 || \
    echo "   (skipped; config/sge-generic does not need it)"
fi

echo
echo "== verifying every environment"
bash "$REPO/scripts/sh/verify_envs.sh" "$PREFIX"

echo
echo "== confirming Snakemake will reuse them rather than try to build"
# --conda-create-envs-only on an air-gapped host is the honest test: if any
# environment is missing or its hash does not match, this fails HERE instead of
# three hours into the first job.
"$PREFIX/launcher/bin/snakemake" --directory "$REPO" \
  --sdm conda --conda-create-envs-only --conda-prefix "$PREFIX" \
  --cores 1 --quiet 2>&1 | tail -3

cat <<EOF

Deployed. Add to your submit-host shell:
    export PATH="$PREFIX/launcher/bin:\$PATH"
    export CONDARC="$REPO/config/condarc"
and always pass --conda-prefix $PREFIX, e.g.

    snakemake --workflow-profile config/sge --sdm conda --conda-prefix $PREFIX -j 200

(the profiles set use-conda: true; --conda-prefix is what points at these).
EOF
