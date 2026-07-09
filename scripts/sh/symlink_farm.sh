#!/usr/bin/env bash
# symlink_farm.sh <paths.tsv> <outdir> <ext>
#
# Mash and Panaroo name each sample after its file's basename, and everything
# downstream of them -- dereplication, the tree's tip labels, the niche map --
# keys on the genome id. When annotation paths come from a manifest, the basename
# is whatever the original file was called and bears no relation to the genome id.
#
# So: given a two-column TSV of "<genome>\t<path>", build a directory of symlinks
# named <genome>.<ext> pointing at the real files, and write the list of them to
# <outdir>/list.txt. Tools then see the genome ids and nothing downstream has to
# know where the file actually lives.
#
# Fails on a missing source file rather than silently producing a dangling link,
# because a genome dropped here is a genome dropped from a diversity estimate.
set -euo pipefail

tsv=${1:?usage: symlink_farm.sh <paths.tsv> <outdir> <ext>}
outdir=${2:?}
ext=${3:?}

mkdir -p "$outdir"
: > "$outdir/list.txt"

n=0
while IFS=$'\t' read -r g p; do
  [ -z "${g:-}" ] && continue
  if [ ! -e "$p" ]; then
    echo "symlink_farm: $g: no such file: $p" >&2
    exit 2
  fi
  ln -sfn "$(realpath "$p")" "$outdir/$g.$ext"
  printf '%s\n' "$outdir/$g.$ext" >> "$outdir/list.txt"
  n=$((n + 1))
done < "$tsv"

if [ "$n" -eq 0 ]; then
  echo "symlink_farm: $tsv contained no genomes" >&2
  exit 2
fi
echo "symlink_farm: linked $n genomes into $outdir as <genome>.$ext" >&2
