#!/usr/bin/env bash
# verify_envs.sh <conda-prefix> [--channels-only]
#
# Two questions, both worth asking before a job is submitted rather than after.
#
#   1. Did any package come from `defaults`/`anaconda`? `conda list --explicit`
#      prints the URL every package was fetched from, so this is checkable rather
#      than assumed. An environment that names `nodefaults` in its YAML can still
#      have been built on a machine whose ~/.condarc overrode the channel list.
#
#   2. Does each environment actually contain the tools the rules invoke? A conda
#      solve can succeed and still leave you without `run_dbcan` on PATH, because
#      the package was renamed or the entry point moved. Every binary below is one
#      a rule calls by name.
set -uo pipefail

PREFIX=${1:?usage: verify_envs.sh <conda-prefix> [--channels-only]}
ONLY=${2:-}
fail=0

banned_re='repo\.anaconda\.com|repo\.continuum\.io|/pkgs/(main|r|free)/'

echo "== channel provenance"
shopt -s nullglob
for env in "$PREFIX"/*/; do
  name=$(basename "$env")
  [ -d "$env/conda-meta" ] || continue
  # every line of an explicit listing is a package URL
  if urls=$(conda list -p "$env" --explicit 2>/dev/null); then
    if bad=$(printf '%s\n' "$urls" | grep -E "$banned_re" | head -5); then
      if [ -n "$bad" ]; then
        echo "  FAIL $name: packages from a banned channel:"
        printf '        %s\n' $bad
        fail=1
      else
        echo "  ok   $name"
      fi
    fi
  else
    echo "  ??   $name: conda could not list it"
  fi
done

[ "$ONLY" = "--channels-only" ] && { [ "$fail" -eq 0 ] && echo "no defaults/anaconda packages"; exit "$fail"; }

echo
echo "== tools each rule invokes by name"
# env-file basename -> binaries that must exist. The env directory is found by
# matching the hash directory that contains the binary, since Snakemake names
# directories by hash rather than by the env's `name:` field.
declare -A NEED=(
  [annotation]="run_dbcan amrfinder antismash hmmscan diamond prodigal seqkit"
  [phylo]="gtdbtk FastTree iqtree2 mash"
  [python]="python"
  [r]="R Rscript"
  [scoary]="scoary"
  [selection]="hyphy mafft FastTree"
  [transition]="panaroo run_gubbins.py iqtree2 mafft"
  [launcher]="snakemake"
)

found_any=0
for envfile in "${!NEED[@]}"; do
  hit=""
  for env in "$PREFIX"/*/; do
    ok=1
    for bin in ${NEED[$envfile]}; do
      [ -x "$env/bin/$bin" ] || { ok=0; break; }
    done
    [ "$ok" -eq 1 ] && { hit=$env; break; }
  done
  if [ -n "$hit" ]; then
    echo "  ok   $envfile -> $(basename "$hit")"
    found_any=1
  else
    echo "  FAIL $envfile: no environment under $PREFIX has all of: ${NEED[$envfile]}"
    fail=1
  fi
done

echo
echo "== R packages every script attaches"
rbin=$(ls -d "$PREFIX"/*/bin/Rscript 2>/dev/null | head -1 || true)
if [ -n "$rbin" ]; then
  "$rbin" --vanilla -e '
    pkgs <- c("yaml","optparse","ape","phylolm","phytools","phangorn","permute",
              "jsonlite","indicspecies","betapart","vegan","picante","caper",
              "arrow","ggplot2","ggtree","svglite","fgsea","treeio","corHMM","castor")
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) { cat("  FAIL missing R packages:", paste(miss, collapse=", "), "\n"); quit(status=1) }
    cat("  ok   all", length(pkgs), "R packages load\n")' || fail=1
else
  echo "  FAIL no Rscript found under $PREFIX"; fail=1
fi

echo
if [ "$fail" -ne 0 ]; then echo "VERIFY FAILED"; exit 1; fi
echo "VERIFY OK"
