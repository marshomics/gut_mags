#!/usr/bin/env python3
"""
resolve_annotation_paths.py
---------------------------
Resolve, once, the on-disk location of every per-genome annotation file, and
write them as one table that every downstream rule reads. Two ways to declare
where the files are, and this script hides the difference from everything else:

  MANIFEST MODE   inputs.annotation_manifest points at a tab-separated table with
                  a genome-id column and one column per annotation kind, each
                  holding a full path. Use this when the files do not follow a
                  regular naming scheme - different directories per study,
                  original basenames, symlink farms, whatever. Column headers are
                  matched case-insensitively against a set of aliases, so
                  'faa', 'prokka_faa' and 'proteins' all mean the same thing.

  TEMPLATE MODE   inputs.annotations holds path templates containing {genome}.
                  Use this when the layout is regular.

A manifest may be partial: any kind it does not supply falls back to the
template for that kind, so you can list the awkward files and template the rest.
The pipeline's own de-novo outputs (dbCAN, antiSMASH, AMRFinder) come from the
templates unless the manifest overrides them, since the workflow writes those.

Validation happens here, before any cluster time is spent: every QC-passed genome
must have a path for each REQUIRED kind, and `inputs.validate_paths` controls
whether those paths are checked for existence (all / sample / none). A genome
missing from the manifest is an error, not a silently skipped genome.

Output: annotation_paths.tsv  (genome + one column per kind, absolute or as given)
"""
import argparse
import os
import sys

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("annot-paths")

# canonical kind -> config key in inputs.annotations (template mode)
TEMPLATE_KEY = {
    "assembly": "assemblies_fasta",
    "prokka_gff": "prokka_gff",
    "prokka_faa": "prokka_faa",
    "prokka_ffn": "prokka_ffn",
    "eggnog": "eggnog_annotations",
    "kofam": "kofam_annotations",
    "dbcan_overview": "dbcan_overview",
    "antismash_json": "antismash_json",
    "amrfinder_tsv": "amrfinder_tsv",
}

# accepted manifest column headers (lowercased) -> canonical kind
ALIASES = {
    "assembly": "assembly", "assemblies_fasta": "assembly", "assembly_fasta": "assembly",
    "fna": "assembly", "fasta": "assembly", "genome_fasta": "assembly",
    "prokka_gff": "prokka_gff", "gff": "prokka_gff",
    "prokka_faa": "prokka_faa", "faa": "prokka_faa", "proteins": "prokka_faa",
    "prokka_ffn": "prokka_ffn", "ffn": "prokka_ffn", "cds": "prokka_ffn",
    "eggnog": "eggnog", "eggnog_annotations": "eggnog", "emapper": "eggnog",
    "emapper_annotations": "eggnog",
    "kofam": "kofam", "kofam_annotations": "kofam", "kofamscan": "kofam",
    "dbcan_overview": "dbcan_overview", "dbcan": "dbcan_overview",
    "antismash_json": "antismash_json", "antismash": "antismash_json",
    "amrfinder_tsv": "amrfinder_tsv", "amrfinder": "amrfinder_tsv",
}
GENOME_ALIASES = {"genome", "genome_id", "id", "file", "sample", "assembly_id"}

# kinds the pipeline reads from disk (as opposed to writing itself)
PROVIDED = ["assembly", "prokka_gff", "prokka_faa", "prokka_ffn", "eggnog", "kofam"]
GENERATED = ["dbcan_overview", "antismash_json", "amrfinder_tsv"]
ALL_KINDS = PROVIDED + GENERATED


def read_manifest(path):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        raise SystemExit(f"annotation manifest {path} is empty")
    gcol, mapping = None, {}
    for c in df.columns:
        key = c.strip().lower()
        if key in GENOME_ALIASES and gcol is None:
            gcol = c
        elif key in ALIASES:
            mapping[c] = ALIASES[key]
    if gcol is None:
        gcol = df.columns[0]
        log.warning("No recognised genome-id column in %s; using the first column %r.",
                    path, gcol)
    unknown = [c for c in df.columns if c != gcol and c not in mapping]
    if unknown:
        log.warning("Ignoring unrecognised manifest column(s): %s", ", ".join(unknown))
    if not mapping:
        raise SystemExit(
            f"annotation manifest {path} has no recognised annotation column. "
            f"Expected one or more of: {', '.join(sorted(set(ALIASES.values())))}")
    out = df[[gcol] + list(mapping)].rename(columns={gcol: "genome", **mapping})
    out["genome"] = out["genome"].astype(str).str.strip()
    dup = out["genome"].duplicated()
    if dup.any():
        raise SystemExit(f"annotation manifest {path} has {int(dup.sum())} duplicate "
                         f"genome ids, e.g. {out.loc[dup, 'genome'].iloc[0]!r}")
    log.info("Manifest %s: %d genomes, kinds supplied: %s",
             path, len(out), ", ".join(sorted(mapping.values())))
    return out.set_index("genome")


def check_existence(paths, mode, n_sample, seed):
    """paths: DataFrame of kind columns. Returns list of (genome, kind, path)."""
    if mode in (None, False, "none"):
        return []
    sub = paths
    if mode == "sample":
        n = min(n_sample, len(paths))
        sub = paths.sample(n=n, random_state=seed) if n else paths
    missing = []
    for kind in [c for c in PROVIDED if c in sub.columns]:
        for g, p in sub[kind].items():
            if pd.isna(p) or not str(p):
                continue                       # unresolved: reported elsewhere
            if not os.path.exists(str(p)):
                missing.append((g, kind, p))
    return missing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True, help="samples.parquet (QC-passed genomes)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    inp = cfg["inputs"]
    ann = inp.get("annotations", {}) or {}
    required = list(inp.get("required_annotations",
                            ["prokka_faa", "eggnog", "kofam"]))

    genomes = pd.read_parquet(args.samples)["genome"].astype(str).tolist()
    log.info("%d QC-passed genomes to resolve", len(genomes))

    man = None
    mpath = inp.get("annotation_manifest")
    if mpath and str(mpath).strip() and not str(mpath).startswith("CHANGE_ME"):
        man = read_manifest(mpath)
        missing_ids = [g for g in genomes if g not in man.index]
        if missing_ids:
            raise SystemExit(
                f"{len(missing_ids)} QC-passed genomes are absent from the annotation "
                f"manifest {mpath} (e.g. {missing_ids[:5]}). Add them, or remove them "
                f"from the metadata, or drop them with the qc thresholds. Silently "
                f"skipping genomes would bias every prevalence estimate.")

    out = pd.DataFrame({"genome": genomes}).set_index("genome")
    for kind in ALL_KINDS:
        col = None
        if man is not None and kind in man.columns:
            col = man.reindex(out.index)[kind].astype(str).str.strip().replace("", pd.NA)
        tmpl = ann.get(TEMPLATE_KEY[kind])
        if tmpl and not str(tmpl).startswith("CHANGE_ME"):
            filled = pd.Series([str(tmpl).format(genome=g) for g in out.index],
                               index=out.index)
            col = filled if col is None else col.fillna(filled)
        elif kind == "prokka_ffn" and col is None:
            faa = ann.get("prokka_faa")
            if faa and not str(faa).startswith("CHANGE_ME"):
                col = pd.Series([str(faa).format(genome=g).replace(".faa", ".ffn")
                                 for g in out.index], index=out.index)
        out[kind] = col if col is not None else pd.NA

    for kind in required:
        if kind not in out.columns or out[kind].isna().any():
            n = int(out[kind].isna().sum()) if kind in out.columns else len(out)
            raise SystemExit(
                f"annotation kind {kind!r} is required but unresolved for {n} genomes. "
                f"Supply it as a column in inputs.annotation_manifest, or as a "
                f"{{genome}} template in inputs.annotations.{TEMPLATE_KEY[kind]}.")

    mode = inp.get("validate_paths", "sample")
    missing = check_existence(out, mode, int(inp.get("validate_paths_n", 200)),
                              cfg.get("seed", 0))
    if missing:
        head = "\n  ".join(f"{g}\t{k}\t{p}" for g, k, p in missing[:10])
        raise SystemExit(
            f"{len(missing)} declared annotation file(s) do not exist "
            f"(validate_paths: {mode}). First few:\n  {head}\n"
            f"Fix the paths, or set inputs.validate_paths: none to skip this check.")
    log.info("Path existence check (%s): OK", mode)

    out.reset_index().fillna("").to_csv(args.out, sep="\t", index=False)
    log.info("Wrote %s: %d genomes x %d kinds (source: %s)", args.out, len(out),
             len(ALL_KINDS), "manifest+templates" if man is not None else "templates")


if __name__ == "__main__":
    main()
