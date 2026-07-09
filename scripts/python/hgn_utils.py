"""
hgn_utils.py -- shared helpers for the human-gut-niche pipeline.

Centralises config loading, logging, deterministic seeding, and I/O so that
every script behaves identically. Importing this module is the only sanctioned
way scripts obtain parameters; nothing scientific is hard-coded downstream.
"""
from __future__ import annotations
import logging
import os
import sys
import random
import hashlib
from typing import Any

import numpy as np
import pandas as pd
import yaml


def load_config(path: str) -> dict:
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    return cfg


def get_logger(name: str = "hgn") -> logging.Logger:
    log = logging.getLogger(name)
    if not log.handlers:
        h = logging.StreamHandler(sys.stderr)
        h.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            "%Y-%m-%d %H:%M:%S"))
        log.addHandler(h)
    log.setLevel(logging.INFO)
    return log


def set_global_seed(seed: int) -> None:
    """Seed every RNG the pipeline can touch. Call at the top of each script."""
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def derive_seed(base_seed: int, *tags: Any) -> int:
    """Deterministically derive a child seed from the base seed and string tags.

    Used so that, e.g., resampling iteration 7 of layer 'cazyme' always draws
    the same genomes regardless of execution order or parallelism.
    """
    key = "|".join([str(base_seed)] + [str(t) for t in tags]).encode()
    return int(hashlib.sha256(key).hexdigest(), 16) % (2**31 - 1)


def load_annotation_paths(path: str) -> dict:
    """Read annotation_paths.tsv into {genome: {kind: path}}.

    The single point at which a script learns where a genome's annotation files
    live. resolve_annotation_paths.py produces the table from either a manifest
    or the {genome} templates, so nothing downstream needs to know which was
    used. An empty cell means "not supplied" and is returned as None rather than
    "", so a caller cannot accidentally open the empty path.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    kinds = [c for c in df.columns if c != "genome"]
    return {str(r["genome"]): {k: (r[k] or None) for k in kinds}
            for _, r in df.iterrows()}


def read_table(path: str, **kw) -> pd.DataFrame:
    if str(path).endswith(".parquet"):
        return pd.read_parquet(path)
    sep = "\t"
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, **kw)


def write_table(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if path.endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Coerce listed columns to numeric, leaving NaN for blanks."""
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c].replace("", np.nan), errors="coerce")
    return df


def is_placeholder_species(name: str) -> bool:
    """GTDB placeholder / de-novo cluster labels contain bracketed tokens,
    e.g. 's__Methanocorpusculum [species 2466]' or 's__[Genus] [species 5474]'.
    They are valid species *clusters* but are not cross-database named species;
    flag them so overlap interpretation can treat them appropriately.
    """
    return ("[" in name) or ("]" in name)


def provenance_stamp(cfg: dict, extra: dict | None = None) -> dict:
    """Minimal provenance block written beside every major output."""
    import datetime
    stamp = {
        "project": cfg.get("project_name"),
        "seed": cfg.get("seed"),
        "utc": datetime.datetime.utcnow().isoformat() + "Z",
        "python": sys.version.split()[0],
    }
    if extra:
        stamp.update(extra)
    return stamp
