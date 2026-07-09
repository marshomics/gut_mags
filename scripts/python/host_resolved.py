#!/usr/bin/env python3
"""
host_resolved.py
----------------
The animal niche is ~80% mouse, so "animal gut" findings could be "mouse gut"
findings. This resolves the animal niche by host species and asks how much of
the animal signal depends on mouse.

Produces:
  host_summary.tsv          per host: genomes, species, share of animal species,
                            and richness rarefied to a common genome effort (CI),
                            so heavily-sampled hosts are compared fairly.
  animal_with_without_mouse.tsv   animal species richness / novelty / phylum
                            composition with mouse included vs excluded.
  host_phylum_composition.tsv     host x phylum species-weighted composition
                            (top hosts), for the host-resolved figure.
  mouse_dependence_by_family.tsv  for each family in the animal niche, the share
                            of its animal species seen only in mouse (flags taxa
                            whose 'animal' status is really 'mouse').
"""
import argparse

import numpy as np
import pandas as pd

from hgn_utils import load_config, get_logger, set_global_seed, derive_seed

log = get_logger("host")


def rarefied_species(genomes_species, depth, n_boot, seed, tag):
    """E[distinct species] when sampling `depth` genomes, bootstrap CI."""
    arr = np.asarray(genomes_species)
    vals = []
    for b in range(n_boot):
        rng = np.random.default_rng(derive_seed(seed, "host", tag, b))
        vals.append(len(set(rng.choice(arr, size=depth, replace=False))))
    v = np.array(vals)
    return v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-boot", type=int, default=300)
    args = ap.parse_args()

    cfg = load_config(args.config)
    set_global_seed(cfg["seed"])
    hcfg = cfg["taxonomy"]["host"]
    host_col = cfg["inputs"]["columns"]["host_common_name"]
    min_g = hcfg["min_genomes_per_host"]

    import os
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    s = pd.read_parquet(args.samples)
    a = s[(s["niche"] == "animal") & (s[host_col].fillna("") != "")].copy()
    a = a.rename(columns={host_col: "host"})
    animal_species = set(a["species"].unique())

    # --- per-host summary, richness rarefied to min qualifying host depth -------
    counts = a.groupby("host").agg(n_genomes=("genome", "nunique"),
                                   n_species=("species", "nunique")).reset_index()
    qual = counts[counts["n_genomes"] >= min_g]["host"].tolist()
    depth = int(a[a["host"].isin(qual)].groupby("host").size().min()) if qual else 0
    rows = []
    for _, r in counts.iterrows():
        host = r["host"]
        rk = (np.nan, np.nan, np.nan)
        if host in qual and depth > 0:
            rk = rarefied_species(a.loc[a["host"] == host, "species"].to_numpy(),
                                  depth, args.n_boot, cfg["seed"], host)
        rows.append({"host": host, "n_genomes": int(r["n_genomes"]),
                     "n_species": int(r["n_species"]),
                     "frac_of_animal_species": round(r["n_species"] / len(animal_species), 4),
                     "rarefied_richness_mean": rk[0], "lo": rk[1], "hi": rk[2],
                     "rarefied_to_genomes": depth if host in qual else np.nan})
    pd.DataFrame(rows).sort_values("n_genomes", ascending=False).to_csv(
        f"{args.out_prefix}_host_summary.tsv", sep="\t", index=False)

    # --- animal with vs without mouse ------------------------------------------
    mouse_mask = a["host"].str.lower() == "mouse"
    no_mouse = a[~mouse_mask]
    def novelty(df):
        sp = df.groupby("species").first()
        return float(sp.index.to_series().str.contains(r"\[").mean())
    wo = {
        "animal_species_richness": [a["species"].nunique(), no_mouse["species"].nunique()],
        "n_hosts": [a["host"].nunique(), no_mouse["host"].nunique()],
        "novel_species_fraction": [round(novelty(a), 4), round(novelty(no_mouse), 4)],
        "species_only_in_mouse": [int(len(animal_species - set(no_mouse["species"]))), np.nan],
    }
    pd.DataFrame(wo, index=["with_mouse", "without_mouse"]).T.reset_index().rename(
        columns={"index": "metric"}).to_csv(
        f"{args.out_prefix}_animal_with_without_mouse.tsv", sep="\t", index=False)

    # --- host x phylum composition (top hosts, species-weighted) ---------------
    top_hosts = counts.sort_values("n_genomes", ascending=False).head(10)["host"].tolist()
    hp = (a[a["host"].isin(top_hosts)].drop_duplicates(["host", "species"])
          .groupby(["host", "gtdb_phylum"])["species"].nunique().rename("n_species").reset_index())
    hp["proportion"] = hp["n_species"] / hp.groupby("host")["n_species"].transform("sum")
    hp.to_csv(f"{args.out_prefix}_host_phylum_composition.tsv", sep="\t", index=False)

    # --- mouse dependence per family -------------------------------------------
    fam = a.drop_duplicates(["species"])[["species", "gtdb_family"]]
    sp_hosts = a.groupby("species")["host"].apply(set)
    fam = fam.assign(only_mouse=fam["species"].map(
        lambda sp: sp_hosts.get(sp, set()) == {"Mouse"}))
    md = fam.groupby("gtdb_family").agg(n_animal_species=("species", "nunique"),
                                        n_only_mouse=("only_mouse", "sum")).reset_index()
    md = md[md["n_animal_species"] >= 3]
    md["frac_only_mouse"] = (md["n_only_mouse"] / md["n_animal_species"]).round(3)
    md.sort_values("frac_only_mouse", ascending=False).to_csv(
        f"{args.out_prefix}_mouse_dependence_by_family.tsv", sep="\t", index=False)

    log.info("Host-resolved: %d hosts, mouse = %.0f%% of animal genomes; "
             "%d animal species only in mouse",
             counts["host"].nunique(), 100 * mouse_mask.mean(),
             int(len(animal_species - set(no_mouse['species']))))


if __name__ == "__main__":
    main()
