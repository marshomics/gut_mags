#!/usr/bin/env python3
"""
make_sensitivity_configs.py
---------------------------
Generate derived config files for the sensitivity sweep declared under
`sensitivity:` in the base config. Each derived config changes exactly one
choice and writes to its own results_dir, so the main result can be compared
against each alternative cleanly. Writing real config files (rather than fragile
command-line overrides of nested keys) keeps every sweep fully reproducible and
self-documenting.

Generates, into config/sensitivity/:
  prev01.yaml / prev09.yaml  species.presence_threshold = 0.1 / 0.9
  hqonly.yaml                qc tightened to HQ (completeness>90, contam<5)
  nomouse.yaml               host 'Mouse' dropped from the animal niche
  nophylo.yaml               phylogenetic_control off (informational flag)
"""
import argparse
import copy
import os

import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--config", required=True)
ap.add_argument("--out-dir", default="config/sensitivity")
a = ap.parse_args()

base = yaml.safe_load(open(a.config))
os.makedirs(a.out_dir, exist_ok=True)


def dump(tag, mutate):
    c = copy.deepcopy(base)
    c["results_dir"] = f"results_{tag}"
    mutate(c)
    path = os.path.join(a.out_dir, f"{tag}.yaml")
    yaml.safe_dump(c, open(path, "w"), sort_keys=False)
    print("wrote", path)


def set_prev(v):
    return lambda c: c["species"].__setitem__("presence_threshold", v)

dump("prev01", set_prev(0.1))
dump("prev09", set_prev(0.9))

def hq(c):
    c["qc"]["min_completeness"] = c["qc"]["hq_completeness"]
    c["qc"]["max_contamination"] = c["qc"]["hq_contamination"]
dump("hqonly", hq)

def nomouse(c):
    c.setdefault("filters", {})["drop_hosts"] = ["Mouse"]
dump("nomouse", nomouse)

def nophylo(c):
    c.setdefault("flags", {})["phylogenetic_control"] = "off"
dump("nophylo", nophylo)
