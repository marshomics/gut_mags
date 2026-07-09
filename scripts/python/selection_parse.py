#!/usr/bin/env python3
"""
selection_parse.py
------------------
Parse HyPhy BUSTED and RELAX JSON outputs for one clade's families.
  BUSTED  -> p-value for episodic positive selection in the family.
  RELAX   -> K (relaxation/intensification parameter; K>1 = intensified
             selection on the human-associated branches, K<1 = relaxed) + p.

Output: selection_<clade>.tsv (clade, family, busted_p, busted_sig, relax_K,
relax_p, relax_class).
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, get_logger

log = get_logger("sel-parse")


def jget(path, *keys, default=None):
    try:
        d = json.load(open(path))
        for k in keys:
            d = d[k]
        return d
    except Exception:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--sel-dir", required=True)
    ap.add_argument("--clade", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rows = []
    for famdir in sorted(glob.glob(f"{args.sel_dir}/families/*")):
        if not os.path.isdir(famdir):
            continue
        fam = os.path.basename(famdir)
        bp = jget(f"{famdir}/BUSTED.json", "test results", "p-value")
        rk = jget(f"{famdir}/RELAX.json", "test results",
                  "relaxation or intensification parameter")
        rp = jget(f"{famdir}/RELAX.json", "test results", "p-value")
        if bp is None and rk is None:
            continue
        relax_class = ("intensified" if (rk is not None and rk > 1) else
                       "relaxed" if (rk is not None and rk < 1) else "ns")
        rows.append({"clade": args.clade, "family": fam,
                     "busted_p": bp, "busted_sig": (bp is not None and bp < 0.05),
                     "relax_K": rk, "relax_p": rp,
                     "relax_class": relax_class if (rp is not None and rp < 0.05) else "ns"})
    pd.DataFrame(rows, columns=["clade", "family", "busted_p", "busted_sig",
                                "relax_K", "relax_p", "relax_class"]).to_csv(
        args.out, sep="\t", index=False)
    log.info("%s selection: %d families parsed", args.clade, len(rows))


if __name__ == "__main__":
    main()
