#!/usr/bin/env python3
"""
make_redundancy_report.py
-------------------------
Report for functional redundancy and the Western vs non-Western comparison.
States how redundant the human gut's functions are, and whether the functional
repertoire is conserved across populations despite species turnover, with the
divergent functions (from the western_vs_nonwestern consensus, when enabled).
"""
import argparse
import glob
import os

import pandas as pd

from hgn_utils import load_config, get_logger, provenance_stamp

log = get_logger("redun-report")


def load(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    R = args.results
    red = load(f"{R}/10_redundancy/redundancy_summary.tsv")
    turn = load(f"{R}/11_population/population_turnover.tsv")
    figs = sorted(glob.glob(f"{R}/figures/redundancy/*.png") +
                  glob.glob(f"{R}/figures/population/*.png"))

    # divergent functions western vs non-western (consensus across layers), if present
    popc = cfg.get("population", {}).get("contrast", "western_vs_nonwestern")
    div = load(f"{R}/05_diff/combined/{popc}_signatures_all.tsv")
    if len(div) and "consensus_signature" in div.columns:
        div = div[div["consensus_signature"].astype(str).isin(["True", "TRUE", "1"])]

    def tbl(df, n=20):
        return df.head(n).to_markdown(index=False) if len(df) else "(none / not run)"

    md = [f"# {cfg['project_name']} — functional redundancy & Western vs non-Western\n",
          f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}_\n",
          "## Functional redundancy of the human gut\n",
          "Species is the unit; high relative FR (1 - Q/D) and high core fraction mean "
          "functions are carried by many species (robust); accumulation of functions "
          "saturating while species keep accumulating is the same signal.\n",
          tbl(red),
          "\n## Western vs non-Western: turnover\n",
          "Taxonomic vs functional Sorensen dissimilarity, and per-shared-function "
          "carrier-species overlap (low = same function, different species). A large "
          "taxonomic-minus-functional gap means functions are conserved despite species "
          "turnover.\n", tbl(turn),
          "\n## Western vs non-Western: divergent functions\n",
          f"Consensus {popc} signature features (the functions that DO differ):\n",
          tbl(div[["feature", "layer", "direction", "consensus_log2or", "n_methods_support"]]
              if len(div) else div, 25),
          "\n## Figures\n"]
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"![{os.path.basename(f)}]({rel})\n")
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write("\n".join(md))

    def htbl(df, n=25):
        return df.head(n).to_html(index=False, border=0) if len(df) else "<p>(none / not run)</p>"
    gallery = "".join(f"<img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}' "
                      f"style='max-width:100%;border:1px solid #ddd'>" for f in figs)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>{cfg['project_name']} redundancy</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.5;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;font-size:12px}}th,td{{padding:3px 8px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>{cfg['project_name']} — functional redundancy &amp; Western vs non-Western</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}</i></p>
<h2>Functional redundancy of the human gut</h2>{htbl(red)}
<h2>Western vs non-Western: turnover</h2>{htbl(turn)}
<h2>Western vs non-Western: divergent functions</h2>
{htbl(div[['feature','layer','direction','consensus_log2or','n_methods_support']] if len(div) else div,25)}
<h2>Figures</h2>{gallery}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Redundancy report: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
