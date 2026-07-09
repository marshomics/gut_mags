#!/usr/bin/env python3
"""
make_community_report.py
------------------------
Report for the community / metabolic-interaction module. Two levels:
  microbial   which auxotrophies, carbon uses, traits and modules are human-
              specific (from the 4-method consensus on the community layers).
  community   niche-pool metrics (rarefied): auxotrophy dependency, byproduct
              cross-feeding, module division-of-labour, trait composition.
States the two confounder controls (species collapse + species-per-niche
rarefaction).
"""
import argparse
import glob
import os

import pandas as pd

from hgn_utils import load_config, provenance_stamp


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
    CM = f"{R}/12_community"
    focal = cfg["inputs"]["focal_niche"]
    contrast = cfg["synthesis"]["focal_contrast"]
    clayers = set(cfg["community"]["differential_layers"])

    summ = load(f"{CM}/community_summary.tsv")
    aa = load(f"{CM}/auxotrophy_by_aa.tsv")
    cf = load(f"{CM}/crossfeeding_by_metabolite.tsv")
    tr = load(f"{CM}/trait_composition.tsv")
    figs = sorted(glob.glob(f"{R}/figures/community/*.png"))

    # microbial: human-specific community-layer features from the consensus
    div = load(f"{R}/05_diff/combined/{contrast}_signatures_all.tsv")
    micro = pd.DataFrame()
    if len(div) and "consensus_signature" in div.columns:
        div = div[div["layer"].isin(clayers)]
        div = div[div["consensus_signature"].astype(str).isin(["True", "TRUE", "1"]) &
                  (div["direction"] == f"{focal}_enriched")]
        micro = div[["feature", "layer", "consensus_log2or", "n_methods_support"]] \
            .sort_values(["layer", "consensus_log2or"], ascending=[True, False]) if len(div) else div

    def tbl(df, n=25):
        return df.head(n).to_markdown(index=False) if len(df) else "(none / not run)"

    controls = ("Both confounders are handled consistently: each genome is a strain "
                "and capabilities are collapsed to species-level prevalence (present "
                "if in >= half of a species' strains); every community-pool metric is "
                "rarefied to equal species number per niche.")

    md = [f"# {cfg['project_name']} — community & metabolic interaction\n",
          f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}_\n",
          f"{controls}\n",
          "## Microbial level: human-specific metabolic features\n",
          "Auxotrophies, carbon uses, traits and modules enriched in human-gut "
          "species by the 4-method consensus (auxotrophy present = the species cannot "
          "synthesise it, i.e. depends on host/diet/community):\n", tbl(micro),
          "\n## Community level (rarefied to equal species per niche)\n",
          tbl(summ),
          "\n### Community trait composition\n", tbl(tr, 40),
          "\n### Cross-feeding potential by metabolite\n", tbl(cf, 40),
          "\n### Auxotrophy by amino acid\n", tbl(aa, 60),
          "\n## Figures\n"]
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"![{os.path.basename(f)}]({rel})\n")
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write("\n".join(md))

    def htbl(df, n=30):
        return df.head(n).to_html(index=False, border=0) if len(df) else "<p>(none / not run)</p>"
    gallery = "".join(f"<img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}' "
                      f"style='max-width:100%;border:1px solid #ddd'>" for f in figs)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>{cfg['project_name']} community</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.5;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;font-size:12px}}th,td{{padding:3px 8px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>{cfg['project_name']} — community &amp; metabolic interaction</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}</i></p>
<p>{controls}</p>
<h2>Microbial: human-specific metabolic features</h2>{htbl(micro)}
<h2>Community level (rarefied)</h2>{htbl(summ)}
<h3>Trait composition</h3>{htbl(tr,40)}
<h3>Cross-feeding by metabolite</h3>{htbl(cf,40)}
<h3>Auxotrophy by amino acid</h3>{htbl(aa,60)}
<h2>Figures</h2>{gallery}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Community report: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
