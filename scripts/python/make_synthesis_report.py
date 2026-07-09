#!/usr/bin/env python3
"""
make_synthesis_report.py
------------------------
The capstone report: answers each of the six manuscript questions with the
integrated evidence, the catalogues and the master figure. One artefact that
states what is human-specific and what drives it.
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, provenance_stamp


def load(p):
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def jload(p):
    try:
        return json.load(open(p))
    except Exception:
        return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    R = args.results
    S = f"{R}/09_synthesis"
    contrast = cfg["synthesis"]["focal_contrast"]

    counts = jload(f"{S}/catalogues/catalogue_counts.json")
    species = load(f"{S}/catalogues/human_specific_species.tsv")
    genes = load(f"{S}/catalogues/human_specific_genes.tsv")
    funcs = load(f"{S}/catalogues/human_specific_functions.tsv")
    drivers = jload(f"{S}/adaptation/adaptation_mode_drivers.json")
    divtest = jload(f"{S}/diversification/diversification_test.json")
    seln = load(f"{S}/selection/selection_summary.tsv")
    eco = load(f"{S}/ecological/ecological_pressure_{contrast}.tsv")
    figs = sorted(glob.glob(f"{R}/figures/synthesis/*.png"))

    def md_tbl(df, n=20):
        return df.head(n).to_markdown(index=False) if len(df) else "(none)"

    md = [f"# {cfg['project_name']} — what makes the human gut human\n",
          f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']} · contrast {contrast}_\n",
          "## Q1. Which species are human-specific?\n",
          f"{counts.get('n_human_specific_species','NA')} human-niche specialists "
          f"({counts.get('n_high_tier_species','NA')} in genera that are significant human "
          f"indicator taxa). Top of the catalogue:\n", md_tbl(species, 15),
          "\n## Q2. Which genes are human-specific?\n",
          f"{counts.get('n_human_specific_genes','NA')} consensus signature features "
          f"enriched in human across layers ({counts.get('genes_per_layer')}). Top:\n",
          md_tbl(genes, 15),
          "\n## Q3. Which functions are human-specific?\n",
          f"{counts.get('n_human_specific_functions','NA')} enriched functional categories / "
          f"pressures.\n", md_tbl(funcs, 20),
          "\n## Q4. What drives niche adaptation?\n",
          f"Adaptation-mode breakdown of the human signatures: {drivers.get('mode_counts')}. "
          f"Genome architecture (PGLS niche term): {drivers.get('genome_architecture_pgls')}. "
          f"Families under positive selection (BUSTED): "
          f"{seln.iloc[0]['n_positive_selection_busted'] if len(seln) else 'NA'}; "
          f"intensified (RELAX K>1): "
          f"{seln.iloc[0]['n_intensified_relax'] if len(seln) else 'NA'}.\n",
          "## Q5. What drives functional adaptation and/or speciation?\n",
          f"Diversification by niche: human-minus-free median DR = "
          f"{divtest.get('human_minus_free_median_DR')}, Kruskal p = {divtest.get('kruskal_p')}, "
          f"phylogeny-aware permutation p = {divtest.get('permutation_p')}. Read with the "
          f"adaptation-mode breakdown (Q4) and the variation partitioning (ordination stage).\n",
          "## Q6. Why have these species adapted to the human gut?\n",
          "Curated human-gut selective pressures tested for enrichment among the human "
          "signatures (the enriched ones are the candidate reasons):\n", md_tbl(eco, 12),
          "\n## Master figure & catalogues\n"]
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"![{os.path.basename(f)}]({rel})\n")
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write("\n".join(md))

    def htbl(df, n=20):
        return df.head(n).to_html(index=False, border=0) if len(df) else "<p>(none)</p>"
    gallery = "".join(
        f"<img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}' "
        f"style='max-width:100%;border:1px solid #ddd'>" for f in figs)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>{cfg['project_name']} synthesis</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.5;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;font-size:12px}}th,td{{padding:3px 8px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}h2{{margin-top:1.6rem}}</style></head><body>
<h1>{cfg['project_name']} — what makes the human gut human</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']} · contrast {contrast}</i></p>
<h2>Q1. Which species are human-specific?</h2>
<p>{counts.get('n_human_specific_species','NA')} human specialists; {counts.get('n_high_tier_species','NA')} in human-indicator genera.</p>{htbl(species,15)}
<h2>Q2. Which genes are human-specific?</h2>
<p>{counts.get('n_human_specific_genes','NA')} consensus signature features ({counts.get('genes_per_layer')}).</p>{htbl(genes,15)}
<h2>Q3. Which functions are human-specific?</h2>{htbl(funcs,20)}
<h2>Q4. What drives niche adaptation?</h2>
<p>Adaptation modes: {drivers.get('mode_counts')}.<br>Genome architecture (PGLS): {drivers.get('genome_architecture_pgls')}.<br>
Positive selection (BUSTED): {seln.iloc[0]['n_positive_selection_busted'] if len(seln) else 'NA'};
intensified (RELAX K&gt;1): {seln.iloc[0]['n_intensified_relax'] if len(seln) else 'NA'}.</p>
<h2>Q5. What drives functional adaptation / speciation?</h2>
<p>Human−free median DR = {divtest.get('human_minus_free_median_DR')}; Kruskal p = {divtest.get('kruskal_p')}; permutation p = {divtest.get('permutation_p')}.</p>
<h2>Q6. Why the human gut specifically?</h2>{htbl(eco,12)}
<h2>Master figure</h2>{gallery}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Synthesis report: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
