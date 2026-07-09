#!/usr/bin/env python3
"""
make_taxonomy_report.py
-----------------------
Assemble the species/taxonomy-by-niche report: dataset, sampling-corrected
diversity, niche specificity and indicator taxa, enrichment, cross-niche overlap
vs null, phylogenetic community structure, novelty, and the host-resolved animal
analysis. Each section pairs the headline numbers with the figure. This is the
artefact for the taxonomy stage of the manuscript.
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, provenance_stamp


def jload(p):
    try:
        return json.load(open(p))
    except Exception:
        return {}


def tload(p):
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
    TX = f"{R}/taxonomy"

    qc = jload(f"{R}/00_ingest/samples_qc_report.json")
    cov = jload(f"{R}/04_phylo/scaffold_coverage.json")
    est = tload(f"{R}/01_ecology/div_richness_estimators.tsv")
    spec = tload(f"{TX}/02_specificity/specialist_summary.tsv")
    overlap = tload(f"{TX}/04_overlap/overlap_nullmodel.tsv")
    phylo = tload(f"{TX}/06_phylo_community/phylo_community.tsv")
    beta = tload(f"{TX}/07_beta/beta_pairwise.tsv")
    nov = tload(f"{TX}/05_novelty/novelty_by_niche.tsv")
    host = tload(f"{TX}/08_host/host_animal_with_without_mouse.tsv")

    controls = [
        ("Strains per species", "Every taxonomic count is species-weighted (each species once); the genome-weighted view is shown alongside only to expose the bias."),
        ("Species per niche", "Richness compared by rarefaction and Chao1/ACE at matched effort, not by raw counts; specialist and overlap tests use null models, not absolute numbers."),
        ("Host imbalance (mouse)", "The animal niche is resolved by host; richness is rarefied to a common effort and an explicit with/without-mouse comparison is reported."),
        ("Phylogeny", "Niche specificity is tested against a label-shuffling null; phylogenetic community structure (PD/NRI/NTI) is computed against tree-based nulls."),
        ("Sampling vs novelty", "Novelty (undescribed taxa) is compared rarefied to equal effort, so it is not just a function of how deeply a niche was sequenced."),
        ("Multiple testing", "Indicator and enrichment p-values are FDR-corrected within rank; overlap and specificity use empirical permutation p-values."),
    ]

    figs = sorted(glob.glob(f"{R}/figures/taxonomy/**/*.png", recursive=True))

    def num(x):
        """Thousands-separated number, or 'NA'. f'{"NA":,}' raises ValueError."""
        return f"{x:,}" if isinstance(x, (int, float)) else "NA"

    def md_tbl(df, n=12):
        return df.head(n).to_markdown(index=False) if not df.empty else "(none)"

    md = [f"# {cfg['project_name']} — taxonomy by niche\n",
          f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}_\n",
          "## 1. Dataset & QC\n"]
    if qc:
        md += [f"- Genomes passing QC: **{num(qc.get('n_passed'))} / {num(qc.get('n_input'))}**",
               f"- Per niche: {qc.get('passed_per_niche')}\n"]
    md.append("## 2. How sampling/ancestry confounders are handled\n")
    md += [f"- **{n}** — {h}" for n, h in controls]
    md += ["\n## 3. Diversity (sampling-corrected)\n", md_tbl(est),
           "\n## 4. Niche-specific taxa (FDR-significant, by rank)\n", md_tbl(spec, 30),
           "\n## 5. Cross-niche overlap vs null\n", md_tbl(overlap),
           "\n## 6. Phylogenetic community structure\n", md_tbl(phylo),
           "\n## 7. Beta diversity (turnover vs nestedness)\n", md_tbl(beta),
           "\n## 8. Novelty per niche (rarefied)\n", md_tbl(nov),
           "\n## 9. Animal niche with vs without mouse\n", md_tbl(host),
           "\n## 10. Phylogenetic scaffold coverage\n",
           f"{ {k: v['fraction_placed'] for k, v in cov.get('by_primary_niche', {}).items()} if cov else 'NA'}\n",
           "\n## 11. Figures\n"]
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"### {os.path.basename(f)}\n\n![{os.path.basename(f)}]({rel})\n")
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write("\n".join(md))

    def html_tbl(df, n=15):
        return df.head(n).to_html(index=False, border=0) if not df.empty else "<p>(none)</p>"
    gallery = "".join(
        f"<h3>{os.path.basename(f)}</h3><img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}'"
        f" style='max-width:100%;border:1px solid #ddd'>" for f in figs)
    ctrl = "".join(f"<li><b>{n}</b> — {h}</li>" for n, h in controls)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>{cfg['project_name']} taxonomy</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.45;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;font-size:13px}}th,td{{padding:3px 9px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>{cfg['project_name']} — taxonomy by niche</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}</i></p>
<h2>1. Dataset &amp; QC</h2><p>Passed QC: <b>{num(qc.get('n_passed'))}</b> / {num(qc.get('n_input'))} · per niche {qc.get('passed_per_niche')}</p>
<h2>2. Confounder handling</h2><ul>{ctrl}</ul>
<h2>3. Diversity (sampling-corrected)</h2>{html_tbl(est)}
<h2>4. Niche-specific taxa per rank</h2>{html_tbl(spec,30)}
<h2>5. Cross-niche overlap vs null</h2>{html_tbl(overlap)}
<h2>6. Phylogenetic community structure</h2>{html_tbl(phylo)}
<h2>7. Beta diversity</h2>{html_tbl(beta)}
<h2>8. Novelty per niche</h2>{html_tbl(nov)}
<h2>9. Animal with/without mouse</h2>{html_tbl(host)}
<h2>10. Figures</h2>{gallery}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Taxonomy report: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
