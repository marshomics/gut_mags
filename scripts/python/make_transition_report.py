#!/usr/bin/env python3
"""
make_transition_report.py
-------------------------
Assemble the within-species niche-transition report: which species were tested,
the cross-species directional result (the headline), the per-species verdicts and
the supporting evidence, plus the figure gallery. States the confounder controls
that make the "recent acquisition" calls defensible.
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, provenance_stamp


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
    TR = f"{R}/transition"

    manifest = tload(f"{TR}/selection/manifest.tsv")
    meta_dir = tload(f"{TR}/meta_directionality.tsv")
    try:
        meta_sum = json.load(open(f"{TR}/meta_summary.json"))
    except Exception:
        meta_sum = {}
    verdicts = [tload(f) for f in glob.glob(f"{TR}/work/*/verdict/transition_verdict.tsv")]
    verdicts = [v for v in verdicts if not v.empty]
    allv = pd.concat(verdicts, ignore_index=True) if verdicts else pd.DataFrame()
    figs = sorted(glob.glob(f"{R}/figures/transition/**/*.png", recursive=True))

    controls = [
        ("Sampling per niche", "Only species with enough near-complete strains in each of >=2 niches are tested; diversity and SFS are computed on equal-n subsamples with bootstraps."),
        ("Clonal oversampling", "Strains are dereplicated within each niche (Mash) so an outbreak or over-sequenced lineage cannot fake low diversity; cross-niche near-identical strains are kept (recent-transfer signal)."),
        ("Recombination", "Core alignments are masked for recombination (Gubbins) before the tree, the diversity statistics and the SFS."),
        ("Rooting", "Trees are rooted with a congeneric outgroup; directionality is never asserted from an unrooted tree."),
        ("Convergent evidence", "A recent-acquisition call requires phylogenetic structure plus agreement among transition depth, diversity reduction, Tajima's D, nestedness and gene gain - not any single statistic."),
        ("Cross-species replication", "The headline is the consistency of direction across many independent species (binomial test), not a single example."),
    ]

    md = [f"# {cfg['project_name']} — within-species niche transitions\n",
          f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}_\n",
          "## 1. Species tested\n",
          (manifest.to_markdown(index=False) if not manifest.empty else "(none qualified)"),
          "\n## 2. Cross-species directional result\n",
          (f"Species analysed: {meta_sum.get('n_species_analysed')}; "
           f"supported acquisitions: {meta_sum.get('n_supported_acquisitions')}; "
           f"direction counts: {meta_sum.get('direction_counts')}\n"),
          (meta_dir.to_markdown(index=False) if not meta_dir.empty else "(no directional calls)"),
          "\n## 3. Confounder controls\n"]
    md += [f"- **{n}** — {h}" for n, h in controls]
    md += ["\n## 4. Per-species verdicts\n",
           (allv.to_markdown(index=False) if not allv.empty else "(none)"),
           "\n## 5. Figures\n"]
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"### {os.path.basename(f)}\n\n![{os.path.basename(f)}]({rel})\n")
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write("\n".join(md))

    def htbl(df, n=40):
        return df.head(n).to_html(index=False, border=0) if not df.empty else "<p>(none)</p>"
    gallery = "".join(
        f"<h3>{os.path.basename(f)}</h3><img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}'"
        f" style='max-width:100%;border:1px solid #ddd'>" for f in figs)
    ctrl = "".join(f"<li><b>{n}</b> — {h}</li>" for n, h in controls)
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>{cfg['project_name']} transitions</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;line-height:1.45;color:#222}}
table{{border-collapse:collapse;margin:1rem 0;font-size:12px}}th,td{{padding:3px 8px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}</style></head><body>
<h1>{cfg['project_name']} — within-species niche transitions</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}</i></p>
<h2>1. Species tested</h2>{htbl(manifest)}
<h2>2. Cross-species directional result</h2>
<p>Species analysed: {meta_sum.get('n_species_analysed')} · supported acquisitions: {meta_sum.get('n_supported_acquisitions')} · directions: {meta_sum.get('direction_counts')}</p>
{htbl(meta_dir)}
<h2>3. Confounder controls</h2><ul>{ctrl}</ul>
<h2>4. Per-species verdicts</h2>{htbl(allv)}
<h2>5. Figures</h2>{gallery}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Transition report: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
