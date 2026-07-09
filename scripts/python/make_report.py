#!/usr/bin/env python3
"""
make_report.py
--------------
Assemble the run report from the JSON/TSV summaries and the figure gallery. The
report states the headline numbers, recaps how each confounder was controlled,
tabulates consensus signatures per layer x contrast, reports the permutation
empirical FDR and positive-control recovery, and embeds the figures. It is the
single artefact to read after a run.
"""
import argparse
import glob
import json
import os

import pandas as pd

from hgn_utils import load_config, provenance_stamp


def jload(path):
    try:
        return json.load(open(path))
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
    contrasts = cfg["stats"]["contrasts"]
    layers = ["ko", "module", "pfam", "cog", "cazyme", "bgc", "amr"]

    qc = jload(f"{R}/00_ingest/samples_qc_report.json")
    cov = jload(f"{R}/04_phylo/scaffold_coverage.json")

    # signature + permutation tables
    sig_rows, perm_rows = [], []
    for c in contrasts:
        for l in layers:
            s = jload(f"{R}/05_diff/{l}/{c}_signatures.json")
            if s:
                sig_rows.append({"contrast": c, "layer": l,
                                 "tested": s.get("n_tested"),
                                 "consensus": s.get("n_consensus"),
                                 "enriched": s.get("n_enriched", s.get("n_human_enriched")),
                                 "depleted": s.get("n_depleted", s.get("n_human_depleted")),
                                 "methods": ",".join(s.get("methods_used", [])),
                                 "scoary_hits": s.get("n_scoary_sig"),
                                 "pos_ctrl_missing": ", ".join(s.get("positive_controls_missing", []))})
            p = jload(f"{R}/05_diff/{l}/{c}_permnull.json")
            if p and p.get("observed_n_sig") is not None:
                perm_rows.append({"contrast": c, "layer": l,
                                  "observed": p.get("observed_n_sig"),
                                  "null_mean": round(p.get("null_mean", float("nan")), 2),
                                  "empirical_fdr": (round(p["empirical_fdr"], 4)
                                                    if p.get("empirical_fdr") is not None else None)})
    sig_df = pd.DataFrame(sig_rows)
    perm_df = pd.DataFrame(perm_rows)

    try:
        clade_sum = pd.read_csv(f"{R}/scoary/clade/clade_scoary_summary.tsv", sep="\t")
    except Exception:
        clade_sum = pd.DataFrame()
    try:
        enrich_top = pd.read_csv(f"{R}/08_enrichment/enrichment_top_by_contrast.tsv", sep="\t")
    except Exception:
        enrich_top = pd.DataFrame()

    figs = sorted(glob.glob(f"{R}/figures/**/*.png", recursive=True))

    def num(x):
        """Thousands-separated number, or 'NA'. f'{"NA":,}' raises ValueError."""
        return f"{x:,}" if isinstance(x, (int, float)) else "NA"

    confounders = [
        ("Strains per species (max 9,606; 65% singletons)",
         "Species is the unit of analysis; genomes collapsed to species-level prevalence profiles. Balanced resampling draws one genome-equivalent per species."),
        ("Species per niche (human 6k vs free 46k)",
         "Per-niche prevalence reported as proportions; coverage-based rarefaction & Hill numbers compare diversity at matched effort; balanced resampling equalises species per niche each iteration."),
        ("Host imbalance (animal ~80% mouse)",
         "Animal species drawn with inverse-host-frequency weights; host as random effect where genome-level; results reported with and without mouse in the sensitivity sweep."),
        ("Phylogenetic non-independence",
         "phyloglm models phylogenetic covariance; CMH matches within genus; phylogenetic signal (D, lambda, K) quantified; ancestral-state reconstruction tests convergence."),
        ("Genome quality / completeness",
         "Completeness and contamination entered as covariates in every model; HQ-only subset run as sensitivity; species-level prevalence buffers single-genome incompleteness."),
        ("Genome size",
         "log10 genome size entered as a covariate in every model and partialled in variation partitioning."),
        ("Multiple testing",
         "Benjamini-Hochberg FDR applied per method; consensus requires agreement of three methods; permutation null gives an empirical FDR."),
    ]

    # ---- markdown ----
    md = []
    md.append(f"# {cfg['project_name']} — run report\n")
    md.append(f"_Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}_\n")
    md.append("## 1. Dataset & QC\n")
    if qc:
        md.append(f"- Genomes passing QC: **{num(qc.get('n_passed'))} / {num(qc.get('n_input'))}**")
        md.append(f"- Per niche: {qc.get('passed_per_niche')}")
        md.append(f"- HQ subset per niche: {qc.get('hq_per_niche')}")
        md.append(f"- QC thresholds: {qc.get('qc_thresholds')}\n")
    md.append("## 2. Confounders and how each is controlled\n")
    for name, how in confounders:
        md.append(f"- **{name}** — {how}")
    md.append("\n## 3. Phylogenetic scaffold coverage\n")
    if cov:
        md.append(f"- Coverage by primary niche: "
                  f"{ {k: v['fraction_placed'] for k, v in cov.get('by_primary_niche', {}).items()} }")
        md.append(f"- Unplaced-species policy: {cov.get('unplaced_policy')}\n")
    md.append("## 4. Consensus signatures (per layer × contrast)\n")
    if not sig_df.empty:
        md.append(sig_df.to_markdown(index=False))
    md.append("\n## 5. Permutation empirical FDR\n")
    if not perm_df.empty:
        md.append(perm_df.to_markdown(index=False))
    md.append("\n## 6. Scoary2 pan-GWAS (4th method)\n")
    md.append("Per-contrast Scoary2 hit counts and the methods entering each "
              "consensus are in the signatures table (columns `scoary_hits`, "
              "`methods`). Clade-stratified (per-genus) ortholog-family results:\n")
    md.append(clade_sum.to_markdown(index=False) if not clade_sum.empty else "(clade run not present)")
    md.append("\n## 7. Functional enrichment of niche signatures\n")
    md.append("Categories over-represented among the consensus niche-signature "
              "features, significant by BOTH over-representation (hypergeometric, "
              "tested-feature background) and preranked GSEA:\n")
    md.append(enrich_top.to_markdown(index=False) if not enrich_top.empty else "(enrichment not present)")
    md.append("\n## 8. Figures\n")
    for f in figs:
        rel = os.path.relpath(f, os.path.dirname(args.out_prefix))
        md.append(f"### {os.path.basename(f)}\n\n![{os.path.basename(f)}]({rel})\n")
    md_text = "\n".join(md)
    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True)
    open(f"{args.out_prefix}.md", "w").write(md_text)

    # ---- minimal HTML ----
    def tbl(df):
        return df.to_html(index=False, border=0) if not df.empty else "<p>(none)</p>"
    rows = "".join(
        f"<h3>{os.path.basename(f)}</h3>"
        f"<img src='{os.path.relpath(f, os.path.dirname(args.out_prefix))}' "
        f"style='max-width:100%;border:1px solid #ddd'>"
        for f in figs)
    conf = "".join(f"<li><b>{n}</b> — {h}</li>" for n, h in confounders)
    html = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>{cfg['project_name']} report</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;max-width:1000px;margin:2rem auto;
line-height:1.45;color:#222}}table{{border-collapse:collapse;margin:1rem 0}}
th,td{{padding:4px 10px;border-bottom:1px solid #eee;text-align:right}}
th:first-child,td:first-child{{text-align:left}}code{{background:#f4f4f4;padding:1px 4px}}</style>
</head><body>
<h1>{cfg['project_name']} — run report</h1>
<p><i>Generated {provenance_stamp(cfg)['utc']} · seed {cfg['seed']}</i></p>
<h2>1. Dataset &amp; QC</h2>
<p>Genomes passing QC: <b>{num(qc.get('n_passed'))}</b> of {num(qc.get('n_input'))}<br>
Per niche: {qc.get('passed_per_niche')}<br>HQ per niche: {qc.get('hq_per_niche')}</p>
<h2>2. Confounders &amp; controls</h2><ul>{conf}</ul>
<h2>3. Scaffold coverage</h2>
<p>{ {k: v['fraction_placed'] for k, v in cov.get('by_primary_niche', {}).items()} if cov else 'NA' }</p>
<h2>4. Consensus signatures</h2>{tbl(sig_df)}
<h2>5. Permutation empirical FDR</h2>{tbl(perm_df)}
<h2>6. Scoary2 pan-GWAS (4th method)</h2>
<p>Per-contrast Scoary2 hits and the methods in each consensus are in the
signatures table (<code>scoary_hits</code>, <code>methods</code>). Clade-stratified
(per-genus) ortholog-family results:</p>{tbl(clade_sum)}
<h2>7. Functional enrichment of niche signatures</h2>
<p>Categories over-represented among the consensus signatures, significant by both
over-representation (hypergeometric, tested-feature background) and preranked
GSEA:</p>{tbl(enrich_top)}
<h2>8. Figures</h2>{rows}
</body></html>"""
    open(f"{args.out_prefix}.html", "w").write(html)
    print(f"Report written: {args.out_prefix}.html / .md")


if __name__ == "__main__":
    main()
