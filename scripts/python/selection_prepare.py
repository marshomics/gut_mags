#!/usr/bin/env python3
"""
selection_prepare.py
--------------------
Prepare HyPhy inputs for the human-associated ortholog families of one clade:
codon-aware alignments and a gene tree with the human-associated (foreground)
branches labelled. This is the substrate for testing sequence-level selection,
the third mode of adaptation (the gene is shared but evolves differently in the
human lineage).

For each top family from the clade's Scoary2 result:
  * gather member CDS from the Prokka nucleotide files (gene ids in the Panaroo
    gene_presence_absence.csv);
  * align the translated proteins (MAFFT) and back-translate to a codon
    alignment (preserves reading frame for dN/dS);
  * build a FastTree gene tree and label foreground tips ({Foreground}) =
    genomes of the foreground niche (human by default).

Outputs per family under families/<family>/: codon.fasta, tree.nwk
plus families.tsv (manifest).
"""
import argparse
import os
import re
import subprocess
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from hgn_utils import load_config, get_logger, load_annotation_paths

log = get_logger("sel-prep")


def read_ffn(path):
    if not os.path.exists(path):
        return {}
    return {r.id: str(r.seq) for r in SeqIO.parse(path, "fasta")}


def mafft_protein(prots):
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
        for name, seq in prots.items():
            fh.write(f">{name}\n{seq}\n")
        infile = fh.name
    try:
        out = subprocess.run(["mafft", "--auto", "--quiet", infile],
                             capture_output=True, text=True, check=True).stdout
    except Exception as e:
        log.warning("mafft failed: %s", e); os.unlink(infile); return {}
    os.unlink(infile)
    aln, name, buf = {}, None, []
    for line in out.splitlines():
        if line.startswith(">"):
            if name:
                aln[name] = "".join(buf)
            name = line[1:].strip(); buf = []
        else:
            buf.append(line.strip())
    if name:
        aln[name] = "".join(buf)
    return aln


def backtranslate(prot_aln, nt):
    """Thread codons from nt onto the gapped protein alignment."""
    out = {}
    for name, pa in prot_aln.items():
        seq = nt.get(name, "")
        codons, i = [], 0
        for aa in pa:
            if aa == "-":
                codons.append("---")
            else:
                codons.append(seq[i:i + 3] if i + 3 <= len(seq) else "---")
                i += 3
        out[name] = "".join(codons)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--gpa-csv", required=True, help="Panaroo gene_presence_absence.csv")
    ap.add_argument("--scoary", required=True, help="clade Scoary2 parsed table")
    ap.add_argument("--species-genomes", required=True)
    ap.add_argument("--annotation-paths", required=True,
                    help="annotation_paths.tsv written by resolve_annotation_paths.py")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    sel = cfg["synthesis"]["selection"]
    fg_niche = sel["foreground_niche"]
    apaths = load_annotation_paths(args.annotation_paths)
    os.makedirs(f"{args.out_dir}/families", exist_ok=True)

    sg = pd.read_csv(args.species_genomes, sep="\t")
    sg["genome"] = sg["genome"].astype(str)
    fg_genomes = set(sg.loc[sg["niche"] == fg_niche, "genome"])

    sc = pd.read_csv(args.scoary, sep="\t") if os.path.getsize(args.scoary) else pd.DataFrame()
    top = []
    if len(sc) and "scoary_sig" in sc.columns:
        sc2 = sc[sc["scoary_sig"].astype(str).isin(["True", "TRUE", "1"])]
        sc2 = sc2.sort_values("scoary_fisher_q") if "scoary_fisher_q" in sc2.columns else sc2
        top = sc2["feature"].astype(str).head(sel["top_families_per_clade"]).tolist()
    if not top:
        pd.DataFrame(columns=["family", "n_seqs", "n_foreground"]).to_csv(
            f"{args.out_dir}/families.tsv", sep="\t", index=False)
        log.warning("No significant families for selection in this clade."); return

    gpa = pd.read_csv(args.gpa_csv, sep="\t" if "\t" in open(args.gpa_csv).readline() else ",",
                      dtype=str, keep_default_na=False)
    gene_col = gpa.columns[0]
    genome_cols = [c for c in gpa.columns if c in set(sg["genome"])]
    ffn_cache = {}
    manifest = []
    for fam in top:
        row = gpa[gpa[gene_col] == fam]
        if row.empty:
            continue
        row = row.iloc[0]
        nt = {}
        for g in genome_cols:
            gid = str(row[g]).split(";")[0].strip()
            if not gid or gid in ("", "nan"):
                continue
            if g not in ffn_cache:
                fp = apaths.get(g, {}).get("prokka_ffn")
                if not fp:
                    log.warning("no prokka_ffn path for genome %s; skipping it", g)
                    ffn_cache[g] = {}
                else:
                    ffn_cache[g] = read_ffn(fp)
            seq = ffn_cache[g].get(gid)
            if seq and len(seq) >= 60 and len(seq) % 3 == 0:
                nt[g] = seq
        if len(nt) < 6:
            continue
        prots = {g: str(Seq(s).translate(table=11)).rstrip("*") for g, s in nt.items()}
        prot_aln = mafft_protein(prots)
        if not prot_aln:
            continue
        codon = backtranslate(prot_aln, nt)
        famdir = f"{args.out_dir}/families/{re.sub(r'[^A-Za-z0-9]+','_',fam)}"
        os.makedirs(famdir, exist_ok=True)
        with open(f"{famdir}/codon.fasta", "w") as fh:
            for g, s in codon.items():
                fh.write(f">{g}\n{s}\n")
        # FastTree (nucleotide) then label foreground tips
        try:
            tre = subprocess.run(["fasttree", "-nt", "-gtr", "-quiet", f"{famdir}/codon.fasta"],
                                 capture_output=True, text=True, check=True).stdout
        except Exception as e:
            log.warning("fasttree failed for %s: %s", fam, e); continue
        for g in nt:
            if g in fg_genomes:
                tre = re.sub(rf"(\b{re.escape(g)}\b)(?![^,()]*\{{)", r"\1{Foreground}", tre, count=1)
        open(f"{famdir}/tree.nwk", "w").write(tre)
        manifest.append({"family": fam, "n_seqs": len(nt),
                         "n_foreground": sum(g in fg_genomes for g in nt)})
    pd.DataFrame(manifest).to_csv(f"{args.out_dir}/families.tsv", sep="\t", index=False)
    log.info("Prepared %d families for HyPhy (foreground=%s)", len(manifest), fg_niche)


if __name__ == "__main__":
    main()
