#!/usr/bin/env python3
"""
filter_annotsv_panel.py

Derives the panel-only AnnotSV view by filtering the already-computed
genome-wide AnnotSV TSV on the candidate gene list, instead of re-running
AnnotSV a second time with -candidateGenesFiltering 1.

AnnotSV's gene symbol column name has varied slightly across versions
(commonly "Gene_name"); this script checks a small list of known aliases
and fails loudly if none is found, rather than silently emitting an
(almost) empty file.
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

GENE_COL_CANDIDATES = ["Gene_name", "Gene_name(s)", "GeneName", "gene_name", "SYMBOL"]


def load_genes(path):
    genes = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith("#"):
                genes.add(g)
    return genes


def find_gene_column(fieldnames):
    for cand in GENE_COL_CANDIDATES:
        if cand in fieldnames:
            return cand
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Filter genome-wide AnnotSV TSV to a candidate gene panel.")
    ap.add_argument("--annotsv", required=True, help="Genome-wide AnnotSV output TSV")
    ap.add_argument("--genes", required=True, help="Candidate gene list, one symbol per line")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    panel_genes = load_genes(a.genes)

    with open(a.annotsv, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        gene_col = find_gene_column(fieldnames)
        if gene_col is None:
            print(f"[ERROR] none of the expected gene-symbol columns "
                  f"{GENE_COL_CANDIDATES} found in {a.annotsv}; header was: "
                  f"{fieldnames}", file=sys.stderr)
            return 1

        matched = []
        for row in reader:
            genes_in_row = {g.strip() for g in row.get(gene_col, "").replace(";", "/").split("/") if g.strip()}
            if genes_in_row & panel_genes:
                matched.append(row)

    out = Path(a.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(matched)

    print(f"[OK] gene_column={gene_col} panel_genes={len(panel_genes)} "
          f"rows_matched={len(matched)} output={out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
