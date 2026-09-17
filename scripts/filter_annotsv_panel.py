#!/usr/bin/env python3
"""Derive a panel-only view from the genome-wide AnnotSV TSV.

AnnotSV is executed only once genome-wide. This script retains every annotation
row whose gene field intersects the candidate-gene list.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

GENE_COL_CANDIDATES = [
    "Gene_name", "Gene_name(s)", "GeneName", "gene_name",
    "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene",
]


def load_genes(path: str) -> set[str]:
    genes = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene.upper())
    return genes


def find_gene_column(fieldnames: list[str]) -> str | None:
    for candidate in GENE_COL_CANDIDATES:
        if candidate in fieldnames:
            return candidate
    lowered = {x.lower(): x for x in fieldnames}
    for candidate in GENE_COL_CANDIDATES:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def split_genes(value: str) -> set[str]:
    return {
        gene.strip().upper()
        for gene in re.split(r"[;,|/]", value or "")
        if gene.strip() and gene.strip() not in {".", "NA", "N/A"}
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="Filter genome-wide AnnotSV TSV to a candidate gene panel")
    ap.add_argument("--annotsv", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    panel = load_genes(args.genes)
    with open(args.annotsv, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        gene_col = find_gene_column(fieldnames)
        if gene_col is None:
            print(
                f"[ERROR] no supported gene-symbol column found in {args.annotsv}; "
                f"header={fieldnames}",
                file=sys.stderr,
            )
            return 1

        matched = [row for row in reader if split_genes(row.get(gene_col, "")) & panel]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(matched)

    print(
        f"[OK] gene_column={gene_col} panel_genes={len(panel)} "
        f"rows_matched={len(matched)} output={out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
