#!/usr/bin/env python3
"""Derive a panel-only view from the genome-wide AnnotSV TSV.

AnnotSV is run once genome-wide.  This script creates a post-hoc panel view by
matching gene symbols case-insensitively, including rows that contain multiple
genes separated by comma/semicolon/slash/pipe delimiters.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)

GENE_COL_CANDIDATES = [
    "Gene_name",
    "Gene_name(s)",
    "GeneName",
    "gene_name",
    "Gene",
    "GENE",
    "Genes",
    "SYMBOL",
]


def load_genes(path: str) -> set[str]:
    genes: set[str] = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene.upper())
    return genes


def find_gene_column(fieldnames: list[str]) -> str | None:
    exact = set(fieldnames)
    for candidate in GENE_COL_CANDIDATES:
        if candidate in exact:
            return candidate
    lower = {name.lower(): name for name in fieldnames}
    for candidate in GENE_COL_CANDIDATES:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Filter genome-wide AnnotSV TSV to the configured candidate-gene panel"
    )
    parser.add_argument("--annotsv", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    panel_genes = load_genes(args.genes)

    with open(
        args.annotsv,
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        gene_col = find_gene_column(fieldnames)
        if gene_col is None:
            print(
                f"[ERROR] no supported gene column found in {args.annotsv}; "
                f"expected one of {GENE_COL_CANDIDATES}; header={fieldnames}",
                file=sys.stderr,
            )
            return 1

        matched: list[dict] = []
        for row in reader:
            value = row.get(gene_col) or ""

            if value in {"", ".", "NA", "N/A", "NONE"}:
                continue

            genes_in_row = {
                gene.strip().upper()
                for gene in re.split(r"[,;/|]", value)
                if gene.strip()
            }

            if genes_in_row & panel_genes:
                matched.append(row)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as fh:
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
        f"[OK] gene_column={gene_col} panel_genes={len(panel_genes)} "
        f"rows_matched={len(matched)} output={output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
