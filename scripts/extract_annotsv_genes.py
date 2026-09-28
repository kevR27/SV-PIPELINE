#!/usr/bin/env python3
"""Extract unique gene symbols from a genome-wide AnnotSV TSV."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from sv_evidence_common import gene_symbols, invalid_gene_labels

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)

GENE_COLUMNS = [
    "Gene_name",
    "Gene_name(s)",
    "GeneName",
    "gene_name",
    "Gene",
    "GENE",
    "Genes",
    "SYMBOL",
]


def find_gene_column(fieldnames: list[str]) -> str | None:
    exact = set(fieldnames)
    for name in GENE_COLUMNS:
        if name in exact:
            return name
    lower = {name.lower(): name for name in fieldnames}
    for candidate in GENE_COLUMNS:
        if candidate.lower() in lower:
            return lower[candidate.lower()]
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotsv", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genes: set[str] = set()
    invalid_rows = 0
    with open(args.annotsv, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        gene_col = find_gene_column(fieldnames)
        if gene_col is None:
            raise ValueError(
                f"No supported AnnotSV gene column found in {args.annotsv}. "
                f"Expected one of {GENE_COLUMNS}; header={fieldnames}"
            )

        for row in reader:
            value = row.get(gene_col, "")
            if not value or value == ".":
                continue
            invalid_rows += bool(invalid_gene_labels(value))
            genes.update(gene_symbols(value))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        for gene in sorted(genes):
            out.write(gene + "\n")

    print(f"[OK] gene_column={gene_col} genes={len(genes)} output={output}")
    if invalid_rows:
        print(f"[WARN] {invalid_rows} AnnotSV rows contain cmpl/incmpl as gene names. "
              "These labels were excluded from the gene list. Check the AnnotSV gene resource; "
              "the corresponding gene names are unresolved.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
