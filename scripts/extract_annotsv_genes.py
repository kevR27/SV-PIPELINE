#!/usr/bin/env python3
"""Extract unique gene symbols from an AnnotSV TSV.

AnnotSV column names can vary slightly by version/output mode, so the script
selects the first supported gene column present rather than silently returning
an empty list when `Gene_name` is not available.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

MISSING = {"", ".", "NA", "N/A", "None"}
GENE_COLUMNS = ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotsv", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    try:
        csv.field_size_limit(sys.maxsize)
    except OverflowError:
        csv.field_size_limit(2**31 - 1)

    genes = set()
    with open(args.annotsv, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        gene_col = next((c for c in GENE_COLUMNS if c in fields), None)
        if gene_col is None:
            raise ValueError(
                "No supported AnnotSV gene column found. Available columns: "
                + ", ".join(fields)
            )

        for row in reader:
            value = str(row.get(gene_col, "")).strip()
            if value in MISSING:
                continue
            for gene in re.split(r"[;,|/]", value):
                gene = gene.strip()
                if gene not in MISSING:
                    genes.add(gene.upper())

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as fh:
        for gene in sorted(genes):
            fh.write(gene + "\n")

    print(f"[OK] genes={len(genes)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
