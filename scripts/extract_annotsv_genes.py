#!/usr/bin/env python3
"""Extract unique gene symbols from the Gene_name column of an AnnotSV TSV."""

import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("--annotsv", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

genes = set()

with open(args.annotsv, newline="", encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    if "Gene_name" not in (reader.fieldnames or []):
        raise SystemExit("ERROR: AnnotSV output has no Gene_name column.")

    for row in reader:
        value = row.get("Gene_name", "").strip()
        if not value or value == ".":
            continue

        # AnnotSV may represent multiple genes in one field.
        for gene in re.split(r"[;,|]", value):
            gene = gene.strip()
            if gene and gene not in {".", "NA", "N/A"}:
                genes.add(gene)

with open(args.output, "w", encoding="utf-8") as out:
    for gene in sorted(genes):
        out.write(gene + "\n")
