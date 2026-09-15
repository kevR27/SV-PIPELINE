#!/usr/bin/env python3
"""Remove optic-neuropathy panel genes from a genome-wide SV gene list."""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--genes", required=True)
parser.add_argument("--panel", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

def read_genes(path):
    genes = set()
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            gene = line.strip().split()[0] if line.strip() else ""
            if gene and not gene.startswith("#"):
                genes.add(gene.upper())
    return genes

all_genes = read_genes(args.genes)
panel_genes = read_genes(args.panel)

with open(args.output, "w", encoding="utf-8") as out:
    for gene in sorted(all_genes - panel_genes):
        out.write(gene + "\n")
