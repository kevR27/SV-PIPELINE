#!/usr/bin/env python3
"""Descriptive rank agreement, including unsolved samples; not diagnostic accuracy."""
import argparse
import csv
import json
from pathlib import Path
from benchmark_genomiser import (read_tsv, custom_global_gene_ranks,
                                 genomiser_global_gene_ranks, TOP_K, gene, get)


def compare(custom, genomiser, ks=TOP_K):
    shared = set(custom) & set(genomiser)
    metrics = {}
    for k in ks:
        left = {g for g, r in custom.items() if r <= k}
        right = {g for g, r in genomiser.items() if r <= k}
        overlap = len(left & right)
        union = len(left | right)
        metrics[str(k)] = {"custom_n": len(left), "genomiser_n": len(right),
                           "shared_n": overlap, "jaccard": overlap / union if union else None}
    return {"custom_genes": len(custom), "genomiser_genes": len(genomiser),
            "shared_genes": len(shared), "top_k": metrics,
            "interpretation": "Descriptive agreement only; no truth-based sensitivity, specificity or false-positive rate."}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--custom-genes", required=True)
    p.add_argument("--genomiser-genes", required=True)
    p.add_argument("--details", required=True)
    p.add_argument("--summary", required=True)
    a = p.parse_args()
    custom_rows, custom_fields = read_tsv(a.custom_genes)
    genome_rows, genome_fields = read_tsv(a.genomiser_genes)
    custom = custom_global_gene_ranks(custom_rows, custom_fields)
    genome = genomiser_global_gene_ranks(genome_rows, genome_fields)
    panels = {gene(get(row, "GENE")): get(row, "PANEL_STATUS") for row in custom_rows}
    fields = ["GENE", "CUSTOM_GENE_RANK_GLOBAL", "GENOMISER_GENE_RANK", "CUSTOM_PANEL_STATUS"]
    Path(a.details).parent.mkdir(parents=True, exist_ok=True)
    Path(a.summary).parent.mkdir(parents=True, exist_ok=True)
    with open(a.details, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for symbol in sorted(set(custom) | set(genome), key=lambda g: (custom.get(g, 10**12), genome.get(g, 10**12), g)):
            writer.writerow(dict(zip(fields, [symbol, custom.get(symbol, "."), genome.get(symbol, "."), panels.get(symbol, ".")])))
    result = compare(custom, genome)
    result["sample_id"] = a.sample
    result["candidate_universe_note"] = "Custom ranking uses SV-associated genes; Genomiser includes supported coding and noncoding sequence variants as well as SVs."
    Path(a.summary).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
