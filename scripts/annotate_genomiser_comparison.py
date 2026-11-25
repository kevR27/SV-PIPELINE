#!/usr/bin/env python3
"""Independent Genomiser gene/model and exact-coordinate symbolic SV comparison."""
from collections import defaultdict
import pandas as pd
from benchmark_genomiser import read_tsv, get, gene, chromosome, rank_number, variant_type, genomiser_global_gene_ranks

def annotate_genomiser(table, genes_path, variants_path):
    genes, fields = read_tsv(genes_path)
    variants, variant_fields = read_tsv(variants_path)
    ranks = genomiser_global_gene_ranks(genes, fields)
    if not {"GENE_SYMBOL", "CONTIG", "START", "END", "ALT"}.issubset(variant_fields):
        raise ValueError("Genomiser variant comparison requires GENE_SYMBOL/CONTIG/START/END/ALT")
    models = {}
    for row in sorted(genes, key=lambda r: rank_number(r, "#RANK") or 10**12):
        symbol = gene(get(row, "GENE_SYMBOL"))
        if symbol and rank_number(row, "#RANK") is not None:
            models.setdefault(symbol, row)
    events = defaultdict(list)
    for row in variants:
        kind = variant_type(row, "genomiser")
        if kind in {"DEL", "DUP", "INV", "INS"}:
            key = (gene(get(row, "GENE_SYMBOL")), chromosome(get(row, "CONTIG")),
                   rank_number(row, "START"), rank_number(row, "END"), kind)
            if key[2] is not None and key[3] is not None:
                events[key].append(row)
    out = table.copy()
    annotations = []
    for row in out.fillna(".").to_dict("records"):
        symbol = gene(get(row, "GENE", "GENES"))
        model = models.get(symbol, {})
        kind = get(row, "SVTYPE").upper()
        key = (symbol, chromosome(get(row, "CHROM")), rank_number(row, "START"),
               rank_number(row, "END"), kind)
        matches = events.get(key, [])
        status = ("NOT_EVALUATED_SV_TYPE" if kind not in {"DEL", "DUP", "INV", "INS"}
                  else "EXACT_COORDINATE_CANDIDATE" if len(matches) == 1
                  else "AMBIGUOUS_COORDINATE_MATCH" if matches else "NOT_REPORTED_AT_EXACT_COORDINATES")
        annotations.append({
            "GENOMISER_GENE_RANK": ranks.get(symbol, "."),
            "GENOMISER_MODEL_RANK": get(model, "#RANK") or ".",
            "GENOMISER_GENE_SCORE": get(model, "GENE_SCORE", "EXOMISER_GENE_SCORE") or ".",
            "GENOMISER_GENE_PHENOTYPE_SCORE": get(model, "PHENOTYPE_SCORE", "EXOMISER_GENE_PHENO_SCORE") or ".",
            "GENOMISER_SV_MATCH_STATUS": status,
            "GENOMISER_SV_VARIANT_SCORE": (get(matches[0], "VARIANT_SCORE", "EXOMISER_VARIANT_SCORE") or ".") if len(matches) == 1 else ".",
        })
    columns = ["GENOMISER_GENE_RANK", "GENOMISER_MODEL_RANK", "GENOMISER_GENE_SCORE",
               "GENOMISER_GENE_PHENOTYPE_SCORE", "GENOMISER_SV_MATCH_STATUS", "GENOMISER_SV_VARIANT_SCORE"]
    for col in columns:
        out[col] = [row[col] for row in annotations]
    return out
