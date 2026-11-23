#!/usr/bin/env python3
"""Independent, per-sample benchmark of custom SV ranking vs Genomiser.

A truth set must be independently validated and prepared before inspecting
either ranking. This program does not optimize ranking weights.
"""
import argparse
import csv
import json
from pathlib import Path

CUTS = (1, 5, 10, 20, 50)


def rows(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def value(row, *names):
    for name in names:
        result = row.get(name)
        if result is not None and str(result).strip() not in {"", ".", "NA", "nan"}:
            return str(result).strip()
    return ""


def gene(value_):
    return value_.strip().upper()


def chrom(value_):
    return value_.strip().upper().removeprefix("CHR").replace("M", "MT")


def positive_int(row, *fields):
    item = value(row, *fields)
    try:
        return int(float(item))
    except (ValueError, TypeError):
        return None


def gene_ranks(data, gene_col):
    ranks = {}
    # Use the file's ranked order. Multiple inheritance modes for one gene
    # count as one gene in a fair gene-level top-K comparison.
    for row in data:
        name = gene(value(row, gene_col))
        if name and name not in ranks:
            ranks[name] = len(ranks) + 1
    return ranks


def type_of(row, method):
    if method == "custom":
        return value(row, "SVTYPE").upper()
    alt = value(row, "ALT")
    if alt.startswith("<") and alt.endswith(">"):
        return alt.strip("<>").split(":")[0].upper()
    if "[" in alt or "]" in alt:
        return "BND"
    return "SMALL"


def sv_match(row, truth, method, tolerance):
    if gene(value(row, "GENE" if method == "custom" else "GENE_SYMBOL")) != gene(truth["gene"]):
        return False
    start = positive_int(row, "START")
    end = positive_int(row, "END")
    if start is None or end is None:
        return False
    if chrom(value(row, "CHROM" if method == "custom" else "CONTIG")) != chrom(truth["chrom"]):
        return False
    svtype = truth["svtype"].upper()
    observed = type_of(row, method)
    # A BND/complex allele is not benchmarked with imprecise interval matching.
    if svtype not in {"DEL", "DUP", "INV", "INS"} or observed != svtype:
        return False
    p = int(truth["pos"])
    q = int(truth["end"])
    if abs(start - p) > tolerance:
        return False
    return svtype == "INS" or abs(end - q) <= tolerance


def write_tsv(path, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "sample_id", "gene", "truth_status", "variant_class",
            "custom_gene_rank", "genomiser_gene_rank",
            "custom_sv_event_rank", "genomiser_sv_parent_gene_rank",
            "variant_comparison_status"
        ], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--custom-genes", required=True)
    parser.add_argument("--custom-sv", required=True)
    parser.add_argument("--genomiser-genes", required=True)
    parser.add_argument("--genomiser-variants", required=True)
    parser.add_argument("--detail-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--sv-tolerance", type=int, default=500)
    args = parser.parse_args()

    truth = [
        row for row in rows(args.truth)
        if value(row, "sample_id") == args.sample
        and value(row, "truth_status").upper() == "CONFIRMED"
    ]
    if not truth:
        raise ValueError(f"No CONFIRMED truth records for sample {args.sample}")
    for row in truth:
        if not value(row, "gene"):
            raise ValueError("Every truth record needs a gene")
    custom_genes = rows(args.custom_genes)
    custom_sv = rows(args.custom_sv)
    genomiser_genes = rows(args.genomiser_genes)
    genomiser_variants = rows(args.genomiser_variants)
    if custom_genes and "GENE" not in custom_genes[0]:
        raise ValueError("Custom gene ranking must have GENE")
    if genomiser_genes and "GENE_SYMBOL" not in genomiser_genes[0]:
        raise ValueError("Genomiser gene ranking must have GENE_SYMBOL")
    custom_rank = gene_ranks(custom_genes, "GENE")
    genomiser_rank = gene_ranks(genomiser_genes, "GENE_SYMBOL")
    output = []
    for t in truth:
        target_gene = gene(t["gene"])
        cr = custom_rank.get(target_gene)
        gr = genomiser_rank.get(target_gene)
        category = value(t, "variant_class").upper() or "GENE_ONLY"
        event_rank = None
        genomiser_event = None
        match_status = "GENE_LEVEL_ONLY"
        if category == "SV":
            required = [value(t, x) for x in ("chrom", "pos", "end", "svtype")]
            if not all(required):
                match_status = "SV_COORDINATES_MISSING"
            elif value(t, "svtype").upper() not in {"DEL", "DUP", "INV", "INS"}:
                match_status = "SV_TYPE_NOT_SUPPORTED_BY_BENCHMARK"
            else:
                candidates = [
                    positive_int(r, "FINAL_EVENT_RANK_GLOBAL", "EVENT_RANK_GLOBAL")
                    for r in custom_sv if sv_match(r, t, "custom", args.sv_tolerance)
                ]
                event_rank = min((rank for rank in candidates if rank is not None), default=None)
                candidates = [
                    positive_int(r, "#RANK", "RANK")
                    for r in genomiser_variants if sv_match(r, t, "genomiser", args.sv_tolerance)
                ]
                genomiser_event = min((rank for rank in candidates if rank is not None), default=None)
                match_status = "SV_COORDINATE_MATCH_EVALUATED"
        elif category != "GENE_ONLY":
            match_status = "VARIANT_CLASS_NOT_COMPARABLE_TO_CUSTOM_SV"
        output.append({
            "sample_id": args.sample,
            "gene": target_gene,
            "truth_status": "CONFIRMED",
            "variant_class": category,
            "custom_gene_rank": cr if cr is not None else "NOT_RANKED",
            "genomiser_gene_rank": gr if gr is not None else "NOT_RANKED",
            "custom_sv_event_rank": event_rank if event_rank is not None else "NOT_RANKED",
            "genomiser_sv_parent_gene_rank": genomiser_event if genomiser_event is not None else "NOT_RANKED",
            "variant_comparison_status": match_status
        })

    write_tsv(args.detail_output, output)
    gene_summary = {}
    for method, column in (("custom", "custom_gene_rank"), ("genomiser", "genomiser_gene_rank")):
        gene_summary[method] = {
            "denominator": len(output),
            "top_k_recall": {
                str(k): sum(isinstance(row[column], int) and row[column] <= k for row in output) / len(output)
                for k in CUTS
            }
        }
    sv_truth = [row for row in output if row["variant_comparison_status"] == "SV_COORDINATE_MATCH_EVALUATED"]
    sv_summary = {}
    for method, column in (("custom", "custom_sv_event_rank"),
                           ("genomiser", "genomiser_sv_parent_gene_rank")):
        sv_summary[method] = {
            "denominator": len(sv_truth),
            "top_k_recall": {
                str(k): (sum(isinstance(row[column], int) and row[column] <= k for row in sv_truth)
                         / len(sv_truth) if sv_truth else None)
                for k in CUTS
            }
        }
    summary = {
        "sample_id": args.sample,
        "cohort_note": "Aggregate across independently confirmed samples before reporting sensitivity.",
        "truth_records": len(output),
        "gene_level": gene_summary,
        "structural_variant_level": sv_summary,
        "sv_matching_tolerance_bp": args.sv_tolerance,
        "genomiser_variant_rank_interpretation": "#RANK is the parent gene/MOI score rank, not an independent allele score rank",
        "caution": "Methods have differing candidate universes; report callability and missing variants separately."
    }
    Path(args.summary_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_output).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
