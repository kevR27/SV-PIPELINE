#!/usr/bin/env python3
"""Evaluate custom SV-gene prioritisation against Genomiser using solved cases.

The custom final TSV is organised by panel status. Reconstruct a genome-wide
gene order from its documented final sorting keys, ignoring panel status.
Never treat a Genomiser gene/MOI rank as a separate structural-variant rank.
"""
import argparse
import csv
import json
from pathlib import Path

TOP_K = (1, 5, 10, 20, 50)
TIER_PRIORITY = {"HIGH": 3, "MODERATE": 2, "SUPPORTING": 1, "LIMITED": 0}


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), set(reader.fieldnames or [])


def get(row, *names):
    for name in names:
        value = str(row.get(name, "") or "").strip()
        if value and value.upper() not in {".", "NA", "N/A", "NAN", "NONE"}:
            return value
    return ""


def rank_number(row, *fields):
    try:
        number = int(float(get(row, *fields)))
        return number if number > 0 else None
    except (ValueError, TypeError):
        return None


def gene(value):
    return value.strip().upper()


def chromosome(value):
    chrom = value.strip().upper()
    if chrom.startswith("CHR"):
        chrom = chrom[3:]
    return "MT" if chrom in {"MT", "M"} else chrom


def sequential_gene_ranks(data, column):
    """Remove multiple mode-of-inheritance rows for the same gene."""
    positions = {}
    for row in data:
        symbol = gene(get(row, column))
        if symbol and symbol not in positions:
            positions[symbol] = len(positions) + 1
    return positions


def custom_global_gene_ranks(data, fields):
    """Reconstruct the final tier/event-based order without panel grouping."""
    required = {"GENE", "FINAL_GENE_RELEVANCE_TIER",
                "BEST_EVENT_FINAL_RANK_GLOBAL", "FINAL_GENE_RELEVANCE_DISPLAY_SCORE"}
    missing = required - fields
    if missing:
        raise ValueError("Custom ranked genes need global comparison fields: "
                         + ", ".join(sorted(missing)))

    def order(row):
        tier = TIER_PRIORITY.get(get(row, "FINAL_GENE_RELEVANCE_TIER").upper(), 0)
        best_event = rank_number(row, "BEST_EVENT_FINAL_RANK_GLOBAL") or 10**12
        try:
            score = float(get(row, "FINAL_GENE_RELEVANCE_DISPLAY_SCORE"))
        except ValueError:
            score = 0.0
        return (-tier, best_event, -score, gene(get(row, "GENE")))

    return sequential_gene_ranks(sorted(data, key=order), "GENE")


def genomiser_global_gene_ranks(data, fields):
    if not {"GENE_SYMBOL", "#RANK"}.issubset(fields):
        raise ValueError("Genomiser genes.tsv must contain GENE_SYMBOL and #RANK")
    return sequential_gene_ranks(
        sorted(data, key=lambda row: rank_number(row, "#RANK") or 10**12),
        "GENE_SYMBOL",
    )


def variant_type(row, method):
    if method == "custom":
        return get(row, "SVTYPE").upper()
    alt = get(row, "ALT")
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].split(":")[0].upper()
    if "[" in alt or "]" in alt:
        return "BND"
    return "SMALL"


def sv_matches(row, truth, method, tolerance):
    """Exploratory symbolic-SV matching; inspect confirmed events manually."""
    symbol = "GENE" if method == "custom" else "GENE_SYMBOL"
    contig = "CHROM" if method == "custom" else "CONTIG"
    if gene(get(row, symbol)) != gene(truth["gene"]):
        return False
    if chromosome(get(row, contig)) != chromosome(truth["chrom"]):
        return False
    requested_type = truth["svtype"].upper()
    if requested_type not in {"DEL", "DUP", "INV", "INS"}:
        return False
    if variant_type(row, method) != requested_type:
        return False
    start = rank_number(row, "START")
    end = rank_number(row, "END")
    if start is None or end is None:
        return False
    try:
        reference_start = int(truth["pos"])
        reference_end = int(truth["end"])
    except ValueError:
        return False
    if abs(start - reference_start) > tolerance:
        return False
    return requested_type == "INS" or abs(end - reference_end) <= tolerance


def output_table(path, data):
    columns = [
        "sample_id", "gene", "variant_class", "truth_status",
        "custom_gene_rank_global_reconstructed", "genomiser_gene_rank",
        "custom_sv_reported", "genomiser_sv_reported",
        "custom_sv_event_rank", "sv_matching_status"
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data)


def fraction(data, predicate):
    return sum(predicate(row) for row in data) / len(data) if data else None


def top_k(data, rank_column, require_sv=None):
    return {
        str(k): fraction(data, lambda row: (
            isinstance(row[rank_column], int)
            and row[rank_column] <= k
            and (require_sv is None or row[require_sv] == "YES")
        ))
        for k in TOP_K
    }


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

    truth, _ = read_tsv(args.truth)
    truth = [
        row for row in truth if get(row, "sample_id") == args.sample
        and get(row, "truth_status").upper() == "CONFIRMED"
    ]
    if not truth:
        raise ValueError("No independently confirmed truth entries for " + args.sample)
    if any(not get(row, "gene") for row in truth):
        raise ValueError("Truth entries must include the confirmed gene")

    custom_genes, custom_fields = read_tsv(args.custom_genes)
    custom_svs, custom_sv_fields = read_tsv(args.custom_sv)
    genomiser_genes, genomiser_fields = read_tsv(args.genomiser_genes)
    genomiser_variants, genomiser_variant_fields = read_tsv(args.genomiser_variants)
    if not {"GENE", "CHROM", "START", "END", "SVTYPE"}.issubset(custom_sv_fields):
        raise ValueError("Custom SV results need GENE, CHROM, START, END and SVTYPE")
    if not {"GENE_SYMBOL", "CONTIG", "START", "END", "ALT"}.issubset(genomiser_variant_fields):
        raise ValueError("Genomiser variants.tsv is missing required SV comparison columns")

    custom_ranks = custom_global_gene_ranks(custom_genes, custom_fields)
    genomiser_ranks = genomiser_global_gene_ranks(genomiser_genes, genomiser_fields)
    details = []

    for confirmed in truth:
        symbol = gene(get(confirmed, "gene"))
        variant_class = get(confirmed, "variant_class").upper() or "GENE_ONLY"
        custom_present = "NOT_EVALUATED"
        genomiser_present = "NOT_EVALUATED"
        event_rank = None
        status = "GENE_LEVEL_ONLY"

        if variant_class == "SV":
            required = ("chrom", "pos", "end", "svtype")
            if not all(get(confirmed, field) for field in required):
                status = "SV_TRUTH_COORDINATES_MISSING"
            elif get(confirmed, "svtype").upper() not in {"DEL", "DUP", "INV", "INS"}:
                status = "SV_TYPE_UNSUPPORTED"
            else:
                custom_hits = [
                    row for row in custom_svs
                    if sv_matches(row, confirmed, "custom", args.sv_tolerance)
                ]
                genomiser_hits = [
                    row for row in genomiser_variants
                    if sv_matches(row, confirmed, "genomiser", args.sv_tolerance)
                ]
                custom_present = "YES" if custom_hits else "NO"
                genomiser_present = "YES" if genomiser_hits else "NO"
                event_positions = [
                    rank_number(row, "FINAL_EVENT_RANK_GLOBAL", "EVENT_RANK_GLOBAL")
                    for row in custom_hits
                ]
                event_rank = min((n for n in event_positions if n is not None),
                                 default=None)
                status = "SYMBOLIC_SV_COORDINATE_MATCH"
        elif variant_class != "GENE_ONLY":
            status = "CUSTOM_PIPELINE_VARIANT_CLASS_NOT_COMPARABLE"

        details.append({
            "sample_id": args.sample,
            "gene": symbol,
            "variant_class": variant_class,
            "truth_status": "CONFIRMED",
            "custom_gene_rank_global_reconstructed": custom_ranks.get(symbol) or "NOT_RANKED",
            "genomiser_gene_rank": genomiser_ranks.get(symbol) or "NOT_RANKED",
            "custom_sv_reported": custom_present,
            "genomiser_sv_reported": genomiser_present,
            "custom_sv_event_rank": event_rank or "NOT_RANKED",
            "sv_matching_status": status
        })

    output_table(args.detail_output, details)

    # Gene-level recovery is meaningful for all independently confirmed
    # gene diagnoses, but a custom SV-only ranking has a different candidate
    # universe from Genomiser. Analyse solved SVs separately as primary cohort.
    comparable_svs = [
        row for row in details
        if row["sv_matching_status"] == "SYMBOLIC_SV_COORDINATE_MATCH"
    ]
    results = {
        "sample_id": args.sample,
        "ground_truth_records": len(details),
        "gene_rank_note": (
            "Custom rank reconstructed without PANEL_STATUS using tier, "
            "best event global rank, phenotype score and gene symbol"
        ),
        "all_confirmed_gene_top_k": {
            "denominator": len(details),
            "custom": top_k(details, "custom_gene_rank_global_reconstructed"),
            "genomiser": top_k(details, "genomiser_gene_rank")
        },
        "solved_sv_top_k": {
            "denominator": len(comparable_svs),
            "custom": top_k(comparable_svs, "custom_gene_rank_global_reconstructed"),
            "genomiser": top_k(comparable_svs, "genomiser_gene_rank")
        },
        "solved_sv_allele_reported_within_top_k_gene": {
            "denominator": len(comparable_svs),
            "custom": top_k(comparable_svs, "custom_gene_rank_global_reconstructed",
                            "custom_sv_reported"),
            "genomiser": top_k(comparable_svs, "genomiser_gene_rank",
                               "genomiser_sv_reported")
        },
        "sv_tolerance_bp": args.sv_tolerance,
        "limitations": (
            "Genomiser ranks gene/MOI models; its #RANK must not be used as an "
            "independent allele rank. Structural variants require manual match "
            "review (especially complex BND and variable breakpoint precision). "
            "Absent Genomiser variants can be filtered rather than uncalled. "
            "Custom ranking prioritises SV-affected genes, not all SNV-only genes."
        )
    }
    Path(args.summary_output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_output).write_text(
        json.dumps(results, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
