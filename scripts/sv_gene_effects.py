#!/usr/bin/env python3
"""Human-readable interpretation of how an SV relates to one gene.

The labels describe genomic geometry. They are not pathogenicity classes.
For inversions and breakends, breakpoint effects are kept separate from genes
that are only located inside a rearranged interval.
"""

from __future__ import annotations

import json
from typing import Any

MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def chromosome(value: Any) -> str | None:
    text = str(value or "").strip()
    if text.upper() in MISSING:
        return None
    return text if text.startswith("chr") else f"chr{text}"


def annotation_rows(row: dict[str, Any], gene: str) -> list[dict[str, Any]]:
    raw = row.get("ANNOTSV_GENE_ROWS_JSON", "[]")
    try:
        records = json.loads(raw) if str(raw).strip().upper() not in MISSING else []
    except (TypeError, json.JSONDecodeError):
        records = []

    gene_upper = gene.upper()
    return [
        record
        for record in records
        if str(record.get("Gene_name", "")).strip().upper() == gene_upper
    ]


def transcript_bounds(record: dict[str, Any], fallback_chrom: Any):
    start = number(record.get("Tx_start"))
    end = number(record.get("Tx_end"))
    chrom = chromosome(record.get("SV_chrom", fallback_chrom))
    if start is None or end is None or chrom is None:
        return None
    if end < start:
        start, end = end, start
    return chrom, int(start), int(end)


def sv_breakpoints(row: dict[str, Any]) -> list[tuple[str, int]]:
    svtype = str(row.get("SVTYPE", "")).upper()
    chrom1 = chromosome(row.get("CHROM"))
    chrom2 = chromosome(row.get("CHR2"))
    start = number(row.get("START"))
    end = number(row.get("END"))
    pos2 = number(row.get("POS2"))

    points: list[tuple[str, int]] = []
    if chrom1 and start is not None:
        points.append((chrom1, int(start)))

    if svtype in {"INV", "BND", "TRA"}:
        second = pos2 if pos2 is not None else end
        if second is not None:
            points.append((chrom2 or chrom1, int(second)))
    elif svtype == "INS":
        pass
    elif chrom1 and end is not None:
        points.append((chrom1, int(end)))

    return [point for point in points if point[0] is not None]


def distance_to_transcript(
    breakpoint: tuple[str, int],
    bounds: tuple[str, int, int],
) -> int | None:
    bp_chrom, bp = breakpoint
    tx_chrom, tx_start, tx_end = bounds
    if bp_chrom != tx_chrom:
        return None
    if tx_start <= bp <= tx_end:
        return 0
    return min(abs(bp - tx_start), abs(bp - tx_end))


def breakpoint_location(records: list[dict[str, Any]]) -> str:
    locations = " ".join(str(record.get("Location", "")).lower() for record in records)
    if "exon" in locations:
        return "EXON"
    if "intron" in locations:
        return "INTRON"
    return "TRANSCRIPT"


def get_sv_gene_effect(
    row: dict[str, Any],
    gene: str,
    near_breakpoint_bp: int = 10_000,
) -> tuple[str, int | None]:
    """Return a readable SV-gene relationship and nearest breakpoint distance."""
    svtype = str(row.get("SVTYPE", "")).upper()
    records = annotation_rows(row, gene)
    bounds = [
        item
        for record in records
        if (item := transcript_bounds(record, row.get("CHROM"))) is not None
    ]

    if not bounds:
        return "GENE_EFFECT_UNRESOLVED", None

    breakpoints = sv_breakpoints(row)
    distances = [
        distance
        for bp in breakpoints
        for tx in bounds
        if (distance := distance_to_transcript(bp, tx)) is not None
    ]
    nearest = min(distances) if distances else None

    gene_chrom = bounds[0][0]
    gene_start = min(item[1] for item in bounds)
    gene_end = max(item[2] for item in bounds)
    start = number(row.get("START"))
    end = number(row.get("END"))
    pos2 = number(row.get("POS2"))

    if svtype in {"DEL", "DUP"} and start is not None and end is not None:
        left, right = sorted((start, end))
        whole_gene = left <= gene_start and right >= gene_end
        overlaps_gene = right >= gene_start and left <= gene_end
        word = "DELETION" if svtype == "DEL" else "DUPLICATION"
        if whole_gene:
            return f"WHOLE_GENE_{word}", nearest
        if overlaps_gene:
            return f"PARTIAL_GENE_{word}", nearest
        return f"{word}_NEAR_GENE", nearest

    if svtype == "INS":
        if nearest == 0:
            location = breakpoint_location(records)
            return f"INSERTION_IN_{location}", 0
        if nearest is not None and nearest <= near_breakpoint_bp:
            return "INSERTION_NEAR_GENE", nearest
        return "INSERTION_GENE_EFFECT_UNRESOLVED", nearest

    if svtype in {"INV", "BND", "TRA"}:
        inside_count = sum(distance == 0 for distance in distances)
        prefix = "INVERSION" if svtype == "INV" else "BREAKEND"

        if inside_count >= 2:
            return f"{prefix}_TWO_BREAKPOINTS_IN_GENE", 0
        if inside_count == 1:
            location = breakpoint_location(records)
            return f"{prefix}_BREAKPOINT_IN_{location}", 0
        if nearest is not None and nearest <= near_breakpoint_bp:
            return f"{prefix}_BREAKPOINT_NEAR_GENE", nearest

        if svtype == "INV" and start is not None:
            second = pos2 if pos2 is not None else end
            if second is not None:
                left, right = sorted((start, second))
                if left <= gene_start and right >= gene_end:
                    return "GENE_FULLY_SPANNED_BY_INVERSION", nearest

        return f"{prefix}_GENE_EFFECT_UNRESOLVED", nearest

    return "OTHER_SV_GENE_OVERLAP", nearest


def get_functional_context(svtype: str, gene_effect: str) -> str:
    """Translate genomic geometry into a cautious functional hypothesis."""
    svtype = str(svtype or "").upper()
    effect = str(gene_effect or "")

    if svtype == "INV":
        if "BREAKPOINT_IN_" in effect or "TWO_BREAKPOINTS_IN_GENE" in effect:
            return "DIRECT_TRANSCRIPT_DISRUPTION_POSSIBLE"
        if "BREAKPOINT_NEAR_GENE" in effect:
            return "REGULATORY_OR_POSITION_EFFECT_POSSIBLE_NEAR_BREAKPOINT"
        if effect in {"GENE_FULLY_SPANNED_BY_INVERSION", "GENE_INSIDE_INVERSION"}:
            return "GENE_FULLY_SPANNED_COPY_NEUTRAL_REGULATORY_3D_CONTEXT_POSSIBLE"
        return "INVERSION_FUNCTIONAL_EFFECT_UNRESOLVED"

    if svtype in {"BND", "TRA"}:
        if "BREAKPOINT_IN_" in effect or "TWO_BREAKPOINTS_IN_GENE" in effect:
            return "DIRECT_TRANSCRIPT_DISRUPTION_POSSIBLE"
        if "BREAKPOINT_NEAR_GENE" in effect:
            return "REGULATORY_OR_POSITION_EFFECT_POSSIBLE_NEAR_BREAKPOINT"
        return "BREAKEND_FUNCTIONAL_EFFECT_UNRESOLVED"

    if svtype == "DEL":
        return "COPY_LOSS_GEOMETRIC_CONTEXT"
    if svtype == "DUP":
        return "COPY_GAIN_GEOMETRIC_CONTEXT"
    if svtype == "INS":
        return "INSERTION_SITE_CONTEXT"
    return "FUNCTIONAL_EFFECT_UNRESOLVED"


def get_analysis_group(svtype: str, sv_size: float | None, gene_effect: str) -> str:
    svtype = str(svtype or "").upper()
    size = abs(sv_size) if sv_size is not None else None

    if svtype in {"DEL", "DUP"}:
        if size is not None and size >= 10_000_000:
            return "COPY_NUMBER_SV_GE_10MB"
        if size is not None and size >= 1_000_000:
            return "LARGE_COPY_NUMBER_SV"
        return "COPY_NUMBER_SV"

    if svtype == "INS":
        return "INSERTION"

    if svtype in {"INV", "BND", "TRA"}:
        if "BREAKPOINT_IN_" in gene_effect or "TWO_BREAKPOINTS_IN_GENE" in gene_effect:
            return "BREAKPOINT_IN_GENE"
        if "BREAKPOINT_NEAR_GENE" in gene_effect:
            return "BREAKPOINT_NEAR_GENE"
        if gene_effect in {"GENE_FULLY_SPANNED_BY_INVERSION", "GENE_INSIDE_INVERSION"}:
            return "INVERSION_SPANNED_GENE"
        return "BREAKPOINT_SV"

    return "OTHER_SV"
