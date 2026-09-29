#!/usr/bin/env python3
"""Create mechanism-aware ranking of individual SV-gene events.

This is a research-prioritization layer. Gene relevance, SV mechanism,
population frequency, technical support, event size and focality remain
separate fields. Large and breakpoint-defined events are retained explicitly
instead of being filtered out or treated as equivalent to small intragenic SVs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}

DIRECT_RELATIONSHIPS = {
    "WHOLE_GENE_DOSAGE_CONTEXT",
    "PARTIAL_GENE_OVERLAP",
    "INSERTION_WITHIN_TRANSCRIPT",
    "BREAKPOINT_WITHIN_TRANSCRIPT",
}

PROXIMAL_RELATIONSHIPS = {
    "INSERTION_PROXIMAL_TO_GENE",
    "BREAKPOINT_PROXIMAL_TO_GENE",
}


def first_existing(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def number(value):
    try:
        return float(value)
    except Exception:
        return None


def normalize_chrom(value) -> str | None:
    text = str(value or "").strip()
    if text.upper() in MISSING:
        return None
    return text if text.startswith("chr") else "chr" + text


def transcript_intervals(row, gene: str, json_col: str | None):
    raw = row.get(json_col, ".") if json_col else "."
    try:
        records = (
            json.loads(raw)
            if str(raw).strip().upper() not in MISSING
            else []
        )
    except Exception:
        records = []

    intervals = []
    for record in records:
        if str(record.get("Gene_name", "")).strip().upper() != gene.upper():
            continue

        start = number(record.get("Tx_start"))
        end = number(record.get("Tx_end"))
        chrom = normalize_chrom(
            record.get("SV_chrom", row.get("CHROM", "."))
        )

        if start is None or end is None or chrom is None:
            continue

        if end < start:
            start, end = end, start

        intervals.append((chrom, int(start), int(end)))

    return intervals


def breakpoint_distance(bp, tx):
    chrom, pos = bp
    tx_chrom, start, end = tx

    if chrom != tx_chrom:
        return None

    if start <= pos <= end:
        return 0

    return min(abs(pos - start), abs(pos - end))


def classify_gene_relationship(
    row,
    gene: str,
    json_col: str | None,
    tolerance: int,
):
    svtype = str(row.get("SVTYPE", "")).upper()
    chrom1 = normalize_chrom(row.get("CHROM", "."))
    chrom2 = normalize_chrom(row.get("CHR2", "."))
    start = number(row.get("START"))
    end = number(row.get("END"))
    pos2 = number(row.get("POS2"))

    if pos2 is None:
        pos2 = end

    txs = transcript_intervals(row, gene, json_col)
    if not txs:
        return "GENE_RELATIONSHIP_UNRESOLVED", np.nan

    breakpoints = []
    if chrom1 and start is not None:
        breakpoints.append((chrom1, int(start)))

    if svtype in {"INV", "BND", "TRA"} and pos2 is not None:
        breakpoints.append((chrom2 or chrom1, int(pos2)))

    distances = [
        distance
        for bp in breakpoints
        for tx in txs
        if (distance := breakpoint_distance(bp, tx)) is not None
    ]
    min_distance = min(distances) if distances else np.nan

    tx_chrom = txs[0][0]
    gene_start = min(x[1] for x in txs)
    gene_end = max(x[2] for x in txs)

    if svtype in {"INV", "BND", "TRA"}:
        if min_distance == 0:
            return "BREAKPOINT_WITHIN_TRANSCRIPT", 0

        if pd.notna(min_distance) and min_distance <= tolerance:
            return "BREAKPOINT_PROXIMAL_TO_GENE", min_distance

        if (
            chrom1 == tx_chrom
            and start is not None
            and pos2 is not None
        ):
            lo, hi = sorted((start, pos2))
            if hi >= gene_start and lo <= gene_end:
                return "INTERVAL_CONTEXT_ONLY", min_distance

        return "BREAKPOINT_RELATIONSHIP_UNRESOLVED", min_distance

    if svtype in {"DEL", "DUP", "CNV"}:
        if start is None or end is None:
            return "DOSAGE_RELATIONSHIP_UNRESOLVED", np.nan

        lo, hi = sorted((start, end))

        if lo <= gene_start and hi >= gene_end:
            return (
                "WHOLE_GENE_DOSAGE_CONTEXT",
                min(abs(lo - gene_start), abs(hi - gene_end)),
            )

        if hi >= gene_start and lo <= gene_end:
            return "PARTIAL_GENE_OVERLAP", 0

        return (
            "GENE_PROXIMAL_INTERVAL",
            min(abs(lo - gene_end), abs(hi - gene_start)),
        )

    if svtype == "INS":
        if min_distance == 0:
            return "INSERTION_WITHIN_TRANSCRIPT", 0

        if pd.notna(min_distance) and min_distance <= tolerance:
            return "INSERTION_PROXIMAL_TO_GENE", min_distance

        return "INSERTION_GENE_RELATIONSHIP_UNRESOLVED", min_distance

    return "OTHER_GENE_RELATIONSHIP", min_distance


def classify_size(
    svtype: str,
    svlen,
    chrom1,
    chrom2,
    start,
    pos2,
):
    if svtype in {"BND", "TRA"}:
        if (
            chrom1
            and chrom2
            and chrom1 == chrom2
            and start is not None
            and pos2 is not None
        ):
            return "BREAKEND_SAME_CHROM", abs(pos2 - start)

        return "BREAKEND_INTERCHROM_OR_UNRESOLVED", np.nan

    size = (
        abs(svlen)
        if svlen is not None
        else (
            abs(pos2 - start)
            if start is not None and pos2 is not None
            else np.nan
        )
    )

    if pd.isna(size):
        return "SIZE_UNKNOWN", size
    if size < 50:
        return "LT50_BP", size
    if size < 1_000:
        return "50BP_1KB", size
    if size < 10_000:
        return "1KB_10KB", size
    if size < 100_000:
        return "10KB_100KB", size
    if size < 1_000_000:
        return "100KB_1MB", size
    if size < 10_000_000:
        return "1MB_10MB", size
    return "GE_10MB", size


def focality_class(gene_count: int) -> str:
    if gene_count <= 0:
        return "NO_RESOLVED_GENE"
    if gene_count == 1:
        return "SINGLE_GENE"
    if gene_count <= 10:
        return "FEW_GENES_2_10"
    if gene_count <= 100:
        return "MULTIGENE_11_100"
    return "VERY_MULTIGENE_GT100"


def population_tier(af, status) -> str:
    status = str(status or "").upper()

    if pd.notna(af):
        return (
            "RARE_AF_LE_0.01"
            if af <= 0.01
            else "COMMON_AF_GT_0.01"
        )

    if "NOT_EVALUABLE_GE_10MB" in status:
        return "NOT_EVALUABLE_GE10MB"

    if "NOT_EVALUABLE_BND" in status:
        return "NOT_EVALUABLE_BND"

    if "NO_MATCH" in status:
        return "NO_MATCH"

    return "UNKNOWN_OR_MISSING"


def event_bucket(svtype: str, size, relationship: str) -> str:
    large = pd.notna(size) and size >= 1_000_000
    very_large = pd.notna(size) and size >= 10_000_000

    if svtype in {"INV", "BND", "TRA"}:
        if relationship in DIRECT_RELATIONSHIPS | PROXIMAL_RELATIONSHIPS:
            return "BREAKPOINT_GENE_CANDIDATE"

        if relationship == "INTERVAL_CONTEXT_ONLY":
            return (
                "LARGE_COMPLEX_INTERVAL_CONTEXT"
                if large
                else "COMPLEX_INTERVAL_CONTEXT"
            )

        return "BREAKPOINT_EVENT_UNRESOLVED"

    if svtype in {"DEL", "DUP", "CNV"}:
        if very_large:
            return "VERY_LARGE_CNV_GENE_CONTEXT"
        if large:
            return "LARGE_CNV_GENE_CANDIDATE"
        return "SMALL_MEDIUM_CNV_GENE_CANDIDATE"

    if svtype == "INS":
        return "INSERTION_GENE_CANDIDATE"

    return "OTHER_EVENT"


def interpretation_scope(relationship: str) -> str:
    mapping = {
        "WHOLE_GENE_DOSAGE_CONTEXT": "DOSAGE_CONTEXT_REVIEW",
        "PARTIAL_GENE_OVERLAP": "PARTIAL_GENE_DISRUPTION_OR_DOSAGE_REVIEW",
        "INSERTION_WITHIN_TRANSCRIPT": "INSERTION_SITE_GENE_REVIEW",
        "INSERTION_PROXIMAL_TO_GENE": "INSERTION_PROXIMAL_CONTEXT",
        "BREAKPOINT_WITHIN_TRANSCRIPT": "DIRECT_BREAKPOINT_GENE_REVIEW",
        "BREAKPOINT_PROXIMAL_TO_GENE": "BREAKPOINT_PROXIMAL_GENE_REVIEW",
        "INTERVAL_CONTEXT_ONLY": "INTERVAL_CONTEXT_ONLY_NOT_DIRECT_DISRUPTION",
    }
    return mapping.get(
        relationship,
        "RELATIONSHIP_UNRESOLVED_REVIEW",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--breakpoint-gene-tolerance",
        type=int,
        default=10_000,
    )
    args = parser.parse_args()

    df = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(
        df,
        ["GENES", "ANNotsv_Gene", "Gene", "GENE"],
    )
    json_col = first_existing(
        df,
        ["ANNOTSV_GENE_ROWS_JSON"],
    )

    if id_col is None or gene_col is None:
        raise ValueError(
            "Input requires SV_ID and overlapping-gene columns."
        )

    resolved_gene = ~(
        df[gene_col]
        .fillna(".")
        .astype(str)
        .str.upper()
        .isin(MISSING)
    )

    gene_count = (
        df.loc[resolved_gene]
        .groupby(id_col)[gene_col]
        .nunique()
        .to_dict()
    )

    out = df.copy()

    relationships = []
    distances = []
    size_classes = []
    spans = []
    buckets = []
    scopes = []

    for _, row in out.iterrows():
        gene = str(row.get(gene_col, ".")).strip()

        if gene.upper() in MISSING:
            relationship = "NO_RESOLVED_GENE"
            distance = np.nan
        else:
            relationship, distance = classify_gene_relationship(
                row,
                gene,
                json_col,
                args.breakpoint_gene_tolerance,
            )

        svtype = str(row.get("SVTYPE", ".")).upper()
        svlen = number(row.get("SVLEN"))
        chrom1 = normalize_chrom(row.get("CHROM"))
        chrom2 = normalize_chrom(row.get("CHR2"))
        start = number(row.get("START"))
        pos2 = number(row.get("POS2"))

        if pos2 is None:
            pos2 = number(row.get("END"))

        size_class, span = classify_size(
            svtype,
            svlen,
            chrom1,
            chrom2,
            start,
            pos2,
        )

        relationships.append(relationship)
        distances.append(distance)
        size_classes.append(size_class)
        spans.append(span)
        buckets.append(
            event_bucket(
                svtype,
                span,
                relationship,
            )
        )
        scopes.append(
            interpretation_scope(relationship)
        )

    out["SV_GENE_RELATIONSHIP"] = relationships
    out["BREAKPOINT_DISTANCE_TO_GENE_BP"] = distances
    out["SV_EVENT_SIZE_CLASS"] = size_classes
    out["SV_EVENT_SPAN_BP"] = spans

    out["SV_GENE_COUNT"] = (
        out[id_col]
        .map(gene_count)
        .fillna(0)
        .astype(int)
    )
    out["SV_FOCALITY_CLASS"] = (
        out["SV_GENE_COUNT"]
        .map(focality_class)
    )

    out["EVENT_REVIEW_BUCKET"] = buckets
    out["EVENT_INTERPRETATION_SCOPE"] = scopes

    out["BREAKPOINT_DEFINED_EVENT"] = (
        out["SVTYPE"]
        .astype(str)
        .str.upper()
        .isin(["INV", "BND", "TRA"])
        .map({True: "YES", False: "NO"})
    )

    out["VERY_LARGE_GE10MB"] = (
        pd.to_numeric(
            out["SV_EVENT_SPAN_BP"],
            errors="coerce",
        )
        .ge(10_000_000)
        .map({True: "YES", False: "NO"})
    )

    af = (
        numeric(out["NEEDLR_AF"])
        if "NEEDLR_AF" in out
        else pd.Series(np.nan, index=out.index)
    )
    status = (
        out["NEEDLR_STATUS"]
        if "NEEDLR_STATUS" in out
        else pd.Series(".", index=out.index)
    )

    out["EVENT_POPULATION_TIER"] = [
        population_tier(a, s)
        for a, s in zip(af, status)
    ]

    callers = (
        numeric(out["CALLER_COUNT"]).fillna(0)
        if "CALLER_COUNT" in out
        else pd.Series(0, index=out.index)
    )
    out["EVENT_TECHNICAL_TIER"] = np.where(
        callers >= 2,
        "MULTI_CALLER",
        "SINGLE_CALLER",
    )

    phenotype = (
        numeric(out["PHENOTYPE_SCORE"]).fillna(0)
        if "PHENOTYPE_SCORE" in out
        else pd.Series(0, index=out.index)
    )
    disease = (
        numeric(
            out["GENE_DISEASE_EVIDENCE_SCORE"]
        ).fillna(0)
        if "GENE_DISEASE_EVIDENCE_SCORE" in out
        else pd.Series(0, index=out.index)
    )

    # Event-level ranking intentionally excludes the gene-level SV-count
    # component. Each specific event is ranked by phenotype + curated
    # gene-disease relevance, then mechanism/population/technical context.
    out["EVENT_GENE_RELEVANCE_SCORE"] = (
        phenotype + disease
    ).round(3)

    mechanism_rank = out["SV_GENE_RELATIONSHIP"].map(
        lambda value: (
            3
            if value in DIRECT_RELATIONSHIPS
            else (
                2
                if value in PROXIMAL_RELATIONSHIPS
                else (
                    1
                    if value == "INTERVAL_CONTEXT_ONLY"
                    else 0
                )
            )
        )
    )

    population_rank = (
        out["EVENT_POPULATION_TIER"]
        .map(
            {
                "RARE_AF_LE_0.01": 3,
                "NOT_EVALUABLE_GE10MB": 2,
                "NOT_EVALUABLE_BND": 2,
                "NO_MATCH": 2,
                "UNKNOWN_OR_MISSING": 1,
                "COMMON_AF_GT_0.01": 0,
            }
        )
        .fillna(0)
    )

    panel_rank = (
        out["PANEL_STATUS"]
        .fillna("")
        .eq("PANEL_GENE")
        .astype(int)
        if "PANEL_STATUS" in out
        else pd.Series(0, index=out.index)
    )

    out["_MECHANISM_RANK"] = mechanism_rank
    out["_POPULATION_RANK"] = population_rank
    out["_CALLER_RANK"] = callers
    out["_PANEL_RANK"] = panel_rank

    out = out.sort_values(
        [
            "EVENT_REVIEW_BUCKET",
            "EVENT_GENE_RELEVANCE_SCORE",
            "_MECHANISM_RANK",
            "_POPULATION_RANK",
            "_CALLER_RANK",
            "_PANEL_RANK",
            id_col,
            gene_col,
        ],
        ascending=[
            True,
            False,
            False,
            False,
            False,
            False,
            True,
            True,
        ],
    )

    out["EVENT_RANK_WITHIN_BUCKET"] = (
        out.groupby("EVENT_REVIEW_BUCKET")
        .cumcount()
        + 1
    )

    panel_group = (
        out["PANEL_STATUS"]
        if "PANEL_STATUS" in out
        else pd.Series("UNSPECIFIED", index=out.index)
    )
    out["_PANEL_GROUP_FOR_RANK"] = panel_group.fillna("UNSPECIFIED").astype(str)
    out["EVENT_RANK_WITHIN_BUCKET_PANEL_STATUS"] = (
        out.groupby(
            ["EVENT_REVIEW_BUCKET", "_PANEL_GROUP_FOR_RANK"]
        )
        .cumcount()
        + 1
    )
    out = out.drop(columns=["_PANEL_GROUP_FOR_RANK"])

    out["EVENT_RANKING_MODEL"] = (
        "geneRelevance_phenoPlusDisease__"
        "mechanism_population_callers__v1"
    )

    out["EVENT_RANKING_INTERPRETATION"] = (
        "Research prioritization only. Event size is retained as context, "
        "not used as a pathogenicity score. INV/BND interval-only gene "
        "overlap is not treated as direct gene disruption."
    )

    out = out.drop(
        columns=[
            "_MECHANISM_RANK",
            "_POPULATION_RANK",
            "_CALLER_RANK",
            "_PANEL_RANK",
        ]
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(
        output,
        sep="\t",
        index=False,
    )

    print(
        f"[OK] ranked_sv_gene_events={len(out)} "
        f"output={output}"
    )
    print(
        out["EVENT_REVIEW_BUCKET"]
        .value_counts()
        .to_string()
    )


if __name__ == "__main__":
    main()
