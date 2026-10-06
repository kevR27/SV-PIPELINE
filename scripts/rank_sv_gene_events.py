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

from ranking_common import (
    DIRECT_RELATIONSHIPS,
    MISSING,
    PROXIMAL_RELATIONSHIPS,
    event_sort_tuple,
    first_existing,
    gene_relevance,
    mechanism_inheritance_summary,
    number,
    population_summary,
    technical_summary,
)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


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
                if svtype == "INV" and lo <= gene_start and hi >= gene_end:
                    return "INVERSION_SPANS_INTACT_GENE", min_distance
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


def technical_review_status(row) -> str:
    """Flag caller-evidence concerns without excluding the discovery event."""
    flags = str(row.get("CALLER_EVIDENCE_FLAGS", "") or "").upper()
    match = str(row.get("CALLER_EVIDENCE_MATCH", "") or "").upper()

    review_tokens = (
        "RESCUED_COV_VAR",
        "LOW_GQ",
        "IMPRECISE",
        "BLACKLIST_REGION",
        "CALLER_RECORD_FAILS_CURRENT_FILTER",
    )
    if "AMBIGUOUS" in match:
        return "REVIEW_REQUIRED_AMBIGUOUS_CALLER_LINK"
    if any(token in flags for token in review_tokens):
        return "REVIEW_REQUIRED_CALLER_QC_FLAG"
    return "NO_REVIEW_FLAG_FROM_CALLER_EVIDENCE"


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


def event_bucket(svtype: str, size, relationship: str) -> str:
    large = pd.notna(size) and size >= 1_000_000
    very_large = pd.notna(size) and size >= 10_000_000

    if svtype in {"INV", "BND", "TRA"}:
        if relationship in DIRECT_RELATIONSHIPS | PROXIMAL_RELATIONSHIPS:
            return "BREAKPOINT_GENE_CANDIDATE"

        if relationship == "INVERSION_SPANS_INTACT_GENE":
            return "INVERSION_SPANNED_GENE_CONTEXT"

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
        "INVERSION_SPANS_INTACT_GENE": (
            "COPY_NEUTRAL_INVERSION_SPANNED_GENE_REGULATORY_3D_CONTEXT_REVIEW"
        ),
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
    parser.add_argument(
        "--rare-af",
        type=float,
        default=0.001,
        help="Very-rare ranking threshold for explicit needLR/gnomAD-SV AF.",
    )
    parser.add_argument(
        "--max-af",
        type=float,
        default=0.01,
        help="Maximum AF above which a rare-disease SV receives population review.",
    )
    parser.add_argument(
        "--compact-output",
        action="store_true",
        help=(
            "After using detailed AnnotSV/caller audit fields for event interpretation, "
            "drop bulky JSON and raw INFO columns that are not needed downstream. "
            "Jasmine ID-list INFO fields are retained for LongPhase matching."
        ),
    )
    args = parser.parse_args()

    if args.compact_output:
        # Inspect only the header first so very large raw evidence columns are
        # never loaded into memory. Keep AnnotSV transcript JSON temporarily
        # because it is required to classify INV/BND/gene relationships.
        header = pd.read_csv(args.input, sep="\t", nrows=0).columns.tolist()
        keep_info_columns = {
            "INFO_IDLIST",
            "INFO_IDLIST_EXT",
            "INFO_INTRASAMPLE_IDLIST",
        }
        skip_on_read = {
            "ANNOTSV_UNRESOLVED_GENE_ROWS_JSON",
            "VEP_ROWS_JSON",
            "CALLER_EVIDENCE_JSON",
        }
        usecols = [
            column
            for column in header
            if column not in skip_on_read
            and (
                not column.startswith("INFO_")
                or column in keep_info_columns
            )
        ]
        df = pd.read_csv(
            args.input,
            sep="\t",
            dtype=str,
            low_memory=False,
            usecols=usecols,
        )
        print(
            f"[INFO] compact_input loaded_columns={len(usecols)} "
            f"skipped_columns={len(header) - len(usecols)}"
        )
    else:
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

    if args.compact_output:
        # The gene-relationship calculation above is the last downstream step
        # that needs the verbose AnnotSV transcript JSON. Allele assessment has
        # already consumed caller/VEP audit JSON upstream, so remove these large
        # payloads before scoring/sorting to reduce both RAM and disk I/O.
        heavy_audit_columns = {
            "ANNOTSV_GENE_ROWS_JSON",
            "ANNOTSV_UNRESOLVED_GENE_ROWS_JSON",
            "VEP_ROWS_JSON",
            "CALLER_EVIDENCE_JSON",
        }
        # Retain only Jasmine provenance IDs used by LongPhase matching.
        keep_info_columns = {
            "INFO_IDLIST",
            "INFO_IDLIST_EXT",
            "INFO_INTRASAMPLE_IDLIST",
        }
        drop_columns = [
            column
            for column in out.columns
            if column in heavy_audit_columns
            or (
                column.startswith("INFO_")
                and column not in keep_info_columns
            )
        ]
        out = out.drop(columns=drop_columns, errors="ignore")
        print(
            f"[INFO] compact_output dropped_columns={len(drop_columns)} "
            f"retained_columns={len(out.columns)}"
        )

    # ------------------------------------------------------------------
    # Shared ranking evidence
    # ------------------------------------------------------------------
    # Gene relevance uses broad tiers so a tiny numerical difference cannot
    # outrank a much stronger SV-gene mechanism or technical signal.
    if "GENE_RELEVANCE_TIER" not in out.columns:
        semantic = (
            numeric(out["HON_SEMANTIC_SIMILARITY_NORMALIZED"]).fillna(0)
            if "HON_SEMANTIC_SIMILARITY_NORMALIZED" in out.columns
            else (
                numeric(out["PHENOTYPE_SCORE"]).fillna(0).div(10.0)
                if "PHENOTYPE_SCORE" in out.columns
                else pd.Series(0.0, index=out.index)
            )
        )
        disease = (
            numeric(out["GENE_DISEASE_EVIDENCE_SCORE"]).fillna(0)
            if "GENE_DISEASE_EVIDENCE_SCORE" in out.columns
            else pd.Series(0.0, index=out.index)
        )
        relevance_rows = [
            gene_relevance(hon, dis)
            for hon, dis in zip(semantic, disease)
        ]
        out["GENE_RELEVANCE_TIER"] = [x["tier"] for x in relevance_rows]
        out["GENE_RELEVANCE_DISPLAY_SCORE"] = [
            x["display_score"] for x in relevance_rows
        ]

    if "GENE_RELEVANCE_DISPLAY_SCORE" not in out.columns:
        out["GENE_RELEVANCE_DISPLAY_SCORE"] = (
            numeric(out.get("INTEGRATED_DISCOVERY_SCORE", 0)).fillna(0)
        )

    out["EVENT_GENE_RELEVANCE_SCORE"] = numeric(
        out["GENE_RELEVANCE_DISPLAY_SCORE"]
    ).fillna(0).round(3)

    mechanism_rows = [
        mechanism_inheritance_summary(row)
        for _, row in out.iterrows()
    ]
    out["SV_EFFECT_CLASS"] = [x["effect_class"] for x in mechanism_rows]
    out["GENE_MOI_SET"] = [
        x["moi_set"] if x["moi_set"] != "." else row.get("GENE_MOI_SET", ".")
        for x, (_, row) in zip(mechanism_rows, out.iterrows())
    ]
    out["GENE_INHERITANCE_CLASS"] = [
        x["inheritance_class"]
        if x["inheritance_class"] != "UNKNOWN"
        else row.get("GENE_INHERITANCE_CLASS", "UNKNOWN")
        for x, (_, row) in zip(mechanism_rows, out.iterrows())
    ]
    out["INHERITANCE_MECHANISM_CLASS"] = [
        x["category"] for x in mechanism_rows
    ]
    out["INHERITANCE_MECHANISM_DETAIL"] = [
        x["detail"] for x in mechanism_rows
    ]
    out["INHERITANCE_MECHANISM_PRIORITY"] = [
        x["priority"] for x in mechanism_rows
    ]
    out["CLINGEN_HI_SCORE_PARSED"] = [
        x["clingen_hi_score"] for x in mechanism_rows
    ]
    out["CLINGEN_TS_SCORE_PARSED"] = [
        x["clingen_ts_score"] for x in mechanism_rows
    ]

    technical_rows = [
        technical_summary(row)
        for _, row in out.iterrows()
    ]
    out["EVENT_TECHNICAL_TIER_V2"] = [x["tier"] for x in technical_rows]
    out["EVENT_TECHNICAL_TIER"] = out["EVENT_TECHNICAL_TIER_V2"]
    out["EVENT_MAX_CALLER_READ_SUPPORT"] = [
        x["max_read_support"] for x in technical_rows
    ]
    out["EVENT_TECHNICAL_REVIEW"] = [
        "REVIEW_REQUIRED_CALLER_EVIDENCE"
        if x["review"] == "YES"
        else "NO_REVIEW_FLAG_FROM_CALLER_EVIDENCE"
        for x in technical_rows
    ]

    population_rows = [
        population_summary(
            row,
            rare_af=args.rare_af,
            max_af=args.max_af,
        )
        for _, row in out.iterrows()
    ]
    out["EVENT_POPULATION_TIER_V2"] = [x["tier"] for x in population_rows]
    out["EVENT_POPULATION_TIER"] = out["EVENT_POPULATION_TIER_V2"]
    out["EVENT_MAX_EXPLICIT_AF"] = [x["max_af"] for x in population_rows]
    out["EVENT_MIN_EXPLICIT_AF"] = [x["min_af"] for x in population_rows]
    out["EVENT_POPULATION_AF_SOURCES"] = [
        x["sources"] for x in population_rows
    ]
    out["EVENT_POPULATION_AF_CONFLICT"] = [
        x["conflict"] for x in population_rows
    ]

    # The shared sort key is reused by downstream candidate and mitochondrial
    # rankers. Panel membership is never a score: separate ranks are produced
    # for panel and non-panel candidates instead.
    order = sorted(
        range(len(out)),
        key=lambda i: event_sort_tuple(out.iloc[i]),
    )
    out = out.iloc[order].reset_index(drop=True)
    out["EVENT_RANK_GLOBAL"] = np.arange(1, len(out) + 1)

    panel_group = (
        out["PANEL_STATUS"].fillna("UNSPECIFIED").astype(str)
        if "PANEL_STATUS" in out.columns
        else pd.Series("UNSPECIFIED", index=out.index)
    )
    out["_PANEL_GROUP_FOR_RANK"] = panel_group
    out["EVENT_RANK_WITHIN_PANEL_STATUS"] = (
        out.groupby("_PANEL_GROUP_FOR_RANK").cumcount() + 1
    )

    # Keep bucket-specific ranks because they are useful for diverse Samplot
    # selection, but use the same evidence ordering within every bucket.
    out["EVENT_RANK_WITHIN_BUCKET"] = (
        out.groupby("EVENT_REVIEW_BUCKET").cumcount() + 1
    )
    out["EVENT_RANK_WITHIN_BUCKET_PANEL_STATUS"] = (
        out.groupby(
            ["EVENT_REVIEW_BUCKET", "_PANEL_GROUP_FOR_RANK"]
        ).cumcount()
        + 1
    )
    out = out.drop(columns=["_PANEL_GROUP_FOR_RANK"])

    out["EVENT_RANKING_MODEL"] = (
        "sharedTieredRanking__geneTier_mechanismInheritance_"
        "technical_population_constraint__v5"
    )
    out["EVENT_RANKING_INTERPRETATION"] = (
        "Research prioritization only. Primary order uses broad gene-relevance "
        "tier, inheritance/mechanism compatibility, technical evidence and the "
        "maximum explicit AF across exact gnomAD-SV and provisional needLR. "
        "Panel membership is not a score. A heterozygous/direct SV in an AR "
        "gene remains a second-allele-required candidate unless biallelic or "
        "trans evidence is available. ClinGen HI/TS are mechanism-specific; "
        "dominant inheritance alone does not establish triplosensitivity. "
        "Missing population evidence is neutral, and no-match is not AF zero."
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
