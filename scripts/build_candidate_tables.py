#!/usr/bin/env python3
"""Build concise, human-readable gene and SV-gene candidate tables.

Existing pipeline outputs are not replaced. These tables are an additional
interpretation layer for thesis review.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ranking_common import MISSING, event_sort_tuple
from sv_gene_effects import (
    get_analysis_group,
    get_functional_context,
    get_sv_gene_effect,
    number,
)


def first_column(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(col).lower(): col for col in df.columns}
    for name in names:
        if name in df.columns:
            return name
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def text_value(row: pd.Series, names: list[str], default: str = ".") -> str:
    for name in names:
        if name in row.index:
            value = str(row[name]).strip()
            if value.upper() not in MISSING:
                return value
    return default


def numeric_value(row: pd.Series, names: list[str]) -> float | None:
    value = text_value(row, names)
    return number(value)


def population_status(row: pd.Series) -> str:
    af = numeric_value(row, ["NEEDLR_AF"])
    status = text_value(row, ["NEEDLR_STATUS"]).upper()

    if af is not None:
        return (
            "PROVISIONAL_LOW_FREQUENCY"
            if af <= 0.01
            else "PROVISIONAL_HIGH_FREQUENCY"
        )
    if "NOT_EVALUABLE_GE_10MB" in status or "NOT_EVALUABLE_GE10MB" in status:
        return "NOT_EVALUABLE_GE_10MB"
    if "NOT_EVALUABLE_BND" in status:
        return "NOT_EVALUABLE_BREAKEND"
    if "NO_MATCH" in status:
        return "NO_POPULATION_MATCH"
    return "UNKNOWN"


def call_support(row: pd.Series) -> str:
    count = numeric_value(row, ["CALLER_COUNT", "SUPP"])
    if count is None:
        return "UNKNOWN"
    return "MULTI_CALLER" if count >= 2 else "SINGLE_CALLER"


def size_group(size: float | None, svtype: str) -> str:
    if svtype in {"BND", "TRA"}:
        return "BREAKEND"
    if size is None:
        return "UNKNOWN"
    size = abs(size)
    if size < 100_000:
        return "LT_100KB"
    if size < 1_000_000:
        return "100KB_TO_1MB"
    if size < 10_000_000:
        return "1MB_TO_10MB"
    return "GE_10MB"


def relationship_to_effect(relationship: str, svtype: str) -> str:
    """Convert compact event relationship to the legacy readable effect label.

    The compact postprocess table intentionally no longer retains
    ANNOTSV_GENE_ROWS_JSON. rank_sv_gene_events.py has already converted that
    geometry into SV_GENE_RELATIONSHIP, which is the authoritative downstream
    field. These labels preserve compatibility with the existing thesis tables
    and R plots without reconstructing gene geometry from deleted JSON.
    """
    relationship = str(relationship or "").upper()
    svtype = str(svtype or "").upper()

    if relationship == "WHOLE_GENE_DOSAGE_CONTEXT":
        if svtype == "DEL":
            return "WHOLE_GENE_DELETION"
        if svtype == "DUP":
            return "WHOLE_GENE_DUPLICATION"
        return "WHOLE_GENE_DOSAGE_CONTEXT"

    if relationship == "PARTIAL_GENE_OVERLAP":
        if svtype == "DEL":
            return "PARTIAL_GENE_DELETION"
        if svtype == "DUP":
            return "PARTIAL_GENE_DUPLICATION"
        return "PARTIAL_GENE_OVERLAP"

    if relationship == "INSERTION_WITHIN_TRANSCRIPT":
        return "INSERTION_IN_TRANSCRIPT"
    if relationship == "INSERTION_PROXIMAL_TO_GENE":
        return "INSERTION_NEAR_GENE"

    if relationship == "BREAKPOINT_WITHIN_TRANSCRIPT":
        prefix = "INVERSION" if svtype == "INV" else "BREAKEND"
        return f"{prefix}_BREAKPOINT_IN_TRANSCRIPT"

    if relationship == "BREAKPOINT_PROXIMAL_TO_GENE":
        prefix = "INVERSION" if svtype == "INV" else "BREAKEND"
        return f"{prefix}_BREAKPOINT_NEAR_GENE"

    if relationship == "INVERSION_SPANS_INTACT_GENE":
        return "GENE_FULLY_SPANNED_BY_INVERSION"

    if relationship == "INTERVAL_CONTEXT_ONLY":
        return (
            "INVERSION_INTERVAL_CONTEXT_ONLY"
            if svtype == "INV"
            else "BREAKEND_INTERVAL_CONTEXT_ONLY"
        )

    if relationship == "GENE_PROXIMAL_INTERVAL":
        if svtype == "DEL":
            return "DELETION_NEAR_GENE"
        if svtype == "DUP":
            return "DUPLICATION_NEAR_GENE"
        return f"{svtype or 'SV'}_NEAR_GENE"

    return relationship if relationship not in MISSING else "GENE_EFFECT_UNRESOLVED"


def compact_gene_category(row: pd.Series) -> str:
    existing = text_value(
        row,
        ["CANDIDATE_CLASS", "candidate_group", "classification"],
    ).upper()

    existing_map = {
        "PANEL_GENE": "OPTIC_NEUROPATHY_PANEL",
        "NONPANEL_HPO_AND_DISEASE_EVIDENCE": "NON_PANEL_PHENOTYPE_AND_DISEASE",
        "NONPANEL_HPO_OVERLAP": "NON_PANEL_PHENOTYPE",
        "NONPANEL_HUMAN_DISEASE_GENE": "NON_PANEL_DISEASE_GENE",
        "OTHER_NONPANEL_CANDIDATE": "OTHER_NON_PANEL",
    }
    if existing in existing_map:
        return existing_map[existing]

    panel = text_value(row, ["PANEL_STATUS", "panel_gene"]).upper()
    anchors = numeric_value(
        row,
        ["optic_neuropathy_anchor_HPO_count", "MITO_ON_ANCHOR_HPO_COUNT"],
    ) or 0
    disease = numeric_value(
        row,
        ["GENE_DISEASE_EVIDENCE_SCORE", "gene_disease_evidence_score"],
    ) or 0

    if panel in {"PANEL_GENE", "YES"}:
        return "OPTIC_NEUROPATHY_PANEL"
    if anchors > 0 and disease > 0:
        return "NON_PANEL_PHENOTYPE_AND_DISEASE"
    if anchors > 0:
        return "NON_PANEL_PHENOTYPE"
    if disease > 0:
        return "NON_PANEL_DISEASE_GENE"
    return "OTHER_NON_PANEL"


def build_sv_table(events: pd.DataFrame, near_breakpoint_bp: int) -> pd.DataFrame:
    gene_col = first_column(events, ["GENES", "ANNotsv_Gene", "GENE", "Gene"])
    if gene_col is None:
        raise ValueError("SV-gene input has no gene column.")

    rows = []
    for _, source in events.iterrows():
        gene = str(source[gene_col]).strip()
        if gene.upper() in MISSING:
            continue

        row_dict = source.to_dict()
        svtype = text_value(source, ["SVTYPE"]).upper()

        # Compact postprocess outputs already contain the mechanism-aware
        # relationship calculated before ANNOTSV_GENE_ROWS_JSON is removed.
        relationship = text_value(source, ["SV_GENE_RELATIONSHIP"])
        if relationship != ".":
            effect = relationship_to_effect(relationship, svtype)
            distance = numeric_value(
                source,
                ["BREAKPOINT_DISTANCE_TO_GENE_BP"],
            )
        else:
            # Legacy completed runs may still carry AnnotSV transcript JSON.
            effect, distance = get_sv_gene_effect(
                row_dict,
                gene,
                near_breakpoint_bp,
            )
            relationship = "."
        original_svlen = numeric_value(source, ["SVLEN"])
        sv_size = numeric_value(source, ["SV_EVENT_SPAN_BP"])
        if sv_size is None:
            sv_size = original_svlen
        if sv_size is None:
            start = numeric_value(source, ["START"])
            end = numeric_value(source, ["END", "POS2"])
            if start is not None and end is not None:
                sv_size = abs(end - start)

        phenotype_score = numeric_value(
            source,
            [
                "HON_SEMANTIC_SCORE_0_10",
                "PHENOTYPE_SCORE",
                "phenotype_score",
            ],
        ) or 0.0
        disease_score = numeric_value(
            source,
            ["GENE_DISEASE_EVIDENCE_SCORE"],
        ) or 0.0
        gene_score = numeric_value(
            source,
            [
                "GENE_RELEVANCE_DISPLAY_SCORE",
                "EVENT_GENE_RELEVANCE_SCORE",
                "INTEGRATED_DISCOVERY_SCORE",
            ],
        )
        if gene_score is None:
            gene_score = phenotype_score + disease_score

        rows.append({
            "SV_ID": text_value(source, ["SV_ID", "ID"]),
            "GENE": gene,
            "CHROM": text_value(source, ["CHROM"]),
            "START": text_value(source, ["START"]),
            "END": text_value(source, ["END"]),
            "CHR2": text_value(source, ["CHR2"]),
            "POS2": text_value(source, ["POS2"]),
            "SVTYPE": svtype,
            "SVLEN": original_svlen if original_svlen is not None else ".",
            "SV_SPAN_BP": sv_size if sv_size is not None else ".",
            "SV_SIZE_GROUP": size_group(sv_size, svtype),
            "SV_GENE_RELATIONSHIP": relationship,
            "EVENT_INTERPRETATION_SCOPE": text_value(
                source,
                ["EVENT_INTERPRETATION_SCOPE"],
            ),
            "EVENT_REVIEW_BUCKET": text_value(
                source,
                ["EVENT_REVIEW_BUCKET"],
            ),
            "SV_GENE_EFFECT": effect,
            "SV_FUNCTIONAL_CONTEXT": get_functional_context(svtype, effect),
            "BREAKPOINT_DISTANCE_BP": distance if distance is not None else ".",
            "SV_ANALYSIS_GROUP": get_analysis_group(svtype, sv_size, effect),
            "GENES_AFFECTED": text_value(source, ["SV_GENE_COUNT"], default="."),
            "CALLERS": text_value(source, ["CALLERS"]),
            "CALLER_COUNT": text_value(source, ["CALLER_COUNT", "SUPP"]),
            "READ_SUPPORT": text_value(source, ["CALLER_READ_SUPPORT"]),
            "CALL_SUPPORT": call_support(source),
            "CALLER_EVIDENCE_FLAGS": text_value(source, ["CALLER_EVIDENCE_FLAGS"]),
            "CALLER_EVIDENCE_MATCH": text_value(source, ["CALLER_EVIDENCE_MATCH"]),
            "TECHNICAL_EVIDENCE_SCOPE": (
                "CALLER_COUNT_AND_READ_SUPPORT_ARE_SUPPORTING_EVIDENCE; "
                "REVIEW_CALLER_FLAGS_AND_READ_LEVEL_SIGNAL"
            ),
            "NEEDLR_AF": text_value(source, ["NEEDLR_AF"]),
            "POPULATION_STATUS": population_status(source),
            "GNOMAD_SV_OVERLAP": text_value(source, ["SV_DB_GNOMAD_OVERLAP"]),
            "GNOMAD_SV_EXACT_MATCH": text_value(source, ["GNOMAD_SV_EXACT_MATCH"]),
            "GNOMAD_SV_ID": text_value(source, ["GNOMAD_SV_ID"]),
            "GNOMAD_SV_AF": text_value(source, ["GNOMAD_SV_AF"]),
            "GNOMAD_SV_AC": text_value(source, ["GNOMAD_SV_AC"]),
            "GNOMAD_SV_AN": text_value(source, ["GNOMAD_SV_AN"]),
            "GNOMAD_SV_FILTER": text_value(source, ["GNOMAD_SV_FILTER"]),
            "GNOMAD_SV_MATCH_SCOPE": text_value(source, ["GNOMAD_SV_MATCH_SCOPE"]),
            "ANNOTSV_BENIGN_AFMAX": text_value(source, ["SV_BENIGN_DB_AFMAX"]),
            "ANNOTSV_BENIGN_DB_SOURCE": text_value(source, ["SV_BENIGN_DB_SOURCE"]),
            # Backward-compatible aliases. AFmax can reflect multiple benign
            # resources; it is not necessarily a gnomAD-specific exact-allele AF.
            "BENIGN_DB_AFMAX": text_value(source, ["SV_BENIGN_DB_AFMAX"]),
            "BENIGN_DB_SOURCE": text_value(source, ["SV_BENIGN_DB_SOURCE"]),
            "ANNOTSV_CLASSIFICATION": text_value(
                source,
                ["ANNOTSV_GENERAL_CLASSIFICATION", "ANNotsv_Classification"],
            ),
            "ANNOTSV_CLASSIFICATION_SCOPE": text_value(
                source,
                ["ANNOTSV_CLASSIFICATION_SCOPE"],
            ),
            "VEP_MATCH_STATUS": text_value(source, ["VEP_MATCH_STATUS"]),
            "VEP_MATCH_METHOD": text_value(source, ["VEP_MATCH_METHOD"]),
            "VEP_GENE_IDS": text_value(source, ["VEP_GENE_IDS"]),
            "VEP_TRANSCRIPT_COUNT": text_value(source, ["VEP_TRANSCRIPT_COUNT"]),
            "VEP_TRANSCRIPTS": text_value(source, ["VEP_TRANSCRIPTS"]),
            "VEP_WHOLE_TRANSCRIPT_COUNT": text_value(source, ["VEP_WHOLE_TRANSCRIPT_COUNT"]),
            "VEP_GENE_TRANSCRIPT_SCOPE": text_value(source, ["VEP_GENE_TRANSCRIPT_SCOPE"]),
            "VEP_CONSEQUENCES": text_value(source, ["VEP_CONSEQUENCES"]),
            "VEP_IMPACTS": text_value(source, ["VEP_IMPACTS"]),
            "VEP_BIOTYPES": text_value(source, ["VEP_BIOTYPES"]),
            "VEP_EXON": text_value(source, ["VEP_EXON"]),
            "VEP_INTRON": text_value(source, ["VEP_INTRON"]),
            "VEP_CANONICAL_TRANSCRIPTS": text_value(source, ["VEP_CANONICAL_TRANSCRIPTS"]),
            "VEP_PICK_TRANSCRIPTS": text_value(source, ["VEP_PICK_TRANSCRIPTS"]),
            "VEP_OVERLAP_BP_MAX": text_value(source, ["VEP_OVERLAP_BP_MAX"]),
            "VEP_OVERLAP_PC_MAX": text_value(source, ["VEP_OVERLAP_PC_MAX"]),
            "VEP_TRANSCRIPT_REGION_CLASS": text_value(source, ["VEP_TRANSCRIPT_REGION_CLASS"]),
            "VEP_STRUCTURAL_EFFECT": text_value(source, ["VEP_STRUCTURAL_EFFECT"]),
            "VEP_TRANSCRIPT_CONTEXT_STATUS": text_value(
                source,
                ["ALLELE_VEP_TRANSCRIPT_CONTEXT_STATUS"],
            ),
            "VEP_TRANSCRIPT_CONTEXT_DETAIL": text_value(
                source,
                ["ALLELE_VEP_TRANSCRIPT_CONTEXT_DETAIL"],
            ),
            "ACMG_CNV_CLASS": text_value(source, ["ACMG_CNV_CLASS"]),
            "DOSAGE_RELEVANCE": text_value(source, ["DOSAGE_RELEVANCE"]),
            "PANEL_STATUS": text_value(source, ["PANEL_STATUS"]),
            "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED"],
            ),
            "HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED"],
            ),
            "HON_SEMANTIC_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_SEMANTIC_SIMILARITY_NORMALIZED"],
            ),
            "HON_SEMANTIC_METHOD": text_value(
                source,
                ["HON_SEMANTIC_METHOD"],
                default="ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF",
            ),
            "GENE_RELEVANCE_TIER": text_value(
                source,
                ["GENE_RELEVANCE_TIER"],
                default="LIMITED",
            ),
            "GENE_RELEVANCE_DISPLAY_SCORE": round(gene_score, 3),
            "GENE_INHERITANCE_CLASS": text_value(
                source,
                ["GENE_INHERITANCE_CLASS"],
                default="UNKNOWN",
            ),
            "GENE_MOI_SET": text_value(source, ["GENE_MOI_SET"]),
            "INHERITANCE_MECHANISM_CLASS": text_value(
                source,
                ["INHERITANCE_MECHANISM_CLASS"],
                default="UNRESOLVED",
            ),
            "INHERITANCE_MECHANISM_DETAIL": text_value(
                source,
                ["INHERITANCE_MECHANISM_DETAIL"],
            ),
            "EVENT_TECHNICAL_TIER_V2": text_value(
                source,
                ["EVENT_TECHNICAL_TIER_V2", "EVENT_TECHNICAL_TIER"],
                default="REVIEW",
            ),
            "EVENT_MAX_CALLER_READ_SUPPORT": text_value(
                source,
                ["EVENT_MAX_CALLER_READ_SUPPORT"],
            ),
            "EVENT_POPULATION_TIER_V2": text_value(
                source,
                ["EVENT_POPULATION_TIER_V2", "EVENT_POPULATION_TIER"],
                default="UNKNOWN",
            ),
            "EVENT_MAX_EXPLICIT_AF": text_value(
                source,
                ["EVENT_MAX_EXPLICIT_AF"],
            ),
            "EVENT_POPULATION_AF_CONFLICT": text_value(
                source,
                ["EVENT_POPULATION_AF_CONFLICT"],
            ),
            "CLINGEN_HI": text_value(source, ["CLINGEN_HI"]),
            "CLINGEN_TS": text_value(source, ["CLINGEN_TS"]),
            "GNOMAD_LOEUF": text_value(source, ["GNOMAD_LOEUF"]),
            "EVENT_RANK_GLOBAL": text_value(source, ["EVENT_RANK_GLOBAL"]),
            "EVENT_RANK_WITHIN_PANEL_STATUS": text_value(
                source,
                ["EVENT_RANK_WITHIN_PANEL_STATUS"],
            ),
            "HON_CONTEXT_SCORE": round(phenotype_score, 3),
            "PHENOTYPE_RELEVANCE_SCORE": round(phenotype_score, 3),
            "PHENOTYPE_SCORE_SCOPE": (
                "GENERIC_HON_ASYMMETRIC_RESNIK_QUERY_COVERAGE_"
                "CORE_PLUS_CAPPED_MITO_CONTEXT_NOT_PATIENT_SPECIFIC"
            ),
            "GENE_DISEASE_SCORE": round(disease_score, 3),
            "GENE_RELEVANCE_SCORE": round(gene_score, 3),
            "GENE_RELEVANCE_SCOPE": (
                "TIERED_HON_SEMANTIC_PLUS_CURATED_DISEASE_"
                "PRIMARY_ORDER_USES_SHARED_EVENT_RANKING"
            ),
            "GENE_DISEASE_EVIDENCE": text_value(source, ["GENE_DISEASE_EVIDENCE_LEVEL", "GENCC"]),
            "GENE_CATEGORY": compact_gene_category(source),
            "MITOCARTA": text_value(source, ["MITOCARTA_ENCODING", "MITOCARTA_STATUS"]),
            "MITO_PATHWAY": text_value(source, ["MITOCARTA_MITOPATHWAYS"]),
            "LONGPHASE": text_value(source, ["LONGPHASE_MATCH", "LONGPHASE_STATUS"]),
            "STRAGLR": text_value(source, ["STRAGLR_MATCH", "STRAGLR_STATUS"]),
            "STRAGLR_CONTEXT": text_value(source, ["STRAGLR_CONTEXT"]),
            "STRAGLR_CONTEXT_LOCI": text_value(source, ["STRAGLR_CONTEXT_LOCI"]),
            "TLDR": text_value(source, ["TLDR_MATCH", "TLDR_STATUS"]),
            "LONGPHASE_GT": text_value(source, ["LONGPHASE_GT"]),
            "LONGPHASE_PS": text_value(source, ["LONGPHASE_PS"]),
            "LONGPHASE_PHASED": text_value(source, ["LONGPHASE_PHASED"]),
            "BND_ORIENTATION": text_value(source, ["BND_ORIENTATION"]),
            "METHYLATION_CONTEXT": text_value(source, ["METHYLATION_CONTEXT"]),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    # Reuse the exact same ordering logic as rank_sv_gene_events.py.
    order = sorted(
        range(len(out)),
        key=lambda i: event_sort_tuple(out.iloc[i]),
    )
    out = out.iloc[order].reset_index(drop=True)
    out["RANK_IN_GROUP"] = (
        out.groupby("SV_ANALYSIS_GROUP").cumcount() + 1
    )
    out["RANK_WITHIN_PANEL_STATUS"] = (
        out.groupby("PANEL_STATUS").cumcount() + 1
    )
    return out


def build_gene_table(gene_ranking: pd.DataFrame, sv_candidates: pd.DataFrame) -> pd.DataFrame:
    gene_col = first_column(gene_ranking, ["gene", "GENE"])
    if gene_col is None:
        raise ValueError("Gene ranking input has no gene column.")

    mito_by_gene = {}
    vep_by_gene = {}
    if not sv_candidates.empty:
        for gene, group in sv_candidates.groupby("GENE"):
            mito = sorted({
                str(value)
                for value in group["MITOCARTA"]
                if str(value).upper() not in MISSING
            })
            pathways = sorted({
                item.strip()
                for value in group["MITO_PATHWAY"]
                for item in str(value).split(";")
                if item.strip().upper() not in MISSING
            })
            mito_by_gene[gene] = (
                ";".join(mito) if mito else ".",
                ";".join(pathways) if pathways else ".",
            )

            effects = sorted({
                item.strip()
                for value in group["VEP_STRUCTURAL_EFFECT"]
                for item in str(value).split(";")
                if item.strip().upper() not in MISSING
            })
            regions = sorted({
                item.strip()
                for value in group["VEP_TRANSCRIPT_REGION_CLASS"]
                for item in str(value).split(";")
                if item.strip().upper() not in MISSING
            })
            vep_by_gene[gene] = {
                "effects": ";".join(effects) if effects else ".",
                "regions": ";".join(regions) if regions else ".",
                "whole": int(group["VEP_STRUCTURAL_EFFECT"].astype(str).str.contains(
                    "WHOLE_TRANSCRIPT_LOSS|WHOLE_TRANSCRIPT_GAIN", regex=True, na=False
                ).sum()),
                "exonic": int(group["VEP_TRANSCRIPT_REGION_CLASS"].astype(str).str.contains(
                    "EXONIC_OR_SPLICE", regex=False, na=False
                ).sum()),
                "intronic": int(group["VEP_TRANSCRIPT_REGION_CLASS"].astype(str).str.contains(
                    "INTRONIC", regex=False, na=False
                ).sum()),
            }

    rows = []
    for _, source in gene_ranking.iterrows():
        gene = str(source[gene_col]).strip()
        phenotype = numeric_value(
            source,
            ["hon_context_score", "HON_CONTEXT_SCORE", "phenotype_score", "PHENOTYPE_SCORE"],
        ) or 0.0
        disease = numeric_value(
            source,
            ["gene_disease_evidence_score", "GENE_DISEASE_EVIDENCE_SCORE"],
        ) or 0.0
        relevance_display = numeric_value(
            source,
            [
                "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
                "GENE_RELEVANCE_DISPLAY_SCORE",
                "integrated_discovery_score",
            ],
        )
        if relevance_display is None:
            relevance_display = phenotype + disease
        mito, pathways = mito_by_gene.get(gene, (".", "."))
        vep = vep_by_gene.get(
            gene,
            {"effects": ".", "regions": ".", "whole": 0, "exonic": 0, "intronic": 0},
        )

        rows.append({
            "GENE": gene,
            "PANEL_STATUS": (
                "PANEL_GENE"
                if text_value(
                    source,
                    ["panel_gene", "PANEL_STATUS"],
                ).upper() in {"YES", "PANEL_GENE"}
                else "NONPANEL_GENE"
            ),
            "GENE_RELEVANCE_TIER": text_value(
                source,
                ["GENE_RELEVANCE_TIER"],
                default="LIMITED",
            ),
            "GENE_RELEVANCE_DISPLAY_SCORE": text_value(
                source,
                ["GENE_RELEVANCE_DISPLAY_SCORE", "integrated_discovery_score"],
                default="0",
            ),
            "GENE_INHERITANCE_CLASS": text_value(
                source,
                ["GENE_INHERITANCE_CLASS"],
                default="UNKNOWN",
            ),
            "GENE_MOI_SET": text_value(source, ["GENE_MOI_SET"]),
            "GENE_RANK_WITHIN_PANEL_STATUS": text_value(
                source,
                ["GENE_RANK_WITHIN_PANEL_STATUS"],
            ),
            "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED"],
            ),
            "HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED"],
            ),
            "HON_SEMANTIC_SIMILARITY_NORMALIZED": text_value(
                source,
                ["HON_SEMANTIC_SIMILARITY_NORMALIZED"],
            ),
            "HON_SEMANTIC_METHOD": text_value(
                source,
                ["HON_SEMANTIC_METHOD"],
                default="ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF",
            ),
            "HPO_COUNT": text_value(source, ["human_HPO_count"]),
            "OPTIC_NEUROPATHY_HPO_COUNT": text_value(source, ["optic_neuropathy_anchor_HPO_count"]),
            "HON_CONTEXT_SCORE": round(phenotype, 3),
            "PHENOTYPE_RELEVANCE_SCORE": round(phenotype, 3),
            "PHENOTYPE_SCORE_SCOPE": "GENERIC_HON_ANCHOR_CONTEXT_NOT_PATIENT_SPECIFIC",
            "GENE_DISEASE_EVIDENCE": text_value(source, ["gene_disease_evidence_level", "GENE_DISEASE_EVIDENCE_LEVEL"]),
            "GENE_DISEASE_SCORE": round(disease, 3),
            "GENE_RELEVANCE_SCORE": round(relevance_display, 3),
            "GENE_RELEVANCE_SCOPE": (
                "BALANCED_HON_SEMANTIC_AND_CURATED_DISEASE_DISPLAY_SCORE_"
                "PRIMARY_ORDER_USES_RELEVANCE_TIER_AND_BEST_EVENT"
            ),
            "SV_COUNT": text_value(source, ["SV_count"]),
            "SV_TYPES": text_value(source, ["SV_types"]),
            "SV_LT_100KB": text_value(source, ["SV_count_lt100kb"]),
            "SV_100KB_TO_1MB": text_value(source, ["SV_count_100kb_to_1Mb"]),
            "SV_1MB_TO_10MB": text_value(source, ["SV_count_1Mb_to_10Mb"]),
            "SV_GE_10MB": text_value(source, ["SV_count_ge10Mb"]),
            "BREAKPOINT_SV_COUNT": text_value(source, ["breakpoint_defined_INV_BND_count"]),
            "MITOCARTA": mito,
            "MITO_PATHWAY": pathways,
            "VEP_STRUCTURAL_EFFECTS": vep["effects"],
            "VEP_TRANSCRIPT_REGION_CLASSES": vep["regions"],
            "VEP_WHOLE_TRANSCRIPT_EVENT_COUNT": vep["whole"],
            "VEP_EXONIC_OR_SPLICE_EVENT_COUNT": vep["exonic"],
            "VEP_INTRONIC_EVENT_COUNT": vep["intronic"],
            "GENE_CATEGORY": compact_gene_category(source),
        })

    out = pd.DataFrame(rows)
    tier_rank = {
        "HIGH": 3,
        "MODERATE": 2,
        "SUPPORTING": 1,
        "LIMITED": 0,
    }
    out["_tier_rank"] = (
        out["GENE_RELEVANCE_TIER"]
        .fillna("LIMITED")
        .astype(str)
        .str.upper()
        .map(tier_rank)
        .fillna(0)
    )
    out["_display_score"] = pd.to_numeric(
        out["GENE_RELEVANCE_DISPLAY_SCORE"],
        errors="coerce",
    ).fillna(0)
    out = out.sort_values(
        [
            "PANEL_STATUS",
            "_tier_rank",
            "_display_score",
            "GENE_DISEASE_SCORE",
            "GENE",
        ],
        ascending=[True, False, False, False, True],
    )
    out["FINAL_GENE_RANK_WITHIN_PANEL_STATUS"] = (
        out.groupby("PANEL_STATUS").cumcount() + 1
    )
    return out.drop(columns=["_tier_rank", "_display_score"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-ranking", required=True)
    parser.add_argument("--sv-gene-events", required=True)
    parser.add_argument("--gene-output", required=True)
    parser.add_argument("--sv-output", required=True)
    parser.add_argument("--breakpoint-near-bp", type=int, default=10_000)
    args = parser.parse_args()

    genes = pd.read_csv(args.gene_ranking, sep="\t", dtype=str, low_memory=False)
    events = pd.read_csv(args.sv_gene_events, sep="\t", dtype=str, low_memory=False)

    sv_candidates = build_sv_table(events, args.breakpoint_near_bp)
    gene_candidates = build_gene_table(genes, sv_candidates)

    gene_path = Path(args.gene_output)
    sv_path = Path(args.sv_output)
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    sv_path.parent.mkdir(parents=True, exist_ok=True)

    gene_candidates.to_csv(gene_path, sep="\t", index=False)
    sv_candidates.to_csv(sv_path, sep="\t", index=False)

    print(f"[OK] gene_candidates={len(gene_candidates)} output={gene_path}")
    print(f"[OK] sv_gene_candidates={len(sv_candidates)} output={sv_path}")


if __name__ == "__main__":
    main()
