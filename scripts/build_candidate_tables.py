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

from sv_gene_effects import get_analysis_group, get_sv_gene_effect, number


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


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
        return "RARE" if af <= 0.01 else "COMMON"
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
    if svtype in {"BND", "TRA"} and size is None:
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


def compact_gene_category(row: pd.Series) -> str:
    panel = text_value(row, ["PANEL_STATUS", "panel_gene"]).upper()
    anchors = numeric_value(row, ["optic_neuropathy_anchor_HPO_count", "MITO_ON_ANCHOR_HPO_COUNT"]) or 0
    disease = numeric_value(row, ["GENE_DISEASE_EVIDENCE_SCORE", "gene_disease_evidence_score"]) or 0

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
        effect, distance = get_sv_gene_effect(row_dict, gene, near_breakpoint_bp)
        svtype = text_value(source, ["SVTYPE"]).upper()
        sv_size = numeric_value(source, ["SVLEN", "SV_EVENT_SPAN_BP"])
        if sv_size is None:
            start = numeric_value(source, ["START"])
            end = numeric_value(source, ["END", "POS2"])
            if start is not None and end is not None:
                sv_size = abs(end - start)

        phenotype_score = numeric_value(source, ["PHENOTYPE_SCORE"]) or 0.0
        disease_score = numeric_value(source, ["GENE_DISEASE_EVIDENCE_SCORE"]) or 0.0
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
            "SVLEN": sv_size if sv_size is not None else ".",
            "SV_SIZE_GROUP": size_group(sv_size, svtype),
            "SV_GENE_EFFECT": effect,
            "BREAKPOINT_DISTANCE_BP": distance if distance is not None else ".",
            "SV_ANALYSIS_GROUP": get_analysis_group(svtype, sv_size, effect),
            "GENES_AFFECTED": text_value(source, ["SV_GENE_COUNT"], default="."),
            "CALLERS": text_value(source, ["CALLERS"]),
            "CALLER_COUNT": text_value(source, ["CALLER_COUNT", "SUPP"]),
            "READ_SUPPORT": text_value(source, ["CALLER_READ_SUPPORT"]),
            "CALL_SUPPORT": call_support(source),
            "NEEDLR_AF": text_value(source, ["NEEDLR_AF"]),
            "POPULATION_STATUS": population_status(source),
            "PANEL_STATUS": text_value(source, ["PANEL_STATUS"]),
            "PHENOTYPE_RELEVANCE_SCORE": round(phenotype_score, 3),
            "GENE_DISEASE_SCORE": round(disease_score, 3),
            "GENE_RELEVANCE_SCORE": round(gene_score, 3),
            "GENE_DISEASE_EVIDENCE": text_value(source, ["GENE_DISEASE_EVIDENCE_LEVEL", "GENCC"]),
            "GENE_CATEGORY": compact_gene_category(source),
            "MITOCARTA": text_value(source, ["MITOCARTA_ENCODING", "MITOCARTA_STATUS"]),
            "MITO_PATHWAY": text_value(source, ["MITOCARTA_MITOPATHWAYS"]),
            "LONGPHASE": text_value(source, ["LONGPHASE_MATCH", "LONGPHASE_STATUS"]),
            "STRAGLR": text_value(source, ["STRAGLR_MATCH", "STRAGLR_STATUS"]),
            "TLDR": text_value(source, ["TLDR_MATCH", "TLDR_STATUS"]),
            "METHYLATION_CONTEXT": text_value(source, ["METHYLATION_CONTEXT"]),
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    population_order = {
        "RARE": 0,
        "NOT_EVALUABLE_GE_10MB": 1,
        "NOT_EVALUABLE_BREAKEND": 1,
        "NO_POPULATION_MATCH": 1,
        "UNKNOWN": 2,
        "COMMON": 3,
    }
    out["_population_order"] = out["POPULATION_STATUS"].map(population_order).fillna(2)
    out["_caller_count"] = pd.to_numeric(out["CALLER_COUNT"], errors="coerce").fillna(0)

    out = out.sort_values(
        ["SV_ANALYSIS_GROUP", "GENE_RELEVANCE_SCORE", "_population_order", "_caller_count", "GENE"],
        ascending=[True, False, True, False, True],
    )
    out["RANK_IN_GROUP"] = out.groupby("SV_ANALYSIS_GROUP").cumcount() + 1
    return out.drop(columns=["_population_order", "_caller_count"])


def build_gene_table(gene_ranking: pd.DataFrame, sv_candidates: pd.DataFrame) -> pd.DataFrame:
    gene_col = first_column(gene_ranking, ["gene", "GENE"])
    if gene_col is None:
        raise ValueError("Gene ranking input has no gene column.")

    mito_by_gene = {}
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

    rows = []
    for _, source in gene_ranking.iterrows():
        gene = str(source[gene_col]).strip()
        phenotype = numeric_value(source, ["phenotype_score", "PHENOTYPE_SCORE"]) or 0.0
        disease = numeric_value(source, ["gene_disease_evidence_score", "GENE_DISEASE_EVIDENCE_SCORE"]) or 0.0
        mito, pathways = mito_by_gene.get(gene, (".", "."))

        rows.append({
            "GENE": gene,
            "PANEL_STATUS": "PANEL_GENE" if text_value(source, ["panel_gene", "PANEL_STATUS"]).upper() in {"YES", "PANEL_GENE"} else "NON_PANEL",
            "HPO_COUNT": text_value(source, ["human_HPO_count"]),
            "OPTIC_NEUROPATHY_HPO_COUNT": text_value(source, ["optic_neuropathy_anchor_HPO_count"]),
            "PHENOTYPE_RELEVANCE_SCORE": round(phenotype, 3),
            "GENE_DISEASE_EVIDENCE": text_value(source, ["gene_disease_evidence_level", "GENE_DISEASE_EVIDENCE_LEVEL"]),
            "GENE_DISEASE_SCORE": round(disease, 3),
            "GENE_RELEVANCE_SCORE": round(phenotype + disease, 3),
            "SV_COUNT": text_value(source, ["SV_count"]),
            "SV_TYPES": text_value(source, ["SV_types"]),
            "SV_LT_100KB": text_value(source, ["SV_count_lt100kb"]),
            "SV_100KB_TO_1MB": text_value(source, ["SV_count_100kb_to_1Mb"]),
            "SV_1MB_TO_10MB": text_value(source, ["SV_count_1Mb_to_10Mb"]),
            "SV_GE_10MB": text_value(source, ["SV_count_ge10Mb"]),
            "BREAKPOINT_SV_COUNT": text_value(source, ["breakpoint_defined_INV_BND_count"]),
            "MITOCARTA": mito,
            "MITO_PATHWAY": pathways,
            "GENE_CATEGORY": compact_gene_category(source),
        })

    out = pd.DataFrame(rows)
    return out.sort_values(
        ["GENE_RELEVANCE_SCORE", "GENE_DISEASE_SCORE", "OPTIC_NEUROPATHY_HPO_COUNT", "GENE"],
        ascending=[False, False, False, True],
    )


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
