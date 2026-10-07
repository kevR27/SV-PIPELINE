#!/usr/bin/env python3
"""Prepare a readable table of SRS structural-variant candidates.

The order reflects the biological focus of this project:

1. nuclear genes encoding mitochondrial proteins;
2. other established optic-neuropathy genes;
3. remaining genome-wide candidates;
4. mitochondrial-genome findings, kept as a secondary section.

This is a review order, not a pathogenicity classification.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


FOCUS_ORDER = {
    "NUCLEAR_MITOCHONDRIAL_PRIMARY": 1,
    "KNOWN_DISEASE_GENE": 2,
    "GENOMEWIDE_OTHER": 3,
    "MTDNA_SECONDARY": 4,
}

TECHNICAL_ORDER = {
    "MULTI_CALLER_PLUS_GRIDSS": 1,
    "MULTI_CALLER": 2,
    "SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE": 3,
    "SINGLE_CALLER": 4,
    "CNVPYTOR_DEPTH_ONLY": 5,
    "REVIEW_REQUIRED": 6,
}


def find_column(table: pd.DataFrame, possible_names: list[str]) -> str | None:
    """Return the first matching column, ignoring capitalization."""
    lower_names = {str(column).lower(): column for column in table.columns}

    for name in possible_names:
        if name in table.columns:
            return name
        if name.lower() in lower_names:
            return lower_names[name.lower()]

    return None


def positive_panel_status(values: pd.Series) -> pd.Series:
    """Recognize the panel labels used by the shared LRS/SRS scripts."""
    normalized = values.fillna("").astype(str).str.upper().str.strip()
    return normalized.isin({"PANEL_GENE", "YES", "TRUE", "1"})


def add_biological_focus(table: pd.DataFrame) -> pd.DataFrame:
    """Label each row according to the study's biological priority."""
    result = table.copy()
    result["DIAGNOSTIC_FOCUS"] = "GENOMEWIDE_OTHER"

    encoding = result.get(
        "MITOCARTA_ENCODING",
        pd.Series(".", index=result.index),
    ).fillna(".")

    panel_values = result.get(
        "PANEL_STATUS",
        result.get("panel_gene", pd.Series(".", index=result.index)),
    )

    is_panel_gene = positive_panel_status(panel_values)
    is_nuclear_mito = encoding.eq("NUCLEAR_MITOCHONDRIAL_GENE")
    is_mtdna = encoding.eq("MTDNA_ENCODED_GENE")

    result.loc[is_panel_gene, "DIAGNOSTIC_FOCUS"] = "KNOWN_DISEASE_GENE"
    result.loc[is_nuclear_mito, "DIAGNOSTIC_FOCUS"] = (
        "NUCLEAR_MITOCHONDRIAL_PRIMARY"
    )
    result.loc[is_mtdna, "DIAGNOSTIC_FOCUS"] = "MTDNA_SECONDARY"

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    candidates = pd.read_csv(
        args.input,
        sep="\t",
        dtype=str,
        low_memory=False,
    )

    gene_column = find_column(
        candidates,
        ["GENE", "GENES", "Gene", "gene", "ANNotsv_Gene"],
    )
    if gene_column is None:
        raise ValueError("Candidate input has no recognized gene column.")

    candidates = add_biological_focus(candidates)
    candidates["_FOCUS_ORDER"] = (
        candidates["DIAGNOSTIC_FOCUS"].map(FOCUS_ORDER).fillna(9)
    )

    score_column = find_column(
        candidates,
        [
            "ALLELE_RESEARCH_SCORE",
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "PHENOTYPE_SCORE",
        ],
    )
    if score_column:
        candidates["_BIOLOGICAL_SCORE"] = pd.to_numeric(
            candidates[score_column],
            errors="coerce",
        ).fillna(0)
    else:
        candidates["_BIOLOGICAL_SCORE"] = 0

    technical_support = candidates.get(
        "TECHNICAL_SUPPORT",
        pd.Series("REVIEW_REQUIRED", index=candidates.index),
    )
    candidates["_TECHNICAL_ORDER"] = (
        technical_support.map(TECHNICAL_ORDER).fillna(9)
    )

    candidates = candidates.sort_values(
        ["_FOCUS_ORDER", "_BIOLOGICAL_SCORE", "_TECHNICAL_ORDER"],
        ascending=[True, False, True],
    ).head(args.top_n)

    candidates = candidates.copy()
    candidates.insert(
        0,
        "DIAGNOSTIC_REVIEW_RANK",
        range(1, len(candidates) + 1),
    )
    candidates["DIAGNOSTIC_USE_NOTE"] = (
        "Research/diagnostic-support candidate; orthogonal validation and "
        "clinical interpretation are required."
    )

    candidates = candidates.drop(
        columns=["_FOCUS_ORDER", "_BIOLOGICAL_SCORE", "_TECHNICAL_ORDER"]
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output, sep="\t", index=False)

    print(f"[OK] diagnostic_candidates={len(candidates)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
