#!/usr/bin/env python3
"""Check whether known positive controls appear in the candidate table."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def find_column(table: pd.DataFrame, names: list[str]) -> str | None:
    """Find a column while tolerating differences in capitalization."""
    lower_names = {str(column).lower(): column for column in table.columns}

    for name in names:
        if name in table.columns:
            return name
        if name.lower() in lower_names:
            return lower_names[name.lower()]

    return None


def truth_row_matches(
    candidates: pd.DataFrame,
    truth_row: pd.Series,
    candidate_gene_column: str,
    truth_gene_column: str,
    candidate_id_column: str | None,
    truth_id_column: str | None,
) -> pd.DataFrame:
    """Match by gene and, when supplied, by a specific SV identifier."""
    truth_gene = str(truth_row[truth_gene_column]).upper().strip()
    candidate_genes = (
        candidates[candidate_gene_column].fillna("").str.upper().str.strip()
    )
    matches = candidates[candidate_genes.eq(truth_gene)]

    if candidate_id_column and truth_id_column:
        truth_id = str(truth_row.get(truth_id_column, ".")).strip()
        if truth_id not in {"", ".", "nan"}:
            matches = matches[
                matches[candidate_id_column].astype(str).eq(truth_id)
            ]

    return matches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str)
    truth = pd.read_csv(args.truth, sep="\t", dtype=str)

    candidate_gene = find_column(candidates, ["GENE", "GENES", "gene"])
    truth_gene = find_column(truth, ["GENE", "gene"])
    if candidate_gene is None or truth_gene is None:
        raise ValueError("Candidate and truth tables must both contain GENE.")

    candidate_id = find_column(candidates, ["SV_ID", "ID"])
    truth_id = find_column(truth, ["SV_ID", "ID"])

    candidates = candidates.reset_index(drop=True)
    candidates["_REVIEW_RANK"] = candidates.index + 1

    detail_rows = []
    for _, truth_row in truth.iterrows():
        matches = truth_row_matches(
            candidates,
            truth_row,
            candidate_gene,
            truth_gene,
            candidate_id,
            truth_id,
        )

        detail_rows.append(
            {
                "GENE": truth_row[truth_gene],
                "SV_ID": truth_row.get(truth_id, ".") if truth_id else ".",
                "FOUND": "YES" if not matches.empty else "NO",
                "BEST_RANK": (
                    int(matches["_REVIEW_RANK"].min())
                    if not matches.empty
                    else "."
                ),
            }
        )

    detail = pd.DataFrame(detail_rows)
    output_tsv = Path(args.output_tsv)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output_tsv, sep="\t", index=False)

    found_count = sum(row["FOUND"] == "YES" for row in detail_rows)
    summary = {
        "truth_count": len(detail_rows),
        "found_count": found_count,
        "recall": found_count / len(detail_rows) if detail_rows else None,
    }
    Path(args.output_json).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    print(
        f"[OK] positive_controls={len(detail_rows)} "
        f"found={found_count} output={output_tsv}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
