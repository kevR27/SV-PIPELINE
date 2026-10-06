#!/usr/bin/env python3
"""Benchmark the final SV-gene ranking against a predefined truth set.

The truth file must contain GENE and may contain SV_ID. When SV_ID is supplied,
that exact SV-gene pair is evaluated. Otherwise the best-ranked event for the
truth gene is used.

This utility evaluates ranking performance; it does not train or optimize
weights.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def first_existing(df, names):
    lower = {str(c).lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ranking", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--detail-output", required=True)
    p.add_argument("--summary-output", required=True)
    p.add_argument(
        "--top-k",
        default="1,5,10,20,50",
        help="Comma-separated top-k cutoffs.",
    )
    args = p.parse_args()

    ranking = pd.read_csv(args.ranking, sep="\t", dtype=str, low_memory=False)
    truth = pd.read_csv(args.truth, sep="\t", dtype=str, low_memory=False)

    gene_col = first_existing(ranking, ["GENE", "GENES", "gene"])
    sv_col = first_existing(ranking, ["SV_ID", "ID"])
    global_rank_col = first_existing(
        ranking,
        ["FINAL_EVENT_RANK_GLOBAL", "EVENT_RANK_GLOBAL"],
    )
    panel_rank_col = first_existing(
        ranking,
        [
            "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS",
            "EVENT_RANK_WITHIN_PANEL_STATUS",
        ],
    )
    panel_col = first_existing(ranking, ["PANEL_STATUS"])

    if gene_col is None or sv_col is None or global_rank_col is None:
        raise ValueError(
            "Ranking requires GENE, SV_ID and final/global event rank columns."
        )
    if "GENE" not in truth.columns:
        raise ValueError("Truth table requires a GENE column.")

    ranking = ranking.copy()
    ranking["_GLOBAL_RANK"] = pd.to_numeric(
        ranking[global_rank_col],
        errors="coerce",
    )
    ranking["_PANEL_RANK"] = (
        pd.to_numeric(ranking[panel_rank_col], errors="coerce")
        if panel_rank_col
        else np.nan
    )
    ranking["_GENE"] = ranking[gene_col].fillna(".").astype(str).str.upper()
    ranking["_SV"] = ranking[sv_col].fillna(".").astype(str)

    total_ranked = int(ranking["_GLOBAL_RANK"].notna().sum())
    rows = []

    for _, target in truth.iterrows():
        gene = str(target["GENE"]).strip().upper()
        target_sv = str(target.get("SV_ID", "")).strip()

        candidates = ranking[ranking["_GENE"].eq(gene)].copy()
        match_scope = "GENE_BEST_EVENT"

        if target_sv and target_sv not in {".", "NA", "N/A", "nan", "None"}:
            exact = candidates[candidates["_SV"].eq(target_sv)].copy()
            if not exact.empty:
                candidates = exact
                match_scope = "EXACT_SV_GENE"
            else:
                candidates = ranking.iloc[0:0].copy()
                match_scope = "EXACT_SV_GENE_NOT_FOUND"

        if candidates.empty:
            rows.append(
                {
                    "GENE": gene,
                    "TRUTH_SV_ID": target_sv or ".",
                    "MATCH_SCOPE": match_scope,
                    "FOUND": "NO",
                    "RANKED_SV_ID": ".",
                    "GLOBAL_RANK": ".",
                    "PANEL_STATUS": ".",
                    "WITHIN_PANEL_RANK": ".",
                    "GLOBAL_RANK_PERCENTILE": ".",
                }
            )
            continue

        best = candidates.sort_values("_GLOBAL_RANK").iloc[0]
        rank = float(best["_GLOBAL_RANK"])
        percentile = (
            100.0 * rank / total_ranked
            if total_ranked > 0
            else np.nan
        )
        rows.append(
            {
                "GENE": gene,
                "TRUTH_SV_ID": target_sv or ".",
                "MATCH_SCOPE": match_scope,
                "FOUND": "YES",
                "RANKED_SV_ID": best["_SV"],
                "GLOBAL_RANK": int(rank),
                "PANEL_STATUS": (
                    str(best.get(panel_col, "."))
                    if panel_col
                    else "."
                ),
                "WITHIN_PANEL_RANK": (
                    int(best["_PANEL_RANK"])
                    if pd.notna(best["_PANEL_RANK"])
                    else "."
                ),
                "GLOBAL_RANK_PERCENTILE": round(percentile, 4),
            }
        )

    detail = pd.DataFrame(rows)
    cutoffs = sorted({
        int(value)
        for value in args.top_k.split(",")
        if value.strip()
    })

    found = detail["FOUND"].eq("YES")
    ranks = pd.to_numeric(detail["GLOBAL_RANK"], errors="coerce")

    panel_ranks = pd.to_numeric(
        detail["WITHIN_PANEL_RANK"],
        errors="coerce",
    )

    summary_rows = [
        {
            "METRIC": "truth_candidates",
            "VALUE": len(detail),
        },
        {
            "METRIC": "truth_candidates_found",
            "VALUE": int(found.sum()),
        },
        {
            "METRIC": "median_global_rank_found",
            "VALUE": (
                float(ranks[found].median())
                if found.any()
                else "."
            ),
        },
        {
            "METRIC": "median_within_panel_rank_found",
            "VALUE": (
                float(panel_ranks[found].median())
                if panel_ranks[found].notna().any()
                else "."
            ),
        },
        {
            "METRIC": "median_rank_percentile_found",
            "VALUE": (
                float(
                    pd.to_numeric(
                        detail.loc[found, "GLOBAL_RANK_PERCENTILE"],
                        errors="coerce",
                    ).median()
                )
                if found.any()
                else "."
            ),
        },
    ]

    denominator = max(len(detail), 1)
    for k in cutoffs:
        recalled_global = int((ranks <= k).fillna(False).sum())
        recalled_panel = int((panel_ranks <= k).fillna(False).sum())
        summary_rows.extend(
            [
                {
                    "METRIC": f"global_top_{k}_recall",
                    "VALUE": round(recalled_global / denominator, 6),
                },
                {
                    "METRIC": f"within_panel_top_{k}_recall",
                    "VALUE": round(recalled_panel / denominator, 6),
                },
                {
                    "METRIC": f"global_top_{k}_recalled_truth_candidates",
                    "VALUE": recalled_global,
                },
                {
                    "METRIC": f"within_panel_top_{k}_recalled_truth_candidates",
                    "VALUE": recalled_panel,
                },
            ]
        )

    detail_path = Path(args.detail_output)
    summary_path = Path(args.summary_output)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    print(
        f"[OK] truth_candidates={len(detail)} found={int(found.sum())} "
        f"detail={detail_path} summary={summary_path}"
    )


if __name__ == "__main__":
    main()
