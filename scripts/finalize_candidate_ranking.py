#!/usr/bin/env python3
"""Finalize patient-aware SV/gene ranking after phased SNV+SV review.

This is the last small candidate-ranking layer. It does not alter the master SV
callset. It can promote a recessive SV candidate only when an eligible same-gene
small variant is phased in trans; unresolved/cis configurations remain review
states.

Panel and non-panel candidates receive separate ranks. Panel membership is never
used as a score.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from ranking_common import (
    GENE_TIER_PRIORITY,
    event_sort_tuple,
)


PAIR_PRIORITY = {
    "AR_TRANS_SNV_SV_CANDIDATE": 3,
    "AR_SECOND_ALLELE_CANDIDATE_PHASE_UNRESOLVED": 2,
    "AR_CIS_NOT_BIALLELIC_BY_PHASE": 1,
    "NOT_AR_PAIRING_MODEL": 0,
}


def best_pair_status(group: pd.DataFrame) -> str:
    if group.empty or "RECESSIVE_PAIR_STATUS" not in group.columns:
        return "."
    values = (
        group["RECESSIVE_PAIR_STATUS"]
        .fillna(".")
        .astype(str)
        .tolist()
    )
    return max(values, key=lambda value: PAIR_PRIORITY.get(value, -1))


def final_sort_row(row: pd.Series):
    values = row.to_dict()

    if str(values.get("FINAL_GENE_RELEVANCE_TIER", ".")).upper() not in {
        "", ".", "NA", "N/A", "NAN", "NONE", "NULL",
    }:
        values["GENE_RELEVANCE_TIER"] = values[
            "FINAL_GENE_RELEVANCE_TIER"
        ]

    if str(values.get("FINAL_GENE_RELEVANCE_DISPLAY_SCORE", ".")).upper() not in {
        "", ".", "NA", "N/A", "NAN", "NONE", "NULL",
    }:
        values["GENE_RELEVANCE_DISPLAY_SCORE"] = values[
            "FINAL_GENE_RELEVANCE_DISPLAY_SCORE"
        ]

    if str(values.get("FINAL_INHERITANCE_MECHANISM_CLASS", ".")).upper() not in {
        "", ".", "NA", "N/A", "NAN", "NONE", "NULL",
    }:
        values["INHERITANCE_MECHANISM_CLASS"] = values[
            "FINAL_INHERITANCE_MECHANISM_CLASS"
        ]

    return event_sort_tuple(pd.Series(values))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genes", required=True)
    p.add_argument("--sv", required=True)
    p.add_argument("--snv-sv", required=True)
    p.add_argument("--gene-output", required=True)
    p.add_argument("--sv-output", required=True)
    args = p.parse_args()

    genes = pd.read_csv(args.genes, sep="	", dtype=str, low_memory=False)
    sv = pd.read_csv(args.sv, sep="	", dtype=str, low_memory=False)
    pairs = pd.read_csv(args.snv_sv, sep="	", dtype=str, low_memory=False)

    if sv.empty:
        Path(args.sv_output).parent.mkdir(parents=True, exist_ok=True)
        sv.to_csv(args.sv_output, sep="	", index=False)
        genes.to_csv(args.gene_output, sep="	", index=False)
        print(f"[OK] ranked_sv_rows=0 ranked_genes={len(genes)}")
        return

    pair_by_event = {}
    if not pairs.empty and {"GENE", "SV_ID"}.issubset(pairs.columns):
        for key, group in pairs.groupby(["GENE", "SV_ID"], sort=False):
            pair_by_event[key] = best_pair_status(group)

    sv["RECESSIVE_PAIR_STATUS"] = [
        pair_by_event.get((str(gene), str(sv_id)), ".")
        for gene, sv_id in zip(sv["GENE"], sv["SV_ID"])
    ]

    final_class = []
    final_detail = []
    for _, row in sv.iterrows():
        base = str(
            row.get("INHERITANCE_MECHANISM_CLASS", "UNRESOLVED")
        )
        pair = str(row.get("RECESSIVE_PAIR_STATUS", "."))
        if pair == "AR_TRANS_SNV_SV_CANDIDATE":
            final_class.append("AR_TRANS_SECOND_ALLELE_SUPPORTED")
            final_detail.append(
                "A same-gene small variant is phased in trans with the SV; "
                "this supports a biallelic candidate model but does not prove "
                "either allele is pathogenic or that segregation is correct."
            )
        elif pair == "AR_CIS_NOT_BIALLELIC_BY_PHASE":
            final_class.append(base)
            final_detail.append(
                "A same-gene small variant is phased in cis with the SV; this "
                "does not provide the second allele required for an AR model."
            )
        elif pair == "AR_SECOND_ALLELE_CANDIDATE_PHASE_UNRESOLVED":
            final_class.append(base)
            final_detail.append(
                "A same-gene small variant exists but phase is unresolved; "
                "the AR second-allele requirement remains unresolved."
            )
        else:
            final_class.append(base)
            final_detail.append(
                str(row.get("INHERITANCE_MECHANISM_DETAIL", "."))
            )

    sv["FINAL_INHERITANCE_MECHANISM_CLASS"] = final_class
    sv["FINAL_INHERITANCE_MECHANISM_DETAIL"] = final_detail

    order = sorted(
        range(len(sv)),
        key=lambda idx: final_sort_row(sv.iloc[idx]),
    )
    sv = sv.iloc[order].reset_index(drop=True)
    sv["FINAL_EVENT_RANK_GLOBAL"] = range(1, len(sv) + 1)

    if "PANEL_STATUS" not in sv.columns:
        sv["PANEL_STATUS"] = "NONPANEL_GENE"
    sv["PANEL_STATUS"] = (
        sv["PANEL_STATUS"]
        .fillna("NONPANEL_GENE")
        .astype(str)
        .replace({"NON_PANEL": "NONPANEL_GENE"})
    )
    sv["FINAL_EVENT_RANK_WITHIN_PANEL_STATUS"] = (
        sv.groupby("PANEL_STATUS").cumcount() + 1
    )

    # Gene ranking is driven by its final phenotype tier plus the best event.
    if not genes.empty:
        best_event = (
            sv.sort_values("FINAL_EVENT_RANK_GLOBAL")
            .drop_duplicates("GENE")
            [[
                "GENE",
                "FINAL_EVENT_RANK_GLOBAL",
                "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS",
                "FINAL_INHERITANCE_MECHANISM_CLASS",
                "RECESSIVE_PAIR_STATUS",
            ]]
            .rename(columns={
                "FINAL_EVENT_RANK_GLOBAL": "BEST_EVENT_FINAL_RANK_GLOBAL",
                "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS": (
                    "BEST_EVENT_FINAL_RANK_WITHIN_PANEL_STATUS"
                ),
                "FINAL_INHERITANCE_MECHANISM_CLASS": (
                    "BEST_EVENT_INHERITANCE_MECHANISM_CLASS"
                ),
                "RECESSIVE_PAIR_STATUS": "BEST_RECESSIVE_PAIR_STATUS",
            })
        )
        genes = genes.merge(best_event, on="GENE", how="left")

        if "PANEL_STATUS" not in genes.columns:
            genes["PANEL_STATUS"] = "NONPANEL_GENE"
        genes["PANEL_STATUS"] = (
            genes["PANEL_STATUS"]
            .fillna("NONPANEL_GENE")
            .astype(str)
            .replace({"NON_PANEL": "NONPANEL_GENE"})
        )

        tier_col = (
            "FINAL_GENE_RELEVANCE_TIER"
            if "FINAL_GENE_RELEVANCE_TIER" in genes.columns
            else "GENE_RELEVANCE_TIER"
        )
        score_col = (
            "FINAL_GENE_RELEVANCE_DISPLAY_SCORE"
            if "FINAL_GENE_RELEVANCE_DISPLAY_SCORE" in genes.columns
            else "GENE_RELEVANCE_DISPLAY_SCORE"
        )

        genes["_tier_rank"] = (
            genes[tier_col]
            .fillna("LIMITED")
            .astype(str)
            .str.upper()
            .map(GENE_TIER_PRIORITY)
            .fillna(0)
        )
        genes["_score"] = pd.to_numeric(
            genes.get(score_col, 0),
            errors="coerce",
        ).fillna(0)
        genes["_best_event_rank"] = pd.to_numeric(
            genes.get("BEST_EVENT_FINAL_RANK_GLOBAL", "."),
            errors="coerce",
        ).fillna(10**9)

        genes = genes.sort_values(
            [
                "PANEL_STATUS",
                "_tier_rank",
                "_score",
                "_best_event_rank",
                "GENE",
            ],
            ascending=[True, False, False, True, True],
        ).reset_index(drop=True)
        genes["FINAL_GENE_RANK_WITHIN_PANEL_STATUS"] = (
            genes.groupby("PANEL_STATUS").cumcount() + 1
        )
        genes = genes.drop(
            columns=["_tier_rank", "_score", "_best_event_rank"]
        )

    gene_output = Path(args.gene_output)
    sv_output = Path(args.sv_output)
    gene_output.parent.mkdir(parents=True, exist_ok=True)
    sv_output.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(gene_output, sep="	", index=False)
    sv.to_csv(sv_output, sep="	", index=False)

    print(
        f"[OK] ranked_sv_rows={len(sv)} ranked_genes={len(genes)} "
        f"ar_trans_events={(sv['RECESSIVE_PAIR_STATUS'] == 'AR_TRANS_SNV_SV_CANDIDATE').sum()} "
        f"sv_output={sv_output}"
    )


if __name__ == "__main__":
    main()
