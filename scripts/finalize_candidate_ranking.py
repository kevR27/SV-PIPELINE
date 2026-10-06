#!/usr/bin/env python3
"""Finalize patient-aware SV/gene ranking after phased SNV+SV review.

This is the last small candidate-ranking layer. It does not alter the master SV
callset. It can promote a recessive SV candidate only when an eligible same-gene
small variant is phased in trans; unresolved/cis configurations remain review
states.

Panel and non-panel candidates receive separate ranks. Panel membership is never
used as a score. Cohort recurrence is retained only as a late cautionary
tie-break and never removes a candidate.
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
    p.add_argument("--sample", required=True)
    p.add_argument("--recurrence-summary", default=None)
    p.add_argument("--recurrence-members", default=None)
    p.add_argument("--mitochondrial-ranking", default=None)
    p.add_argument("--gene-output", required=True)
    p.add_argument("--sv-output", required=True)
    p.add_argument("--mito-output", default=None)
    args = p.parse_args()

    genes = pd.read_csv(args.genes, sep="	", dtype=str, low_memory=False)
    sv = pd.read_csv(args.sv, sep="	", dtype=str, low_memory=False)
    pairs = pd.read_csv(args.snv_sv, sep="	", dtype=str, low_memory=False)

    mito = (
        pd.read_csv(
            args.mitochondrial_ranking,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        if args.mitochondrial_ranking
        else pd.DataFrame()
    )

    if args.recurrence_summary and args.recurrence_members:
        recurrence_summary = pd.read_csv(
            args.recurrence_summary,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        recurrence_members = pd.read_csv(
            args.recurrence_members,
            sep="\t",
            dtype=str,
            low_memory=False,
        )
        sample_members = recurrence_members[
            recurrence_members["SAMPLE"].astype(str).eq(str(args.sample))
        ].copy()

        if not sample_members.empty:
            sample_members = sample_members[
                ["COHORT_SV_ID", "SV_ID", "COHORT_RECURRENCE"]
            ].drop_duplicates("SV_ID")

            if not recurrence_summary.empty:
                counts = recurrence_summary[
                    ["COHORT_SV_ID", "SAMPLE_COUNT"]
                ].drop_duplicates("COHORT_SV_ID")
                sample_members = sample_members.merge(
                    counts,
                    on="COHORT_SV_ID",
                    how="left",
                )
                sample_members = sample_members.rename(
                    columns={"SAMPLE_COUNT": "COHORT_SAMPLE_COUNT"}
                )

            sv = sv.merge(sample_members, on="SV_ID", how="left")
            sv["COHORT_RECURRENCE"] = sv[
                "COHORT_RECURRENCE"
            ].fillna("NOT_MAPPED")
            if "COHORT_SAMPLE_COUNT" not in sv.columns:
                sv["COHORT_SAMPLE_COUNT"] = "."
            else:
                sv["COHORT_SAMPLE_COUNT"] = sv[
                    "COHORT_SAMPLE_COUNT"
                ].fillna(".")
        else:
            sv["COHORT_RECURRENCE"] = "NOT_MAPPED"
            sv["COHORT_SAMPLE_COUNT"] = "."
    else:
        sv["COHORT_RECURRENCE"] = "NOT_EVALUATED"
        sv["COHORT_SAMPLE_COUNT"] = "."

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

        # Final gene ordering keeps broad relevance tier first, then lets the
        # best inheritance/mechanism-aware SV event outrank small continuous
        # score differences. The continuous score is only a later tie-break.
        genes = genes.sort_values(
            [
                "PANEL_STATUS",
                "_tier_rank",
                "_best_event_rank",
                "_score",
                "GENE",
            ],
            ascending=[True, False, True, False, True],
        ).reset_index(drop=True)
        genes["FINAL_GENE_RANK_WITHIN_PANEL_STATUS"] = (
            genes.groupby("PANEL_STATUS").cumcount() + 1
        )
        genes = genes.drop(
            columns=["_tier_rank", "_score", "_best_event_rank"]
        )

    mito_ranked = mito.copy()
    if not mito_ranked.empty:
        final_gene_cols = [
            col
            for col in [
                "GENE",
                "FINAL_GENE_RELEVANCE_TIER",
                "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
                "FINAL_GENE_RANK_WITHIN_PANEL_STATUS",
                "BEST_EVENT_FINAL_RANK_GLOBAL",
                "BEST_EVENT_FINAL_RANK_WITHIN_PANEL_STATUS",
                "BEST_EVENT_INHERITANCE_MECHANISM_CLASS",
                "BEST_RECESSIVE_PAIR_STATUS",
            ]
            if col in genes.columns
        ]
        if "GENE" in final_gene_cols:
            mito_ranked = mito_ranked.merge(
                genes[final_gene_cols].drop_duplicates("GENE"),
                on="GENE",
                how="left",
            )

        mito_ranked["_final_tier_rank"] = (
            mito_ranked.get(
                "FINAL_GENE_RELEVANCE_TIER",
                mito_ranked.get(
                    "GENE_RELEVANCE_TIER",
                    pd.Series("LIMITED", index=mito_ranked.index),
                ),
            )
            .fillna(
                mito_ranked.get(
                    "GENE_RELEVANCE_TIER",
                    pd.Series("LIMITED", index=mito_ranked.index),
                )
            )
            .astype(str)
            .str.upper()
            .map(GENE_TIER_PRIORITY)
            .fillna(0)
        )
        mito_ranked["_final_score"] = pd.to_numeric(
            mito_ranked.get(
                "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
                mito_ranked.get(
                    "MAX_GENE_RELEVANCE",
                    pd.Series(0, index=mito_ranked.index),
                ),
            ),
            errors="coerce",
        ).fillna(
            pd.to_numeric(
                mito_ranked.get(
                    "MAX_GENE_RELEVANCE",
                    pd.Series(0, index=mito_ranked.index),
                ),
                errors="coerce",
            ).fillna(0)
        )
        mito_ranked["_best_event_rank"] = pd.to_numeric(
            mito_ranked.get(
                "BEST_EVENT_FINAL_RANK_GLOBAL",
                pd.Series(index=mito_ranked.index, dtype=float),
            ),
            errors="coerce",
        ).fillna(10**9)
        mito_ranked["_original_rank"] = pd.to_numeric(
            mito_ranked.get(
                "MITO_RANK_WITHIN_ENCODING",
                pd.Series(index=mito_ranked.index, dtype=float),
            ),
            errors="coerce",
        ).fillna(10**9)

        # Patient-aware gene/event evidence refines nuclear mitochondrial genes.
        # mtDNA rows usually have no nuclear gene-HPO row and therefore retain
        # their specialized original ordering.
        mito_ranked = mito_ranked.sort_values(
            [
                "ENCODING_GENOME",
                "_final_tier_rank",
                "_best_event_rank",
                "_final_score",
                "_original_rank",
                "GENE",
            ],
            ascending=[True, False, True, False, True, True],
        ).reset_index(drop=True)
        mito_ranked["MITO_FINAL_RANK_WITHIN_ENCODING"] = (
            mito_ranked.groupby("ENCODING_GENOME").cumcount() + 1
        )
        mito_ranked["MITO_FINAL_RANK_WITHIN_ENCODING_PANEL_STATUS"] = (
            mito_ranked.groupby(
                ["ENCODING_GENOME", "PANEL_STATUS"]
            ).cumcount()
            + 1
        )
        mito_ranked["MITO_FINAL_RANKING_SCOPE"] = [
            (
                "PATIENT_AWARE_NUCLEAR_MITOCHONDRIAL_RANK"
                if str(encoding).upper() == "NUCLEAR"
                and pd.notna(final_rank)
                else "SPECIALIZED_OR_GENERIC_MITOCHONDRIAL_RANK"
            )
            for encoding, final_rank in zip(
                mito_ranked["ENCODING_GENOME"],
                mito_ranked.get(
                    "FINAL_GENE_RANK_WITHIN_PANEL_STATUS",
                    pd.Series(index=mito_ranked.index, dtype=float),
                ),
            )
        ]
        mito_ranked = mito_ranked.drop(
            columns=[
                "_final_tier_rank",
                "_final_score",
                "_best_event_rank",
                "_original_rank",
            ],
            errors="ignore",
        )

    gene_output = Path(args.gene_output)
    sv_output = Path(args.sv_output)
    gene_output.parent.mkdir(parents=True, exist_ok=True)
    sv_output.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(gene_output, sep="	", index=False)
    sv.to_csv(sv_output, sep="	", index=False)

    if args.mito_output:
        mito_output = Path(args.mito_output)
        mito_output.parent.mkdir(parents=True, exist_ok=True)
        mito_ranked.to_csv(mito_output, sep="	", index=False)

    print(
        f"[OK] ranked_sv_rows={len(sv)} ranked_genes={len(genes)} "
        f"ar_trans_events={(sv['RECESSIVE_PAIR_STATUS'] == 'AR_TRANS_SNV_SV_CANDIDATE').sum()} "
        f"mito_genes={len(mito_ranked)} sv_output={sv_output}"
    )


if __name__ == "__main__":
    main()
