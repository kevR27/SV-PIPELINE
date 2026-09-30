#!/usr/bin/env python3
"""Summarize population-frequency evidence for final SV-gene candidates.

needLR AF is treated as a coordinate-matched ONT control-frequency estimate.
AnnotSV B_*_AFmax is the maximum AF among reported benign overlapping regions
and may combine gnomAD with other sources. It is therefore not renamed or
interpreted as an exact gnomAD allele frequency.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def as_num(series):
    return pd.to_numeric(series, errors="coerce")


def freq_class(value):
    if pd.isna(value):
        return "NO_FREQUENCY"
    value = float(value)
    if value == 0:
        return "NOT_OBSERVED"
    if value <= 0.001:
        return "VERY_RARE_LE_0.001"
    if value <= 0.01:
        return "RARE_LE_0.01"
    if value < 0.05:
        return "COMMON_GT_0.01"
    return "VERY_COMMON_GE_0.05"


def population_interpretation(row):
    needlr = row["_needlr_af"]
    benign = row["_benign_afmax"]

    if pd.notna(needlr) and needlr >= 0.05:
        return "HIGH_FREQUENCY_WEAKENS_RARE_MENDELIAN_CANDIDACY"
    if pd.notna(needlr) and needlr > 0.01:
        return "COMMON_FREQUENCY_WEAKENS_RARE_MENDELIAN_CANDIDACY"

    # AnnotSV B_* AFmax is overlap evidence from benign population resources,
    # not guaranteed exact-allele equivalence.
    if pd.notna(benign) and benign >= 0.05:
        return "HIGH_FREQUENCY_BENIGN_REGION_OVERLAP_REVIEW"
    if pd.notna(benign) and benign > 0.01:
        return "COMMON_BENIGN_REGION_OVERLAP_REVIEW"

    if pd.notna(needlr) and needlr <= 0.01:
        return "LOW_FREQUENCY_DOES_NOT_ESTABLISH_PATHOGENICITY"
    if pd.notna(benign) and benign <= 0.01:
        return "LOW_FREQUENCY_BENIGN_REGION_OVERLAP"

    status = str(row.get("POPULATION_STATUS", "")).upper()
    if "NOT_EVALUABLE" in status:
        return "POPULATION_AF_NOT_EVALUABLE"
    if "NO_POPULATION_MATCH" in status:
        return "NO_NEEDLR_MATCH_NOT_EQUIVALENT_TO_AF_ZERO"
    return "INSUFFICIENT_POPULATION_FREQUENCY_DATA"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--sv-output", required=True)
    p.add_argument("--gene-output", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.candidates, sep="\t", dtype=str, low_memory=False)

    if df.empty:
        Path(args.sv_output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.sv_output, sep="\t", index=False)
        pd.DataFrame().to_csv(args.gene_output, sep="\t", index=False)
        return

    out = df.copy()
    out["_needlr_af"] = as_num(out.get("NEEDLR_AF", pd.Series(index=out.index, dtype=float)))
    out["_benign_afmax"] = as_num(
        out.get(
            "ANNOTSV_BENIGN_AFMAX",
            out.get("BENIGN_DB_AFMAX", pd.Series(index=out.index, dtype=float)),
        )
    )

    source = out.get(
        "ANNOTSV_BENIGN_DB_SOURCE",
        out.get("BENIGN_DB_SOURCE", pd.Series(".", index=out.index)),
    ).fillna(".").astype(str)

    out["NEEDLR_AF_CLASS"] = out["_needlr_af"].map(freq_class)
    out["ANNOTSV_BENIGN_AFMAX_CLASS"] = out["_benign_afmax"].map(freq_class)
    out["GNOMAD_INCLUDED_IN_BENIGN_SOURCE"] = np.where(
        source.str.contains("gnomAD", case=False, regex=False),
        "YES",
        np.where(source.isin([".", "", "nan", "None"]), "UNKNOWN", "NO"),
    )
    out["GNOMAD_AF_SCOPE"] = np.where(
        out["GNOMAD_INCLUDED_IN_BENIGN_SOURCE"].eq("YES"),
        "GNOMAD_OVERLAP_PRESENT_BUT_AFMAX_MAY_INCLUDE_OTHER_BENIGN_SOURCES",
        "NO_GNOMAD_SPECIFIC_AF_AVAILABLE_FROM_CURRENT_ANNOTSV_FIELDS",
    )
    out["POPULATION_EVIDENCE_INTERPRETATION"] = out.apply(
        population_interpretation,
        axis=1,
    )

    out = out.drop(columns=["_needlr_af", "_benign_afmax"])

    group_cols = ["GENE"]
    gene_rows = []
    for gene, group in out.groupby("GENE", sort=False):
        needlr = as_num(group["NEEDLR_AF"])
        benign = as_num(
            group.get(
                "ANNOTSV_BENIGN_AFMAX",
                group.get("BENIGN_DB_AFMAX", pd.Series(index=group.index, dtype=float)),
            )
        )

        interpretations = group["POPULATION_EVIDENCE_INTERPRETATION"].fillna(".")
        gene_rows.append({
            "GENE": gene,
            "SV_COUNT": int(group["SV_ID"].nunique()),
            "RARE_OR_NOT_OBSERVED_NEEDLR_SV_COUNT": int(
                group.loc[
                    group["NEEDLR_AF_CLASS"].isin(
                        ["NOT_OBSERVED", "VERY_RARE_LE_0.001", "RARE_LE_0.01"]
                    ),
                    "SV_ID",
                ].nunique()
            ),
            "COMMON_NEEDLR_SV_COUNT": int(
                group.loc[
                    group["NEEDLR_AF_CLASS"].isin(
                        ["COMMON_GT_0.01", "VERY_COMMON_GE_0.05"]
                    ),
                    "SV_ID",
                ].nunique()
            ),
            "GNOMAD_BENIGN_OVERLAP_SV_COUNT": int(
                group.loc[
                    group["GNOMAD_INCLUDED_IN_BENIGN_SOURCE"].eq("YES"),
                    "SV_ID",
                ].nunique()
            ),
            "MAX_NEEDLR_AF": (
                round(float(needlr.max()), 6)
                if needlr.notna().any()
                else "."
            ),
            "MAX_ANNOTSV_BENIGN_AFMAX": (
                round(float(benign.max()), 6)
                if benign.notna().any()
                else "."
            ),
            "PANEL_STATUS": (
                "PANEL_GENE"
                if group.get("PANEL_STATUS", pd.Series("", index=group.index))
                .fillna("")
                .eq("PANEL_GENE")
                .any()
                else "NON_PANEL"
            ),
            "GENE_RELEVANCE_SCORE": (
                pd.to_numeric(
                    group.get(
                        "GENE_RELEVANCE_SCORE",
                        pd.Series(index=group.index, dtype=float),
                    ),
                    errors="coerce",
                ).max()
            ),
            "POPULATION_INTERPRETATIONS": ";".join(
                sorted(set(x for x in interpretations.astype(str) if x not in {"", "."}))
            ) or ".",
        })

    gene = pd.DataFrame(gene_rows)
    if not gene.empty:
        gene["GENE_RELEVANCE_SCORE"] = pd.to_numeric(
            gene["GENE_RELEVANCE_SCORE"],
            errors="coerce",
        ).fillna(0)
        gene = gene.sort_values(
            [
                "GENE_RELEVANCE_SCORE",
                "RARE_OR_NOT_OBSERVED_NEEDLR_SV_COUNT",
                "COMMON_NEEDLR_SV_COUNT",
                "GENE",
            ],
            ascending=[False, False, True, True],
        )

    sv_path = Path(args.sv_output)
    gene_path = Path(args.gene_output)
    sv_path.parent.mkdir(parents=True, exist_ok=True)
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(sv_path, sep="\t", index=False)
    gene.to_csv(gene_path, sep="\t", index=False)

    print(
        f"[OK] sv_gene_rows={len(out)} genes={len(gene)} "
        f"sv_output={sv_path} gene_output={gene_path}"
    )


if __name__ == "__main__":
    main()
