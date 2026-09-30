#!/usr/bin/env python3
"""Summarize population-frequency context and SV-gene effects.

needLR AF is a direct long-read population annotation. AnnotSV's benign AFmax
can include gnomAD-SV among several benign resources and is therefore retained
as BENIGN_DB_AFMAX rather than mislabeled as a source-specific gnomAD AF.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def gnomad_context(row):
    overlap = str(row.get("GNOMAD_SV_OVERLAP", ".")).upper()
    af = number(row.get("BENIGN_DB_AFMAX"))
    if overlap == "YES":
        if af is None:
            return "GNOMAD_OVERLAP_AF_UNKNOWN"
        return "GNOMAD_CONTEXT_RARE_AFMAX_LE_0.01" if af <= 0.01 else "GNOMAD_CONTEXT_COMMON_AFMAX_GT_0.01"
    if overlap == "NOT_REPORTED":
        return "GNOMAD_NOT_REPORTED_IN_ANNOTSV_OVERLAPS"
    if overlap in {"NOT_APPLICABLE", "UNKNOWN", ".", ""}:
        return "GNOMAD_CONTEXT_UNKNOWN"
    return "GNOMAD_CONTEXT_UNKNOWN"


def combined_class(row):
    needlr = str(row.get("POPULATION_STATUS", "UNKNOWN")).upper()
    gnomad = gnomad_context(row)

    if needlr == "COMMON":
        return "COMMON_NEEDLR"
    if "COMMON_AFMAX" in gnomad:
        if needlr == "RARE":
            return "DISCORDANT_RARE_NEEDLR_COMMON_GNOMAD_CONTEXT"
        return "COMMON_GNOMAD_CONTEXT"
    if needlr == "RARE":
        return "RARE_NEEDLR"
    if needlr == "NO_POPULATION_MATCH":
        return "NO_NEEDLR_MATCH"
    if needlr.startswith("NOT_EVALUABLE"):
        return needlr
    return "POPULATION_FREQUENCY_UNKNOWN"


def effect_group(effect, svtype):
    effect = str(effect or "").upper()
    svtype = str(svtype or "").upper()
    if "BREAKPOINT_IN_" in effect or "TWO_BREAKPOINTS_IN_GENE" in effect:
        return "DIRECT_BREAKPOINT"
    if "WHOLE_GENE_DELETION" in effect or "PARTIAL_GENE_DELETION" in effect:
        return "COPY_LOSS"
    if "WHOLE_GENE_DUPLICATION" in effect or "PARTIAL_GENE_DUPLICATION" in effect:
        return "COPY_GAIN"
    if effect.startswith("INSERTION_IN_"):
        return "INSERTION_IN_GENE"
    if "SPANNED_BY_INVERSION" in effect or "INSIDE_INVERSION" in effect:
        return "INVERSION_SPANNED_GENE"
    if "NEAR_GENE" in effect:
        return "NEAR_GENE"
    if svtype == "INV":
        return "INVERSION_OTHER"
    return "OTHER_OR_UNRESOLVED"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--gene-summary", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.candidates, sep="\t", dtype=str, low_memory=False)
    if df.empty:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output, sep="\t", index=False)
        pd.DataFrame().to_csv(args.gene_summary, sep="\t", index=False)
        return

    out = df.copy()
    out["GNOMAD_CONTEXT_CLASS"] = out.apply(gnomad_context, axis=1)
    out["POPULATION_CLASS"] = out.apply(combined_class, axis=1)
    out["SV_EFFECT_GROUP"] = [
        effect_group(effect, svtype)
        for effect, svtype in zip(out["SV_GENE_EFFECT"], out["SVTYPE"])
    ]
    out["POPULATION_INTERPRETATION"] = (
        "needLR AF and AnnotSV gnomAD/benign-database context are population evidence only; "
        "common frequency can argue against a highly penetrant rare-disease allele, while rarity "
        "alone does not establish pathogenicity. BENIGN_DB_AFMAX is not guaranteed to be a "
        "source-specific gnomAD AF when multiple benign resources overlap."
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    unique = out.drop_duplicates(["GENE", "SV_ID"])
    summary_rows = []
    for gene, group in unique.groupby("GENE", sort=False):
        classes = group["POPULATION_CLASS"].value_counts().to_dict()
        effects = group["SV_EFFECT_GROUP"].value_counts().to_dict()
        score = pd.to_numeric(group["GENE_RELEVANCE_SCORE"], errors="coerce").max()
        summary_rows.append({
            "GENE": gene,
            "GENE_RELEVANCE_SCORE": score if pd.notna(score) else ".",
            "PANEL_STATUS": "PANEL_GENE" if (group["PANEL_STATUS"] == "PANEL_GENE").any() else "NON_PANEL",
            "SV_COUNT": int(group["SV_ID"].nunique()),
            "RARE_NEEDLR_COUNT": int((group["POPULATION_CLASS"] == "RARE_NEEDLR").sum()),
            "COMMON_POPULATION_COUNT": int(group["POPULATION_CLASS"].str.startswith("COMMON").sum()),
            "GNOMAD_OVERLAP_COUNT": int((group["GNOMAD_SV_OVERLAP"] == "YES").sum()),
            "POPULATION_CLASSES": ";".join(f"{k}={v}" for k, v in sorted(classes.items())),
            "SV_EFFECTS": ";".join(f"{k}={v}" for k, v in sorted(effects.items())),
        })

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["_score"] = pd.to_numeric(summary["GENE_RELEVANCE_SCORE"], errors="coerce").fillna(0)
        summary = summary.sort_values(
            ["RARE_NEEDLR_COUNT", "_score", "SV_COUNT", "GENE"],
            ascending=[False, False, False, True],
        ).drop(columns="_score")
    summary.to_csv(args.gene_summary, sep="\t", index=False)

    print(
        f"[OK] population_rows={len(out)} genes={out['GENE'].nunique()} "
        f"output={output} gene_summary={args.gene_summary}"
    )


if __name__ == "__main__":
    main()
