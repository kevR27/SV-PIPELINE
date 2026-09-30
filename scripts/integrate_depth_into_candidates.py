#!/usr/bin/env python3
"""Attach large DEL/DUP read-depth evidence to final SV-gene candidate tables.

Depth evidence is event-level support from the same sequencing data. It does
not validate pathogenicity and it does not change the master SV callset.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def num(series):
    return pd.to_numeric(series, errors="coerce")


def support_class(pattern):
    value = str(pattern or "").strip().upper()
    if value in {"CONSISTENT_WITH_LOSS", "CONSISTENT_WITH_GAIN"}:
        return "SUPPORTS_CALLED_COPY_CHANGE"
    if value == "NOT_CLEAR":
        return "DOES_NOT_SHOW_EXPECTED_COPY_CHANGE"
    if value in MISSING:
        return "NO_DEPTH_RESULT"
    return "REVIEW"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gene-candidates", required=True)
    p.add_argument("--sv-candidates", required=True)
    p.add_argument("--depth-summary", required=True)
    p.add_argument("--hpo-similarity", default=None)
    p.add_argument("--gene-output", required=True)
    p.add_argument("--sv-output", required=True)
    p.add_argument("--min-size", type=int, default=100000)
    args = p.parse_args()

    genes = pd.read_csv(args.gene_candidates, sep="\t", dtype=str, low_memory=False)
    sv = pd.read_csv(args.sv_candidates, sep="\t", dtype=str, low_memory=False)
    depth = pd.read_csv(args.depth_summary, sep="\t", dtype=str, low_memory=False)
    hpo = (
        pd.read_csv(args.hpo_similarity, sep="\t", dtype=str, low_memory=False)
        if args.hpo_similarity
        else pd.DataFrame()
    )

    depth_cols = [
        "MEDIAN_DEPTH_INSIDE",
        "MEDIAN_DEPTH_FLANKS",
        "DEPTH_RATIO",
        "DEPTH_PATTERN",
        "PLOT_BIN_BP",
    ]
    if depth.empty:
        depth_index = pd.DataFrame(columns=["SV_ID", *depth_cols])
    else:
        available = ["SV_ID", *[c for c in depth_cols if c in depth.columns]]
        depth_index = depth[available].drop_duplicates("SV_ID")

    out = sv.merge(depth_index, on="SV_ID", how="left")
    span = num(out["SV_SPAN_BP"]) if "SV_SPAN_BP" in out else num(out["SVLEN"]).abs()
    is_large_copy = (
        out["SVTYPE"].fillna("").isin(["DEL", "DUP"])
        & span.ge(args.min_size)
    )

    matched_depth = out["DEPTH_PATTERN"].notna() if "DEPTH_PATTERN" in out else pd.Series(False, index=out.index)
    out["DEPTH_EVALUATION"] = "NOT_APPLICABLE"
    out.loc[is_large_copy & ~matched_depth, "DEPTH_EVALUATION"] = "NO_DEPTH_SUMMARY"
    out.loc[is_large_copy & matched_depth, "DEPTH_EVALUATION"] = "EVALUATED"

    if "DEPTH_PATTERN" not in out:
        out["DEPTH_PATTERN"] = "."
    out["DEPTH_PATTERN"] = out["DEPTH_PATTERN"].fillna(".")
    out["DEPTH_SUPPORT_CLASS"] = out["DEPTH_PATTERN"].map(support_class)
    out.loc[~is_large_copy, "DEPTH_SUPPORT_CLASS"] = "NOT_APPLICABLE"

    for col in depth_cols:
        if col not in out:
            out[col] = "."
        else:
            out[col] = out[col].fillna(".")

    if not hpo.empty and "GENE" in hpo.columns:
        hpo_cols = [
            "GENE",
            "HPO_SEMANTIC_STATUS",
            "PATIENT_HPO_COUNT",
            "GENE_REFERENCE_HPO_COUNT",
            "HPO_EXACT_MATCH_COUNT",
            "HPO_BMA_RESNIK",
            "HPO_BMA_RESNIK_NORMALIZED",
            "BEST_MATCHED_PATIENT_HPO",
        ]
        hpo_cols = [col for col in hpo_cols if col in hpo.columns]
        hpo_index = hpo[hpo_cols].drop_duplicates("GENE")
        out = out.merge(hpo_index, on="GENE", how="left")
        genes = genes.merge(hpo_index, on="GENE", how="left")

    if not genes.empty and not out.empty:
        unique = out.drop_duplicates(["GENE", "SV_ID"])
        stats = (
            unique.groupby("GENE")
            .agg(
                DEPTH_EVALUATED_SV_COUNT=(
                    "DEPTH_EVALUATION",
                    lambda s: int((s == "EVALUATED").sum()),
                ),
                DEPTH_SUPPORTED_COPY_CHANGE_COUNT=(
                    "DEPTH_SUPPORT_CLASS",
                    lambda s: int((s == "SUPPORTS_CALLED_COPY_CHANGE").sum()),
                ),
                DEPTH_NOT_CLEAR_SV_COUNT=(
                    "DEPTH_SUPPORT_CLASS",
                    lambda s: int((s == "DOES_NOT_SHOW_EXPECTED_COPY_CHANGE").sum()),
                ),
            )
            .reset_index()
        )
        genes = genes.merge(stats, on="GENE", how="left")
    for col in [
        "DEPTH_EVALUATED_SV_COUNT",
        "DEPTH_SUPPORTED_COPY_CHANGE_COUNT",
        "DEPTH_NOT_CLEAR_SV_COUNT",
    ]:
        if col not in genes:
            genes[col] = 0
        genes[col] = pd.to_numeric(genes[col], errors="coerce").fillna(0).astype(int)

    gene_path = Path(args.gene_output)
    sv_path = Path(args.sv_output)
    gene_path.parent.mkdir(parents=True, exist_ok=True)
    sv_path.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(gene_path, sep="\t", index=False)
    out.to_csv(sv_path, sep="\t", index=False)

    print(
        f"[OK] sv_rows={len(out)} depth_evaluated="
        f"{int((out['DEPTH_EVALUATION'] == 'EVALUATED').sum())} "
        f"sv_output={sv_path}"
    )


if __name__ == "__main__":
    main()
