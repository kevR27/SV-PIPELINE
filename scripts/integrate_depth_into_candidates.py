#!/usr/bin/env python3
"""Attach large DEL/DUP read-depth evidence to final SV-gene candidate tables.

Depth evidence is event-level support from the same sequencing data. It does
not validate pathogenicity and it does not change the master SV callset.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from ranking_common import event_sort_tuple, gene_relevance


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
        "DEPTH_INTERPRETATION_SCOPE",
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

    # Patient-specific phenotype ranking is applied only when semantic
    # similarity was actually evaluated. Otherwise retain the generic HON
    # relevance tier from the upstream gene ranking.
    if not out.empty:
        patient_similarity = pd.to_numeric(
            out.get(
                "HPO_BMA_RESNIK_NORMALIZED",
                pd.Series(index=out.index, dtype=float),
            ),
            errors="coerce",
        )
        semantic_status = out.get(
            "HPO_SEMANTIC_STATUS",
            pd.Series("PATIENT_HPO_NOT_AVAILABLE", index=out.index),
        ).fillna("PATIENT_HPO_NOT_AVAILABLE").astype(str)

        disease_score = pd.to_numeric(
            out.get(
                "GENE_DISEASE_SCORE",
                out.get(
                    "GENE_DISEASE_EVIDENCE_SCORE",
                    pd.Series(0, index=out.index),
                ),
            ),
            errors="coerce",
        ).fillna(0)

        final_tiers = []
        final_scores = []
        final_scopes = []
        for idx in out.index:
            if (
                semantic_status.loc[idx] == "EVALUATED"
                and pd.notna(patient_similarity.loc[idx])
            ):
                rel = gene_relevance(
                    float(patient_similarity.loc[idx]),
                    float(disease_score.loc[idx]),
                )
                final_tiers.append(rel["tier"])
                final_scores.append(rel["display_score"])
                final_scopes.append("PATIENT_SPECIFIC_HPO_RESNIK_BMA")
            else:
                final_tiers.append(
                    str(out.loc[idx].get("GENE_RELEVANCE_TIER", "LIMITED"))
                )
                generic = pd.to_numeric(
                    pd.Series([
                        out.loc[idx].get(
                            "GENE_RELEVANCE_DISPLAY_SCORE",
                            out.loc[idx].get("GENE_RELEVANCE_SCORE", 0),
                        )
                    ]),
                    errors="coerce",
                ).iloc[0]
                final_scores.append(float(generic) if pd.notna(generic) else 0.0)
                final_scopes.append("GENERIC_HON_FALLBACK_NO_PATIENT_HPO")

        out["FINAL_GENE_RELEVANCE_TIER"] = final_tiers
        out["FINAL_GENE_RELEVANCE_DISPLAY_SCORE"] = final_scores
        out["FINAL_PHENOTYPE_RANKING_SCOPE"] = final_scopes

        # Reuse the shared event sort by substituting only the phenotype-aware
        # final gene tier/score into a temporary row view.
        order = sorted(
            range(len(out)),
            key=lambda i: event_sort_tuple(
                out.iloc[i].assign()
                if False
                else pd.Series({
                    **out.iloc[i].to_dict(),
                    "GENE_RELEVANCE_TIER": out.iloc[i]["FINAL_GENE_RELEVANCE_TIER"],
                    "GENE_RELEVANCE_DISPLAY_SCORE": out.iloc[i]["FINAL_GENE_RELEVANCE_DISPLAY_SCORE"],
                })
            ),
        )
        out = out.iloc[order].reset_index(drop=True)
        out["FINAL_EVENT_RANK_WITHIN_PANEL_STATUS"] = (
            out.groupby(
                out.get(
                    "PANEL_STATUS",
                    pd.Series("UNSPECIFIED", index=out.index),
                )
            ).cumcount()
            + 1
        )

    if not genes.empty:
        if not hpo.empty and "GENE" in hpo.columns:
            # genes already received the HPO columns above.
            pass

        gene_disease = pd.to_numeric(
            genes.get(
                "GENE_DISEASE_SCORE",
                genes.get(
                    "gene_disease_evidence_score",
                    pd.Series(0, index=genes.index),
                ),
            ),
            errors="coerce",
        ).fillna(0)
        gene_patient = pd.to_numeric(
            genes.get(
                "HPO_BMA_RESNIK_NORMALIZED",
                pd.Series(index=genes.index, dtype=float),
            ),
            errors="coerce",
        )
        gene_status = genes.get(
            "HPO_SEMANTIC_STATUS",
            pd.Series("PATIENT_HPO_NOT_AVAILABLE", index=genes.index),
        ).fillna("PATIENT_HPO_NOT_AVAILABLE").astype(str)

        final_tiers = []
        final_scores = []
        final_scopes = []
        for idx in genes.index:
            if gene_status.loc[idx] == "EVALUATED" and pd.notna(gene_patient.loc[idx]):
                rel = gene_relevance(
                    float(gene_patient.loc[idx]),
                    float(gene_disease.loc[idx]),
                )
                final_tiers.append(rel["tier"])
                final_scores.append(rel["display_score"])
                final_scopes.append("PATIENT_SPECIFIC_HPO_RESNIK_BMA")
            else:
                final_tiers.append(
                    str(genes.loc[idx].get("GENE_RELEVANCE_TIER", "LIMITED"))
                )
                generic = pd.to_numeric(
                    pd.Series([
                        genes.loc[idx].get(
                            "GENE_RELEVANCE_DISPLAY_SCORE",
                            genes.loc[idx].get("GENE_RELEVANCE_SCORE", 0),
                        )
                    ]),
                    errors="coerce",
                ).iloc[0]
                final_scores.append(float(generic) if pd.notna(generic) else 0.0)
                final_scopes.append("GENERIC_HON_FALLBACK_NO_PATIENT_HPO")

        genes["FINAL_GENE_RELEVANCE_TIER"] = final_tiers
        genes["FINAL_GENE_RELEVANCE_DISPLAY_SCORE"] = final_scores
        genes["FINAL_PHENOTYPE_RANKING_SCOPE"] = final_scopes

        tier_rank = {
            "HIGH": 3,
            "MODERATE": 2,
            "SUPPORTING": 1,
            "LIMITED": 0,
        }
        genes["_final_tier_rank"] = (
            genes["FINAL_GENE_RELEVANCE_TIER"]
            .fillna("LIMITED")
            .astype(str)
            .str.upper()
            .map(tier_rank)
            .fillna(0)
        )
        genes["_final_score"] = pd.to_numeric(
            genes["FINAL_GENE_RELEVANCE_DISPLAY_SCORE"],
            errors="coerce",
        ).fillna(0)
        panel_group = genes.get(
            "PANEL_STATUS",
            pd.Series("NONPANEL_GENE", index=genes.index),
        )
        genes["_panel_group"] = panel_group.fillna("NONPANEL_GENE").astype(str)
        genes = genes.sort_values(
            ["_panel_group", "_final_tier_rank", "_final_score", "GENE"],
            ascending=[True, False, False, True],
        ).reset_index(drop=True)
        genes["FINAL_GENE_RANK_WITHIN_PANEL_STATUS"] = (
            genes.groupby("_panel_group").cumcount() + 1
        )
        genes = genes.drop(
            columns=["_final_tier_rank", "_final_score", "_panel_group"]
        )

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
