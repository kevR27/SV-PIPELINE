#!/usr/bin/env python3
"""Summarize population-frequency context and SV-gene effects.

This is a research-prioritization table, not a pathogenicity classifier.

Frequency sources are kept distinct:
- needLR AF: frequency from the matched long-read control resource used by
  needLR, transferred to the Jasmine event through the pipeline's conservative
  coordinate/type/size matching. This is provisional event-level population
  evidence, not exact allele identity.
- AnnotSV benign AFmax: maximum AF among overlapping benign SV regions reported
  by AnnotSV. The source field may include gnomAD-SV together with DGV, 1000G,
  ClinVar or other resources, so AFmax is not relabelled as a gnomAD-specific
  exact-allele frequency.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def needlr_class(row):
    af = number(row.get("NEEDLR_AF"))
    status = str(row.get("POPULATION_STATUS", "UNKNOWN")).upper()

    if af is not None:
        if af == 0:
            return "NOT_OBSERVED_IN_NEEDLR_CONTROLS"
        if af <= 0.001:
            return "VERY_RARE_NEEDLR_LE_0.001"
        if af <= 0.01:
            return "RARE_NEEDLR_LE_0.01"
        if af < 0.05:
            return "COMMON_NEEDLR_GT_0.01"
        return "VERY_COMMON_NEEDLR_GE_0.05"

    if status == "NO_POPULATION_MATCH":
        return "NO_NEEDLR_MATCH_AF_UNKNOWN"
    if status == "NOT_EVALUABLE_GE_10MB":
        return "NEEDLR_NOT_EVALUABLE_GE_10MB"
    if status == "NOT_EVALUABLE_BREAKEND":
        return "NEEDLR_NOT_EVALUABLE_BREAKEND"
    return "NEEDLR_FREQUENCY_UNKNOWN"


def benign_afmax_class(row):
    af = number(
        row.get(
            "ANNOTSV_BENIGN_AFMAX",
            row.get("BENIGN_DB_AFMAX"),
        )
    )
    if af is None:
        return "NO_BENIGN_REGION_AFMAX"
    if af == 0:
        return "BENIGN_REGION_AFMAX_ZERO"
    if af <= 0.001:
        return "BENIGN_REGION_AFMAX_LE_0.001"
    if af <= 0.01:
        return "BENIGN_REGION_AFMAX_LE_0.01"
    if af < 0.05:
        return "BENIGN_REGION_AFMAX_GT_0.01"
    return "BENIGN_REGION_AFMAX_GE_0.05"


def gnomad_exact_class(row):
    match = str(row.get("GNOMAD_SV_EXACT_MATCH", ".")).upper()
    af = number(row.get("GNOMAD_SV_AF"))
    if match == "RESOURCE_NOT_CONFIGURED":
        return "GNOMAD_SV_RESOURCE_NOT_CONFIGURED"
    if match != "YES":
        return "NO_EXACT_GNOMAD_SV_MATCH"
    if af is None:
        return "EXACT_GNOMAD_MATCH_AF_UNAVAILABLE"
    if af == 0:
        return "EXACT_GNOMAD_AF_ZERO"
    if af <= 0.001:
        return "EXACT_GNOMAD_VERY_RARE_LE_0.001"
    if af <= 0.01:
        return "EXACT_GNOMAD_RARE_LE_0.01"
    if af < 0.05:
        return "EXACT_GNOMAD_COMMON_GT_0.01"
    return "EXACT_GNOMAD_VERY_COMMON_GE_0.05"


def gnomad_context(row):
    overlap = str(row.get("GNOMAD_SV_OVERLAP", ".")).upper()
    source = str(
        row.get(
            "ANNOTSV_BENIGN_DB_SOURCE",
            row.get("BENIGN_DB_SOURCE", "."),
        )
    )

    source_mentions = "GNOMAD" in source.upper()
    if overlap == "YES" or source_mentions:
        return "GNOMAD_INCLUDED_IN_ANNOTSV_BENIGN_OVERLAP"
    if overlap == "NOT_REPORTED":
        return "GNOMAD_NOT_REPORTED_IN_ANNOTSV_BENIGN_OVERLAP"
    if overlap == "NOT_APPLICABLE":
        return "GNOMAD_OVERLAP_NOT_APPLICABLE_FOR_SVTYPE"
    return "GNOMAD_OVERLAP_CONTEXT_UNKNOWN"


def population_class(row):
    """Combine needLR and exact-site gnomAD frequency without hiding conflicts."""
    nclass = needlr_class(row)
    gclass = gnomad_exact_class(row)
    bclass = benign_afmax_class(row)

    needlr_common = nclass in {
        "VERY_COMMON_NEEDLR_GE_0.05",
        "COMMON_NEEDLR_GT_0.01",
    }
    needlr_low = nclass in {
        "NOT_OBSERVED_IN_NEEDLR_CONTROLS",
        "VERY_RARE_NEEDLR_LE_0.001",
        "RARE_NEEDLR_LE_0.01",
    }

    gnomad_common = gclass in {
        "EXACT_GNOMAD_COMMON_GT_0.01",
        "EXACT_GNOMAD_VERY_COMMON_GE_0.05",
    }
    gnomad_low = gclass in {
        "EXACT_GNOMAD_AF_ZERO",
        "EXACT_GNOMAD_VERY_RARE_LE_0.001",
        "EXACT_GNOMAD_RARE_LE_0.01",
    }

    if (needlr_common and gnomad_low) or (needlr_low and gnomad_common):
        return "FREQUENCY_SOURCES_CONFLICT"

    if needlr_common or gnomad_common:
        if needlr_common and gnomad_common:
            return "COMMON_SUPPORTED_BY_NEEDLR_AND_GNOMAD"
        if gnomad_common:
            return "COMMON_BY_EXACT_GNOMAD"
        return "COMMON_BY_NEEDLR"

    if needlr_low or gnomad_low:
        if needlr_low and gnomad_low:
            return "LOW_FREQUENCY_SUPPORTED_BY_NEEDLR_AND_GNOMAD"
        if gnomad_low:
            return "LOW_FREQUENCY_BY_EXACT_GNOMAD"
        return "LOW_FREQUENCY_BY_NEEDLR"

    if bclass in {
        "BENIGN_REGION_AFMAX_GE_0.05",
        "BENIGN_REGION_AFMAX_GT_0.01",
    }:
        return "COMMON_BENIGN_REGION_OVERLAP_CONTEXT"

    if bclass in {
        "BENIGN_REGION_AFMAX_ZERO",
        "BENIGN_REGION_AFMAX_LE_0.001",
        "BENIGN_REGION_AFMAX_LE_0.01",
    }:
        return "LOW_AF_BENIGN_REGION_OVERLAP_CONTEXT"

    if "NOT_EVALUABLE" in nclass:
        return "POPULATION_AF_NOT_EVALUABLE"
    if nclass == "NO_NEEDLR_MATCH_AF_UNKNOWN":
        return "NO_NEEDLR_MATCH_AF_UNKNOWN"
    if gclass == "GNOMAD_SV_RESOURCE_NOT_CONFIGURED":
        return "GNOMAD_RESOURCE_NOT_CONFIGURED"
    return "POPULATION_FREQUENCY_UNKNOWN"


def population_interpretation(row):
    cls = population_class(row)

    if cls == "FREQUENCY_SOURCES_CONFLICT":
        return "NEEDLR_AND_GNOMAD_FREQUENCY_EVIDENCE_DISAGREE_REVIEW_EVENT_MATCHING"

    if cls in {
        "COMMON_SUPPORTED_BY_NEEDLR_AND_GNOMAD",
        "COMMON_BY_EXACT_GNOMAD",
        "COMMON_BY_NEEDLR",
    }:
        return "COMMON_FREQUENCY_WEAKENS_CANDIDACY_FOR_A_HIGHLY_PENETRANT_RARE_MENDELIAN_ALLELE"

    if cls in {
        "LOW_FREQUENCY_SUPPORTED_BY_NEEDLR_AND_GNOMAD",
        "LOW_FREQUENCY_BY_EXACT_GNOMAD",
        "LOW_FREQUENCY_BY_NEEDLR",
    }:
        return "LOW_FREQUENCY_SUPPORTS_RARITY_ONLY_NOT_PATHOGENICITY"

    if cls == "COMMON_BENIGN_REGION_OVERLAP_CONTEXT":
        return "COMMON_OVERLAPPING_BENIGN_SV_CONTEXT_REQUIRES_ALLELE_EQUIVALENCE_REVIEW"

    if cls == "LOW_AF_BENIGN_REGION_OVERLAP_CONTEXT":
        return "LOW_AF_OVERLAPPING_BENIGN_SV_CONTEXT_DOES_NOT_ESTABLISH_BENIGNITY"

    if cls == "NO_NEEDLR_MATCH_AF_UNKNOWN":
        return "NO_NEEDLR_MATCH_IS_NOT_EQUIVALENT_TO_AF_ZERO"

    if cls == "POPULATION_AF_NOT_EVALUABLE":
        return "NO_VALID_NEEDLR_AF_FOR_THIS_EVENT_CLASS"

    if cls == "GNOMAD_RESOURCE_NOT_CONFIGURED":
        return "NO_GNOMAD_SPECIFIC_SITE_FREQUENCY_AVAILABLE"

    return "INSUFFICIENT_POPULATION_FREQUENCY_EVIDENCE"


def effect_group(effect, svtype):
    effect = str(effect or "").upper()
    svtype = str(svtype or "").upper()
    if "BREAKPOINT_IN_" in effect or "TWO_BREAKPOINTS_IN_GENE" in effect:
        return "DIRECT_BREAKPOINT"
    if "WHOLE_GENE_DELETION" in effect or "PARTIAL_GENE_DELETION" in effect:
        return "COPY_LOSS_GEOMETRY"
    if "WHOLE_GENE_DUPLICATION" in effect or "PARTIAL_GENE_DUPLICATION" in effect:
        return "COPY_GAIN_GEOMETRY"
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
        pd.DataFrame(
            columns=[
                "GENE", "GENE_RELEVANCE_SCORE", "PANEL_STATUS", "SV_COUNT",
                "LOW_FREQUENCY_NEEDLR_COUNT", "COMMON_NEEDLR_COUNT",
                "GNOMAD_EXACT_MATCH_COUNT", "LOW_FREQUENCY_GNOMAD_EXACT_COUNT",
                "COMMON_GNOMAD_EXACT_COUNT", "GNOMAD_BENIGN_OVERLAP_COUNT",
                "COMMON_BENIGN_REGION_CONTEXT_COUNT",
                "POPULATION_CLASSES", "SV_EFFECTS",
            ]
        ).to_csv(args.gene_summary, sep="\t", index=False)
        return

    out = df.copy()
    out["NEEDLR_FREQUENCY_CLASS"] = out.apply(needlr_class, axis=1)
    out["ANNOTSV_BENIGN_AFMAX_CLASS"] = out.apply(benign_afmax_class, axis=1)
    out["GNOMAD_EXACT_AF_CLASS"] = out.apply(gnomad_exact_class, axis=1)
    out["GNOMAD_CONTEXT_CLASS"] = out.apply(gnomad_context, axis=1)
    out["POPULATION_CLASS"] = out.apply(population_class, axis=1)
    out["POPULATION_INTERPRETATION"] = out.apply(
        population_interpretation,
        axis=1,
    )
    out["SV_EFFECT_GROUP"] = [
        effect_group(effect, svtype)
        for effect, svtype in zip(out["SV_GENE_EFFECT"], out["SVTYPE"])
    ]
    out["POPULATION_EVIDENCE_SCOPE"] = (
        "needLR is provisional coordinate-compatible long-read control-frequency evidence. "
        "gnomAD-SV exact AF is reported only for the conservative exact "
        "coordinate/type match defined by this workflow. AnnotSV benign AFmax "
        "is overlap-region evidence and may combine multiple benign resources."
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    unique = out.drop_duplicates(["GENE", "SV_ID"])
    summary_rows = []
    for gene, group in unique.groupby("GENE", sort=False):
        classes = group["POPULATION_CLASS"].value_counts().to_dict()
        effects = group["SV_EFFECT_GROUP"].value_counts().to_dict()
        score = pd.to_numeric(
            group["GENE_RELEVANCE_SCORE"],
            errors="coerce",
        ).max()

        low_needlr = group["NEEDLR_FREQUENCY_CLASS"].isin({
            "NOT_OBSERVED_IN_NEEDLR_CONTROLS",
            "VERY_RARE_NEEDLR_LE_0.001",
            "RARE_NEEDLR_LE_0.01",
        })
        common_needlr = group["NEEDLR_FREQUENCY_CLASS"].isin({
            "COMMON_NEEDLR_GT_0.01",
            "VERY_COMMON_NEEDLR_GE_0.05",
        })
        exact_gnomad = group["GNOMAD_SV_EXACT_MATCH"].fillna("").eq("YES")
        low_gnomad = group["GNOMAD_EXACT_AF_CLASS"].isin({
            "EXACT_GNOMAD_AF_ZERO",
            "EXACT_GNOMAD_VERY_RARE_LE_0.001",
            "EXACT_GNOMAD_RARE_LE_0.01",
        })
        common_gnomad = group["GNOMAD_EXACT_AF_CLASS"].isin({
            "EXACT_GNOMAD_COMMON_GT_0.01",
            "EXACT_GNOMAD_VERY_COMMON_GE_0.05",
        })

        summary_rows.append({
            "GENE": gene,
            "GENE_RELEVANCE_SCORE": score if pd.notna(score) else ".",
            "PANEL_STATUS": (
                "PANEL_GENE"
                if (group["PANEL_STATUS"] == "PANEL_GENE").any()
                else "NON_PANEL"
            ),
            "SV_COUNT": int(group["SV_ID"].nunique()),
            "LOW_FREQUENCY_NEEDLR_COUNT": int(
                group.loc[low_needlr, "SV_ID"].nunique()
            ),
            "COMMON_NEEDLR_COUNT": int(
                group.loc[common_needlr, "SV_ID"].nunique()
            ),
            "GNOMAD_EXACT_MATCH_COUNT": int(
                group.loc[exact_gnomad, "SV_ID"].nunique()
            ),
            "LOW_FREQUENCY_GNOMAD_EXACT_COUNT": int(
                group.loc[low_gnomad, "SV_ID"].nunique()
            ),
            "COMMON_GNOMAD_EXACT_COUNT": int(
                group.loc[common_gnomad, "SV_ID"].nunique()
            ),
            "GNOMAD_BENIGN_OVERLAP_COUNT": int(
                group.loc[
                    group["GNOMAD_CONTEXT_CLASS"].eq(
                        "GNOMAD_INCLUDED_IN_ANNOTSV_BENIGN_OVERLAP"
                    ),
                    "SV_ID",
                ].nunique()
            ),
            "COMMON_BENIGN_REGION_CONTEXT_COUNT": int(
                group.loc[
                    group["POPULATION_CLASS"].eq(
                        "COMMON_BENIGN_REGION_OVERLAP_CONTEXT"
                    ),
                    "SV_ID",
                ].nunique()
            ),
            "POPULATION_CLASSES": ";".join(
                f"{k}={v}" for k, v in sorted(classes.items())
            ),
            "SV_EFFECTS": ";".join(
                f"{k}={v}" for k, v in sorted(effects.items())
            ),
        })

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["_score"] = pd.to_numeric(
            summary["GENE_RELEVANCE_SCORE"],
            errors="coerce",
        ).fillna(0)
        summary = summary.sort_values(
            [
                "PANEL_STATUS",
                "LOW_FREQUENCY_NEEDLR_COUNT",
                "_score",
                "COMMON_NEEDLR_COUNT",
                "SV_COUNT",
                "GENE",
            ],
            ascending=[True, False, False, True, False, True],
        ).drop(columns="_score")

    summary.to_csv(args.gene_summary, sep="\t", index=False)

    print(
        f"[OK] population_rows={len(out)} genes={out['GENE'].nunique()} "
        f"output={output} gene_summary={args.gene_summary}"
    )


if __name__ == "__main__":
    main()
