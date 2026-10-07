#!/usr/bin/env python3
"""Rank mitochondrial candidate genes affected by structural variants.

The ranking is intentionally split by encoding genome:
  1. nuclear-encoded mitochondrial genes from MitoCarta3.0;
  2. genes encoded by the mitochondrial genome (protein, rRNA and tRNA genes).

This is a research-prioritization layer, not a pathogenicity classifier.
Mechanism/directness, optic-neuropathy context, technical evidence and
population evidence remain explicit columns so the ordering is auditable.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from ranking_common import (
    GENE_TIER_PRIORITY,
    MISSING,
    POPULATION_PRIORITY,
    TECHNICAL_PRIORITY,
    first_existing,
    mechanism_inheritance_summary,
    number,
    population_summary,
    relationship_priority,
    technical_summary,
)

MTDNA_PROTEIN_GENES = {
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
}
MTDNA_RRNA_GENES = {"MT-RNR1", "MT-RNR2"}


def normalize_chrom(value) -> str:
    text = str(value or "").strip()
    if text.upper() in MISSING:
        return "."
    if text in {"MT", "M", "chrMT"}:
        return "chrM"
    return text if text.startswith("chr") else "chr" + text


def is_yes(value) -> bool:
    return str(value or "").strip().upper() in {"YES", "TRUE", "1", "PANEL_GENE"}


def mt_gene_class(gene: str) -> str:
    gene = gene.upper()
    if gene in MTDNA_PROTEIN_GENES:
        return "MTDNA_PROTEIN_CODING"
    if gene in MTDNA_RRNA_GENES:
        return "MTDNA_RRNA"
    if gene.startswith("MT-T"):
        return "MTDNA_TRNA"
    return "MTDNA_OTHER_GENE"


def mt_gene_function(gene: str) -> str:
    """Broad biological role for genes physically encoded by mtDNA."""
    gene = gene.upper()
    if gene in {"MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6"}:
        return "OXPHOS_COMPLEX_I"
    if gene == "MT-CYB":
        return "OXPHOS_COMPLEX_III"
    if gene in {"MT-CO1", "MT-CO2", "MT-CO3"}:
        return "OXPHOS_COMPLEX_IV"
    if gene in {"MT-ATP6", "MT-ATP8"}:
        return "OXPHOS_COMPLEX_V"
    if gene in MTDNA_RRNA_GENES:
        return "MITOCHONDRIAL_TRANSLATION_RRNA"
    if gene.startswith("MT-T"):
        return "MITOCHONDRIAL_TRANSLATION_TRNA"
    return "MTDNA_OTHER_FUNCTION"


def relationship_scope(value: str) -> str:
    rank = relationship_priority(value)
    return {
        3: "DIRECT_STRUCTURAL_EFFECT",
        2: "PROXIMAL_OR_BREAKPOINT_CONTEXT",
        1: "INTERVAL_OR_INVERSION_SPANNED_CONTEXT",
        0: "UNRESOLVED_CONTEXT",
    }.get(rank, "UNRESOLVED_CONTEXT")


def read_gene_set(path: str | None) -> set[str]:
    if not path:
        return set()
    genes = set()
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            gene = line.split("\t", 1)[0].split(",", 1)[0].strip().upper()
            if gene and gene not in MISSING:
                genes.add(gene)
    return genes


def read_gene_bed(path: str) -> pd.DataFrame:
    bed = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        dtype=str,
        low_memory=False,
    )
    if bed.shape[1] < 4:
        raise ValueError("Gene BED must contain chromosome, start, end and gene columns.")
    bed = bed.iloc[:, :4].copy()
    bed.columns = ["CHROM", "GENE_START", "GENE_END", "GENE"]
    bed["CHROM"] = bed["CHROM"].map(normalize_chrom)
    bed["GENE_START"] = pd.to_numeric(bed["GENE_START"], errors="coerce")
    bed["GENE_END"] = pd.to_numeric(bed["GENE_END"], errors="coerce")
    bed["GENE"] = bed["GENE"].fillna("").astype(str).str.upper().str.strip()
    return bed.dropna(subset=["GENE_START", "GENE_END"]).drop_duplicates(
        ["CHROM", "GENE_START", "GENE_END", "GENE"]
    )


def classify_mt_overlap(event: pd.Series, gene_start: int, gene_end: int) -> str:
    svtype = str(event.get("SVTYPE", "") or "").upper()
    start = number(event.get("START"))
    end = number(event.get("END"))
    pos2 = number(event.get("POS2"))

    if start is None:
        return "GENE_RELATIONSHIP_UNRESOLVED"

    # Convert 1-based VCF coordinates to BED-like zero-based positions.
    start0 = max(0, int(start) - 1)
    end0 = int(end) if end is not None else start0 + 1

    if svtype in {"DEL", "DUP", "CNV"}:
        lo, hi = sorted((start0, max(start0 + 1, end0)))
        if lo <= gene_start and hi >= gene_end:
            return "WHOLE_GENE_DOSAGE_CONTEXT"
        if hi > gene_start and lo < gene_end:
            return "PARTIAL_GENE_OVERLAP"
        return "GENE_RELATIONSHIP_UNRESOLVED"

    if svtype == "INS":
        return (
            "INSERTION_WITHIN_TRANSCRIPT"
            if gene_start <= start0 < gene_end
            else "GENE_RELATIONSHIP_UNRESOLVED"
        )

    if svtype in {"INV", "BND", "TRA"}:
        bps = [start0]
        if pos2 is not None:
            bps.append(max(0, int(pos2) - 1))
        elif end is not None:
            bps.append(max(0, int(end) - 1))

        if any(gene_start <= bp < gene_end for bp in bps):
            return "BREAKPOINT_WITHIN_TRANSCRIPT"

        if svtype == "INV" and len(bps) >= 2:
            lo, hi = sorted(bps[:2])
            if lo <= gene_start and hi >= gene_end:
                return "INVERSION_SPANS_INTACT_GENE"
            if hi > gene_start and lo < gene_end:
                return "INTERVAL_CONTEXT_ONLY"

    return "GENE_RELATIONSHIP_UNRESOLVED"


def build_mtdna_overlap_rows(
    df: pd.DataFrame,
    gene_bed: pd.DataFrame,
    panel_genes: set[str],
) -> pd.DataFrame:
    id_col = first_existing(df, ["SV_ID", "ID"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    if id_col is None or chrom_col is None:
        return pd.DataFrame()

    mt_genes = gene_bed[
        gene_bed["CHROM"].eq("chrM")
        & gene_bed["GENE"].str.startswith("MT-")
    ].copy()
    if mt_genes.empty:
        return pd.DataFrame()

    event_cols = []
    for col in [
        id_col, chrom_col, "START", "END", "POS2", "CHR2", "SVTYPE", "SVLEN",
        "CALLER_COUNT", "SUPP", "EVENT_TECHNICAL_TIER",
        "EVENT_TECHNICAL_REVIEW", "EVENT_POPULATION_TIER",
        "GNOMAD_SV_AF", "GNOMAD_SV_EXACT_MATCH", "NEEDLR_AF",
    ]:
        if col in df.columns and col not in event_cols:
            event_cols.append(col)
    events = df[event_cols].drop_duplicates(id_col).copy()
    events["_CHROM_NORM"] = events[chrom_col].map(normalize_chrom)
    events = events[events["_CHROM_NORM"].eq("chrM")]

    rows = []
    for _, event in events.iterrows():
        start = number(event.get("START"))
        end = number(event.get("END"))
        if start is None:
            continue

        start0 = max(0, int(start) - 1)
        if str(event.get("SVTYPE", "")).upper() == "INS":
            end0 = start0 + 1
        else:
            end0 = int(end) if end is not None else start0 + 1
            if end0 <= start0:
                end0 = start0 + 1

        overlapping = mt_genes[
            mt_genes["GENE_END"].gt(start0)
            & mt_genes["GENE_START"].lt(end0)
        ].copy()

        # For breakpoint events, also retain a gene hit by either breakpoint.
        if str(event.get("SVTYPE", "")).upper() in {"INV", "BND", "TRA"}:
            bp_positions = [start0]
            pos2 = number(event.get("POS2"))
            if pos2 is not None:
                bp_positions.append(max(0, int(pos2) - 1))
            for bp in bp_positions:
                extra = mt_genes[
                    mt_genes["GENE_START"].le(bp)
                    & mt_genes["GENE_END"].gt(bp)
                ]
                overlapping = pd.concat([overlapping, extra], ignore_index=True)

        if overlapping.empty:
            continue

        overlapping = overlapping.drop_duplicates(["GENE_START", "GENE_END", "GENE"])
        for _, gene_rec in overlapping.iterrows():
            row = event.to_dict()
            row["GENE"] = gene_rec["GENE"]
            row["SV_GENE_RELATIONSHIP"] = classify_mt_overlap(
                event,
                int(gene_rec["GENE_START"]),
                int(gene_rec["GENE_END"]),
            )
            row["MITOCARTA_STATUS"] = (
                "YES" if gene_rec["GENE"] in MTDNA_PROTEIN_GENES else "NO"
            )
            row["MITOCARTA_ENCODING"] = "MTDNA_ENCODED_GENE"
            row["MITOCARTA_MITOPATHWAYS"] = "."
            row["MITOCARTA_TOP_LEVEL_PATHWAYS"] = "."
            row["MITOCARTA_SUBCOMPARTMENT"] = "."
            is_panel = gene_rec["GENE"] in panel_genes
            row["MITO_ON_CONTEXT"] = "NO"
            row["MITO_ON_ANCHOR_HPO_COUNT"] = "0"
            row["PANEL_STATUS"] = "PANEL_GENE" if is_panel else "NONPANEL_GENE"
            # mtDNA panel membership is not converted into gene-relevance
            # points. mtDNA rows are ranked by their own SV relationship,
            # technical and population context.
            row["GENE_RELEVANCE_TIER"] = "LIMITED"
            row["GENE_RELEVANCE_DISPLAY_SCORE"] = "0"
            row["EVENT_GENE_RELEVANCE_SCORE"] = "0"
            rows.append(row)

    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, help="Final integrated SV-gene table")
    p.add_argument("--gene-bed", required=True, help="BED4+ gene annotation including chrM")
    p.add_argument(
        "--panel",
        default=None,
        help="Optional optic-neuropathy candidate gene list; used for context, never as a discovery filter.",
    )
    p.add_argument("--output", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    id_col = first_existing(df, ["SV_ID", "ID"])
    if gene_col is None or id_col is None:
        raise ValueError("Integrated table needs SV_ID and gene columns.")

    work = df.copy()
    work["_GENE"] = work[gene_col].fillna(".").astype(str).str.upper().str.strip()

    nuclear = work[
        work.get("MITOCARTA_ENCODING", pd.Series("", index=work.index))
        .fillna("")
        .astype(str)
        .str.upper()
        .eq("NUCLEAR_MITOCHONDRIAL_GENE")
    ].copy()
    nuclear["_ENCODING_GENOME"] = "NUCLEAR"
    nuclear["_GENE_CLASS"] = "NUCLEAR_ENCODED_MITOCHONDRIAL"
    nuclear["_MTDNA_FUNCTION"] = "."

    mt_existing = work[
        work.get("MITOCARTA_ENCODING", pd.Series("", index=work.index))
        .fillna("")
        .astype(str)
        .str.upper()
        .eq("MTDNA_ENCODED_GENE")
    ].copy()
    if not mt_existing.empty:
        mt_existing["_ENCODING_GENOME"] = "MTDNA"
        mt_existing["_GENE_CLASS"] = mt_existing["_GENE"].map(mt_gene_class)
        mt_existing["_MTDNA_FUNCTION"] = mt_existing["_GENE"].map(mt_gene_function)

    bed = read_gene_bed(args.gene_bed)
    panel_genes = read_gene_set(args.panel)
    mt_overlap = build_mtdna_overlap_rows(work, bed, panel_genes)
    if not mt_overlap.empty:
        mt_overlap["_GENE"] = mt_overlap["GENE"].astype(str).str.upper().str.strip()
        mt_overlap["_ENCODING_GENOME"] = "MTDNA"
        mt_overlap["_GENE_CLASS"] = mt_overlap["_GENE"].map(mt_gene_class)
        mt_overlap["_MTDNA_FUNCTION"] = mt_overlap["_GENE"].map(mt_gene_function)

    mt = pd.concat([mt_existing, mt_overlap], ignore_index=True, sort=False)
    if not mt.empty:
        mt = mt.drop_duplicates([id_col, "_GENE"], keep="first")

    candidates = pd.concat([nuclear, mt], ignore_index=True, sort=False)
    if candidates.empty:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=[
            "GENE", "ENCODING_GENOME", "GENE_CLASS", "MTDNA_FUNCTION",
            "MITO_RANK_WITHIN_ENCODING", "MITO_PRIORITY_TIER",
            "UNIQUE_SVS", "DIRECT_SVS", "PROXIMAL_SVS",
            "CONTEXT_SVS", "BEST_SV_ID", "BEST_SVTYPE", "BEST_SV_GENE_RELATIONSHIP",
        ]).to_csv(output, sep="\t", index=False)
        print(f"[OK] mitochondrial_ranked_genes=0 output={output}")
        return

    if "SV_GENE_RELATIONSHIP" not in candidates.columns:
        candidates["SV_GENE_RELATIONSHIP"] = "GENE_RELATIONSHIP_UNRESOLVED"

    # Apply the same evidence axes used by the general event ranker.
    mechanism_rows = [
        mechanism_inheritance_summary(row)
        for _, row in candidates.iterrows()
    ]
    technical_rows = [
        technical_summary(row)
        for _, row in candidates.iterrows()
    ]
    population_rows = [
        population_summary(row)
        for _, row in candidates.iterrows()
    ]

    candidates["_GENE_TIER_RANK"] = (
        candidates.get(
            "GENE_RELEVANCE_TIER",
            pd.Series("LIMITED", index=candidates.index),
        )
        .fillna("LIMITED")
        .astype(str)
        .str.upper()
        .map(GENE_TIER_PRIORITY)
        .fillna(0)
    )
    candidates["_GENE_RELEVANCE"] = pd.to_numeric(
        candidates.get(
            "GENE_RELEVANCE_DISPLAY_SCORE",
            candidates.get(
                "EVENT_GENE_RELEVANCE_SCORE",
                pd.Series(0, index=candidates.index),
            ),
        ),
        errors="coerce",
    ).fillna(0)

    candidates["_MECHANISM_CLASS"] = [x["category"] for x in mechanism_rows]
    candidates["_INHERITANCE_CLASS"] = [
        x["inheritance_class"] for x in mechanism_rows
    ]
    candidates["_MECHANISM_DETAIL"] = [x["detail"] for x in mechanism_rows]

    # Nuclear-encoded genes use the shared inheritance/mechanism rank. mtDNA
    # genes are a separate biological class and are ordered by direct/proximal/
    # interval relationship while remaining explicitly marked for specialized
    # heteroplasmy/mtDNA interpretation.
    candidates["_MECHANISM_RANK"] = [
        (
            relationship_priority(row.get("SV_GENE_RELATIONSHIP"))
            if str(row.get("_ENCODING_GENOME", "")).upper() == "MTDNA"
            else int(summary["priority"])
        )
        for summary, (_, row) in zip(mechanism_rows, candidates.iterrows())
    ]

    candidates["_TECHNICAL_TIER"] = [x["tier"] for x in technical_rows]
    candidates["_TECHNICAL_RANK"] = [
        TECHNICAL_PRIORITY.get(x["tier"], 0) for x in technical_rows
    ]
    candidates["_POPULATION_TIER"] = [x["tier"] for x in population_rows]
    candidates["_POPULATION_RANK"] = [
        POPULATION_PRIORITY.get(x["tier"], 1) for x in population_rows
    ]
    candidates["_MAX_AF"] = [x["max_af"] for x in population_rows]

    rows = []
    for (encoding, gene), group in candidates.groupby(
        ["_ENCODING_GENOME", "_GENE"],
        sort=False,
    ):
        group = group.copy().sort_values(
            [
                "_GENE_TIER_RANK",
                "_MECHANISM_RANK",
                "_TECHNICAL_RANK",
                "_POPULATION_RANK",
                "_GENE_RELEVANCE",
                id_col,
            ],
            ascending=[False, False, False, False, False, True],
        )
        best = group.iloc[0]

        relation_rank = group["SV_GENE_RELATIONSHIP"].map(
            relationship_priority
        )
        direct = relation_rank.eq(3)
        proximal = relation_rank.eq(2)
        context = relation_rank.eq(1)

        pathways = sorted({
            item.strip()
            for raw in group.get(
                "MITOCARTA_MITOPATHWAYS",
                pd.Series(".", index=group.index),
            ).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })
        top_pathways = sorted({
            item.strip()
            for raw in group.get(
                "MITOCARTA_TOP_LEVEL_PATHWAYS",
                pd.Series(".", index=group.index),
            ).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })
        compartments = sorted({
            item.strip()
            for raw in group.get(
                "MITOCARTA_SUBCOMPARTMENT",
                pd.Series(".", index=group.index),
            ).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })

        mech_rank = int(best["_MECHANISM_RANK"])
        gene_tier = str(best.get("GENE_RELEVANCE_TIER", "LIMITED")).upper()
        if mech_rank >= 4 and gene_tier in {"HIGH", "MODERATE"}:
            priority = "TIER_1_MECHANISM_AND_GENE_RELEVANCE"
        elif mech_rank >= 3:
            priority = "TIER_2_DIRECT_MITOCHONDRIAL_GENE"
        elif mech_rank >= 2:
            priority = "TIER_3_PROXIMAL_MITOCHONDRIAL_CONTEXT"
        elif mech_rank >= 1:
            priority = "TIER_4_INTERVAL_CONTEXT_ONLY"
        else:
            priority = "TIER_5_UNRESOLVED_OR_SPECIALIZED_REVIEW"

        panel_status = (
            "PANEL_GENE"
            if group.get(
                "PANEL_STATUS",
                pd.Series(".", index=group.index),
            )
            .fillna("")
            .astype(str)
            .str.upper()
            .eq("PANEL_GENE")
            .any()
            else "NONPANEL_GENE"
        )

        row = {
            "GENE": gene,
            "ENCODING_GENOME": encoding,
            "GENE_CLASS": str(best.get("_GENE_CLASS", ".")),
            "MTDNA_FUNCTION": str(best.get("_MTDNA_FUNCTION", ".")),
            "PANEL_STATUS": panel_status,
            "ON_PANEL_GENE": "YES" if panel_status == "PANEL_GENE" else "NO",
            "GENE_RELEVANCE_TIER": gene_tier,
            "MAX_GENE_RELEVANCE": round(
                float(group["_GENE_RELEVANCE"].max()), 3
            ),
            "MITO_PRIORITY_TIER": priority,
            "UNIQUE_SVS": int(group[id_col].nunique()),
            "DIRECT_SVS": int(group.loc[direct, id_col].nunique()),
            "PROXIMAL_SVS": int(group.loc[proximal, id_col].nunique()),
            "CONTEXT_SVS": int(group.loc[context, id_col].nunique()),
            "BEST_SV_ID": str(best.get(id_col, ".")),
            "BEST_SVTYPE": str(best.get("SVTYPE", ".")),
            "BEST_SV_GENE_RELATIONSHIP": str(
                best.get("SV_GENE_RELATIONSHIP", ".")
            ),
            "BEST_RELATIONSHIP_SCOPE": relationship_scope(
                best.get("SV_GENE_RELATIONSHIP", ".")
            ),
            "BEST_INHERITANCE_CLASS": str(
                best.get("_INHERITANCE_CLASS", "UNKNOWN")
            ),
            "BEST_INHERITANCE_MECHANISM_CLASS": str(
                best.get("_MECHANISM_CLASS", "UNRESOLVED")
            ),
            "BEST_INHERITANCE_MECHANISM_DETAIL": str(
                best.get("_MECHANISM_DETAIL", ".")
            ),
            "BEST_TECHNICAL_TIER": str(
                best.get("_TECHNICAL_TIER", ".")
            ),
            "BEST_POPULATION_TIER": str(
                best.get("_POPULATION_TIER", ".")
            ),
            "BEST_MAX_EXPLICIT_AF": str(best.get("_MAX_AF", ".")),
            "BEST_GNOMAD_SV_AF": str(best.get("GNOMAD_SV_AF", ".")),
            "BEST_NEEDLR_AF": str(best.get("NEEDLR_AF", ".")),
            "MITOCARTA_MAESTRO_SCORE": str(
                best.get("MITOCARTA_MAESTRO_SCORE", ".")
            ),
            "MITOCARTA_SUBCOMPARTMENT": (
                ";".join(compartments) if compartments else "."
            ),
            "MITOCARTA_TOP_LEVEL_PATHWAYS": (
                ";".join(top_pathways) if top_pathways else "."
            ),
            "MITOCARTA_MITOPATHWAYS": (
                ";".join(pathways) if pathways else "."
            ),
            "_GENE_TIER_RANK": int(best["_GENE_TIER_RANK"]),
            "_MECHANISM_RANK": mech_rank,
            "_TECHNICAL_RANK": int(best["_TECHNICAL_RANK"]),
            "_POPULATION_RANK": int(best["_POPULATION_RANK"]),
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    # Nuclear-encoded mitochondrial genes are the primary study focus. mtDNA
    # genes remain visible, but they form a clearly labelled secondary section
    # that requires dedicated mtDNA/heteroplasmy interpretation.
    encoding_priority = {"NUCLEAR": 0, "MTDNA": 1}
    out["_ENCODING_PRIORITY"] = (
        out["ENCODING_GENOME"].map(encoding_priority).fillna(9)
    )
    out["MITO_ANALYSIS_ROLE"] = out["ENCODING_GENOME"].map(
        {
            "NUCLEAR": "PRIMARY_NUCLEAR_MITOCHONDRIAL",
            "MTDNA": "SECONDARY_MTDNA",
        }
    ).fillna("SECONDARY_OTHER")

    # Panel membership is descriptive and does not change the biological
    # priority or the evidence score.
    out = out.sort_values(
        [
            "_ENCODING_PRIORITY",
            "_GENE_TIER_RANK",
            "_MECHANISM_RANK",
            "_TECHNICAL_RANK",
            "_POPULATION_RANK",
            "MAX_GENE_RELEVANCE",
            "UNIQUE_SVS",
            "GENE",
        ],
        ascending=[
            True, False, False, False, False, False, False, True
        ],
    )

    out["MITO_REVIEW_ORDER"] = range(1, len(out) + 1)

    out["MITO_RANK_WITHIN_ENCODING"] = (
        out.groupby("ENCODING_GENOME").cumcount() + 1
    )
    out["MITO_RANK_WITHIN_ENCODING_PANEL_STATUS"] = (
        out.groupby(["ENCODING_GENOME", "PANEL_STATUS"]).cumcount() + 1
    )
    out["MITO_RANKING_MODEL"] = (
        "sharedTieredRanking__encoding_panel__geneTier_"
        "inheritanceMechanism_technical_population__v2"
    )
    out["MITO_RANKING_INTERPRETATION"] = (
        "Research prioritization only. Nuclear-encoded mitochondrial genes "
        "reuse the same inheritance/mechanism, technical and population axes "
        "as the general SV-gene ranker. mtDNA genes are ranked separately by "
        "SV-gene relationship and require mtDNA/heteroplasmy-specific review. "
        "Panel and non-panel ranks are separate; panel membership is not a score."
    )

    out = out.drop(
        columns=[
            "_GENE_TIER_RANK",
            "_MECHANISM_RANK",
            "_TECHNICAL_RANK",
            "_POPULATION_RANK",
            "_ENCODING_PRIORITY",
        ],
        errors="ignore",
    )
    ordered = [
        "MITO_REVIEW_ORDER",
        "MITO_ANALYSIS_ROLE",
        "MITO_RANK_WITHIN_ENCODING",
        "MITO_RANK_WITHIN_ENCODING_PANEL_STATUS",
        "GENE", "ENCODING_GENOME", "GENE_CLASS", "MTDNA_FUNCTION",
        "PANEL_STATUS", "ON_PANEL_GENE", "GENE_RELEVANCE_TIER",
        "MITO_PRIORITY_TIER", "UNIQUE_SVS", "DIRECT_SVS", "PROXIMAL_SVS",
        "CONTEXT_SVS", "BEST_SV_ID", "BEST_SVTYPE",
        "BEST_SV_GENE_RELATIONSHIP", "BEST_RELATIONSHIP_SCOPE",
        "BEST_INHERITANCE_CLASS", "BEST_INHERITANCE_MECHANISM_CLASS",
        "BEST_INHERITANCE_MECHANISM_DETAIL", "MAX_GENE_RELEVANCE",
        "BEST_TECHNICAL_TIER", "BEST_POPULATION_TIER",
        "BEST_MAX_EXPLICIT_AF", "BEST_GNOMAD_SV_AF", "BEST_NEEDLR_AF",
        "MITOCARTA_MAESTRO_SCORE", "MITOCARTA_SUBCOMPARTMENT",
        "MITOCARTA_TOP_LEVEL_PATHWAYS", "MITOCARTA_MITOPATHWAYS",
        "MITO_RANKING_MODEL", "MITO_RANKING_INTERPRETATION",
    ]
    out = out[[col for col in ordered if col in out.columns]]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    counts = out["ENCODING_GENOME"].value_counts().to_dict()
    print(
        f"[OK] mitochondrial_ranked_genes={len(out)} "
        f"nuclear={counts.get('NUCLEAR', 0)} mtdna={counts.get('MTDNA', 0)} "
        f"output={output}"
    )


if __name__ == "__main__":
    main()
