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


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}

DIRECT_RELATIONSHIPS = {
    "WHOLE_GENE_DOSAGE_CONTEXT",
    "PARTIAL_GENE_OVERLAP",
    "INSERTION_WITHIN_TRANSCRIPT",
    "BREAKPOINT_WITHIN_TRANSCRIPT",
}

PROXIMAL_RELATIONSHIPS = {
    "INSERTION_PROXIMAL_TO_GENE",
    "BREAKPOINT_PROXIMAL_TO_GENE",
    "GENE_PROXIMAL_INTERVAL",
}

CONTEXT_RELATIONSHIPS = {
    "INVERSION_SPANS_INTACT_GENE",
    "INTERVAL_CONTEXT_ONLY",
}

MTDNA_PROTEIN_GENES = {
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
}
MTDNA_RRNA_GENES = {"MT-RNR1", "MT-RNR2"}


def first_existing(df: pd.DataFrame, names: list[str]) -> str | None:
    lower = {str(c).strip().lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def normalize_chrom(value) -> str:
    text = str(value or "").strip()
    if text.upper() in MISSING:
        return "."
    if text in {"MT", "M", "chrMT"}:
        return "chrM"
    return text if text.startswith("chr") else "chr" + text


def number(value):
    try:
        result = float(value)
        return result if np.isfinite(result) else None
    except Exception:
        return None


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


def relationship_rank(value: str) -> int:
    value = str(value or "").upper()
    if value in DIRECT_RELATIONSHIPS:
        return 4
    if value in PROXIMAL_RELATIONSHIPS:
        return 3
    if value in CONTEXT_RELATIONSHIPS:
        return 2
    if value not in MISSING:
        return 1
    return 0


def relationship_scope(value: str) -> str:
    rank = relationship_rank(value)
    return {
        4: "DIRECT_STRUCTURAL_EFFECT",
        3: "PROXIMAL_OR_BREAKPOINT_CONTEXT",
        2: "INTERVAL_OR_INVERSION_SPANNED_CONTEXT",
        1: "UNRESOLVED_CONTEXT",
        0: "UNRESOLVED_CONTEXT",
    }[rank]


def population_rank(row: pd.Series) -> int:
    afs = []
    for col in ("GNOMAD_SV_AF", "NEEDLR_AF"):
        if col in row.index:
            value = number(row.get(col))
            if value is not None:
                afs.append(value)
    if not afs:
        return 1
    return 2 if min(afs) <= 0.01 else 0


def technical_rank(row: pd.Series) -> tuple[int, int]:
    callers = number(row.get("CALLER_COUNT"))
    if callers is None:
        callers = number(row.get("SUPP")) or 0
    tier = 2 if callers >= 2 else (1 if callers >= 1 else 0)
    review = str(row.get("EVENT_TECHNICAL_REVIEW", "") or "").upper()
    clean = int(review in {"", ".", "NO_REVIEW_FLAG_FROM_CALLER_EVIDENCE"})
    return tier, clean


def gene_relevance(row: pd.Series) -> float:
    for col in (
        "EVENT_GENE_RELEVANCE_SCORE",
        "GENE_RELEVANCE_SCORE",
        "INTEGRATED_DISCOVERY_SCORE",
    ):
        if col in row.index:
            value = number(row.get(col))
            if value is not None:
                return value
    phenotype = number(row.get("PHENOTYPE_SCORE")) or 0.0
    disease = number(row.get("GENE_DISEASE_EVIDENCE_SCORE")) or 0.0
    return phenotype + disease


def disease_context_rank(row: pd.Series) -> int:
    panel = is_yes(row.get("PANEL_STATUS"))
    mito_on = is_yes(row.get("MITO_ON_CONTEXT"))
    anchor = number(row.get("MITO_ON_ANCHOR_HPO_COUNT")) or 0
    if panel and (mito_on or anchor > 0):
        return 3
    if panel:
        return 2
    if mito_on or anchor > 0:
        return 1
    return 0


def priority_tier(mechanism: int, disease: int) -> str:
    if mechanism >= 4 and disease > 0:
        return "TIER_1_DIRECT_WITH_ON_DISEASE_CONTEXT"
    if mechanism >= 4:
        return "TIER_2_DIRECT_MITOCHONDRIAL_GENE"
    if mechanism == 3 and disease > 0:
        return "TIER_2_PROXIMAL_WITH_ON_DISEASE_CONTEXT"
    if mechanism == 3:
        return "TIER_3_PROXIMAL_MITOCHONDRIAL_CONTEXT"
    if mechanism == 2:
        return "TIER_4_INTERVAL_CONTEXT_ONLY"
    return "TIER_5_UNRESOLVED_CONTEXT"


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
            row["MITO_ON_CONTEXT"] = "YES" if is_panel else "NO"
            row["MITO_ON_ANCHOR_HPO_COUNT"] = "0"
            row["PANEL_STATUS"] = "PANEL_GENE" if is_panel else "NON_PANEL"
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

    candidates["_MECHANISM_RANK"] = candidates["SV_GENE_RELATIONSHIP"].map(relationship_rank)
    candidates["_DISEASE_CONTEXT_RANK"] = candidates.apply(disease_context_rank, axis=1)
    candidates["_GENE_RELEVANCE"] = candidates.apply(gene_relevance, axis=1)
    candidates["_POPULATION_RANK"] = candidates.apply(population_rank, axis=1)
    tech = candidates.apply(technical_rank, axis=1)
    candidates["_TECHNICAL_RANK"] = [x[0] for x in tech]
    candidates["_TECHNICAL_CLEAN"] = [x[1] for x in tech]

    rows = []
    for (encoding, gene), group in candidates.groupby(["_ENCODING_GENOME", "_GENE"], sort=False):
        group = group.copy()
        group = group.sort_values(
            [
                "_MECHANISM_RANK", "_DISEASE_CONTEXT_RANK", "_GENE_RELEVANCE",
                "_TECHNICAL_CLEAN", "_TECHNICAL_RANK", "_POPULATION_RANK", id_col,
            ],
            ascending=[False, False, False, False, False, False, True],
        )
        best = group.iloc[0]

        relation = group["SV_GENE_RELATIONSHIP"].fillna(".").astype(str)
        direct = relation.map(relationship_rank).eq(4)
        proximal = relation.map(relationship_rank).eq(3)
        context = relation.map(relationship_rank).eq(2)

        pathways = sorted({
            item.strip()
            for raw in group.get("MITOCARTA_MITOPATHWAYS", pd.Series(".", index=group.index)).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })
        top_pathways = sorted({
            item.strip()
            for raw in group.get("MITOCARTA_TOP_LEVEL_PATHWAYS", pd.Series(".", index=group.index)).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })
        compartments = sorted({
            item.strip()
            for raw in group.get("MITOCARTA_SUBCOMPARTMENT", pd.Series(".", index=group.index)).fillna(".")
            for item in str(raw).split(";")
            if item.strip().upper() not in MISSING
        })

        priority = priority_tier(
            int(best["_MECHANISM_RANK"]),
            int(best["_DISEASE_CONTEXT_RANK"]),
        )
        tier_order = {
            "TIER_1_DIRECT_WITH_ON_DISEASE_CONTEXT": 1,
            "TIER_2_DIRECT_MITOCHONDRIAL_GENE": 2,
            "TIER_2_PROXIMAL_WITH_ON_DISEASE_CONTEXT": 2,
            "TIER_3_PROXIMAL_MITOCHONDRIAL_CONTEXT": 3,
            "TIER_4_INTERVAL_CONTEXT_ONLY": 4,
            "TIER_5_UNRESOLVED_CONTEXT": 5,
        }.get(priority, 5)

        row = {
            "GENE": gene,
            "ENCODING_GENOME": encoding,
            "GENE_CLASS": str(best.get("_GENE_CLASS", ".")),
            "MTDNA_FUNCTION": str(best.get("_MTDNA_FUNCTION", ".")),
            "MITO_PRIORITY_TIER": priority,
            "_PRIORITY_TIER_ORDER": tier_order,
            "UNIQUE_SVS": int(group[id_col].nunique()),
            "DIRECT_SVS": int(group.loc[direct, id_col].nunique()),
            "PROXIMAL_SVS": int(group.loc[proximal, id_col].nunique()),
            "CONTEXT_SVS": int(group.loc[context, id_col].nunique()),
            "BEST_SV_ID": str(best.get(id_col, ".")),
            "BEST_SVTYPE": str(best.get("SVTYPE", ".")),
            "BEST_SV_GENE_RELATIONSHIP": str(best.get("SV_GENE_RELATIONSHIP", ".")),
            "BEST_RELATIONSHIP_SCOPE": relationship_scope(best.get("SV_GENE_RELATIONSHIP", ".")),
            "ON_PANEL_GENE": "YES" if group.get("PANEL_STATUS", pd.Series(".", index=group.index)).map(is_yes).any() else "NO",
            "MITO_ON_CONTEXT": "YES" if group.get("MITO_ON_CONTEXT", pd.Series(".", index=group.index)).map(is_yes).any() else "NO",
            "MAX_ON_ANCHOR_HPO_COUNT": int(max([
                number(v) or 0 for v in group.get("MITO_ON_ANCHOR_HPO_COUNT", pd.Series(0, index=group.index))
            ])),
            "MAX_GENE_RELEVANCE": round(float(group["_GENE_RELEVANCE"].max()), 3),
            "BEST_TECHNICAL_TIER": str(best.get("EVENT_TECHNICAL_TIER", ".")),
            "BEST_TECHNICAL_REVIEW": str(best.get("EVENT_TECHNICAL_REVIEW", ".")),
            "BEST_POPULATION_TIER": str(best.get("EVENT_POPULATION_TIER", ".")),
            "BEST_GNOMAD_SV_AF": str(best.get("GNOMAD_SV_AF", ".")),
            "BEST_NEEDLR_AF": str(best.get("NEEDLR_AF", ".")),
            "MITOCARTA_MAESTRO_SCORE": str(best.get("MITOCARTA_MAESTRO_SCORE", ".")),
            "MITOCARTA_SUBCOMPARTMENT": ";".join(compartments) if compartments else ".",
            "MITOCARTA_TOP_LEVEL_PATHWAYS": ";".join(top_pathways) if top_pathways else ".",
            "MITOCARTA_MITOPATHWAYS": ";".join(pathways) if pathways else ".",
            "_MECHANISM_RANK": int(best["_MECHANISM_RANK"]),
            "_DISEASE_CONTEXT_RANK": int(best["_DISEASE_CONTEXT_RANK"]),
            "_TECHNICAL_CLEAN": int(best["_TECHNICAL_CLEAN"]),
            "_TECHNICAL_RANK": int(best["_TECHNICAL_RANK"]),
            "_POPULATION_RANK": int(best["_POPULATION_RANK"]),
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    out = out.sort_values(
        [
            "ENCODING_GENOME", "_PRIORITY_TIER_ORDER",
            "_DISEASE_CONTEXT_RANK", "MAX_GENE_RELEVANCE",
            "_MECHANISM_RANK", "_TECHNICAL_CLEAN", "_TECHNICAL_RANK",
            "_POPULATION_RANK", "UNIQUE_SVS", "GENE",
        ],
        ascending=[True, True, False, False, False, False, False, False, False, True],
    )
    out["MITO_RANK_WITHIN_ENCODING"] = out.groupby("ENCODING_GENOME").cumcount() + 1
    out["MITO_RANKING_MODEL"] = (
        "mechanism_directness__ON_disease_context__gene_relevance__"
        "technical_support__explicit_population_evidence__v1"
    )
    out["MITO_RANKING_INTERPRETATION"] = (
        "Research prioritization only. Nuclear-encoded mitochondrial genes and "
        "mtDNA-encoded genes are ranked separately. MitoCarta membership and "
        "pathway assignment indicate mitochondrial biology, not pathogenicity. "
        "Missing or no-match population evidence is neutral and is not treated "
        "as proof of rarity."
    )

    out = out.drop(columns=[
        "_PRIORITY_TIER_ORDER", "_MECHANISM_RANK",
        "_DISEASE_CONTEXT_RANK", "_TECHNICAL_CLEAN",
        "_TECHNICAL_RANK", "_POPULATION_RANK",
    ])
    ordered = [
        "MITO_RANK_WITHIN_ENCODING", "GENE", "ENCODING_GENOME", "GENE_CLASS",
        "MTDNA_FUNCTION", "MITO_PRIORITY_TIER", "UNIQUE_SVS", "DIRECT_SVS", "PROXIMAL_SVS",
        "CONTEXT_SVS", "BEST_SV_ID", "BEST_SVTYPE", "BEST_SV_GENE_RELATIONSHIP",
        "BEST_RELATIONSHIP_SCOPE", "ON_PANEL_GENE", "MITO_ON_CONTEXT",
        "MAX_ON_ANCHOR_HPO_COUNT", "MAX_GENE_RELEVANCE", "BEST_TECHNICAL_TIER",
        "BEST_TECHNICAL_REVIEW", "BEST_POPULATION_TIER", "BEST_GNOMAD_SV_AF",
        "BEST_NEEDLR_AF", "MITOCARTA_MAESTRO_SCORE", "MITOCARTA_SUBCOMPARTMENT",
        "MITOCARTA_TOP_LEVEL_PATHWAYS", "MITOCARTA_MITOPATHWAYS",
        "MITO_RANKING_MODEL", "MITO_RANKING_INTERPRETATION",
    ]
    out = out[[c for c in ordered if c in out.columns]]

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
