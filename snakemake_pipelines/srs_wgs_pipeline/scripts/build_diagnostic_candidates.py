#!/usr/bin/env python3
"""Prioritize SRS SV/gene calls for diagnostic-support review, without classifying pathogenicity."""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def column(df, names):
    lower = {str(c).lower(): c for c in df.columns}
    return next((n if n in df.columns else lower[n.lower()] for n in names if n in df.columns or n.lower() in lower), None)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("--input", required=True)
    p.add_argument("--top-n", type=int, default=100); p.add_argument("--output", required=True); a = p.parse_args()
    df = pd.read_csv(a.input, sep="\t", dtype=str, low_memory=False)
    gene = column(df, ["GENE", "GENES", "Gene", "gene", "ANNotsv_Gene"])
    if gene is None: raise ValueError("Candidate input has no gene column")
    mito = df.get("MITOCARTA_ENCODING", pd.Series(".", index=df.index)).eq("NUCLEAR_MITOCHONDRIAL_GENE")
    mtdna = df.get("MITOCARTA_ENCODING", pd.Series(".", index=df.index)).eq("MTDNA_ENCODED_GENE")
    panel_text = df.get("PANEL_STATUS", df.get("panel_gene", pd.Series(".", index=df.index))).fillna("").str.upper()
    panel = panel_text.isin({"PANEL_GENE", "YES", "TRUE", "1"})
    df["DIAGNOSTIC_FOCUS"] = "GENOMEWIDE_OTHER"
    df.loc[mtdna, "DIAGNOSTIC_FOCUS"] = "MTDNA_SECONDARY"
    df.loc[panel, "DIAGNOSTIC_FOCUS"] = "KNOWN_DISEASE_GENE"
    df.loc[mito, "DIAGNOSTIC_FOCUS"] = "NUCLEAR_MITOCHONDRIAL_PRIMARY"
    order = {"NUCLEAR_MITOCHONDRIAL_PRIMARY": 0, "KNOWN_DISEASE_GENE": 1, "GENOMEWIDE_OTHER": 2, "MTDNA_SECONDARY": 3}
    df["_focus"] = df["DIAGNOSTIC_FOCUS"].map(order).fillna(9)
    score_col = column(df, ["ALLELE_RESEARCH_SCORE", "EVENT_GENE_RELEVANCE_SCORE", "INTEGRATED_DISCOVERY_SCORE", "PHENOTYPE_SCORE"])
    df["_score"] = pd.to_numeric(df[score_col], errors="coerce").fillna(0) if score_col else 0
    technical = df.get("TECHNICAL_SUPPORT", pd.Series("REVIEW_REQUIRED", index=df.index))
    tech_order = {"MULTI_CALLER_PLUS_GRIDSS": 0, "MULTI_CALLER": 1, "SINGLE_CALLER_PLUS_SUPPORTING_EVIDENCE": 2, "SINGLE_CALLER": 3, "CNVPYTOR_DEPTH_ONLY": 4, "REVIEW_REQUIRED": 5}
    df["_technical"] = technical.map(tech_order).fillna(9)
    df = df.sort_values(["_focus", "_score", "_technical"], ascending=[True, False, True]).head(a.top_n).copy()
    df.insert(0, "DIAGNOSTIC_REVIEW_RANK", range(1, len(df) + 1))
    df["DIAGNOSTIC_USE_NOTE"] = "Research/diagnostic-support candidate; orthogonal validation and clinical interpretation required."
    df.drop(columns=["_focus", "_score", "_technical"], inplace=True)
    out = Path(a.output); out.parent.mkdir(parents=True, exist_ok=True); df.to_csv(out, sep="\t", index=False)
    print(f"[OK] diagnostic_candidates={len(df)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
