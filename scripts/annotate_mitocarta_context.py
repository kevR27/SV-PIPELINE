#!/usr/bin/env python3
"""Annotate SV-overlapping genes with MitoCarta3.0 mitochondrial context.

This is a gene-context annotation layer. It does not infer that a structural
variant is pathogenic merely because the affected gene is mitochondrial or
phenotype-associated.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MTDNA_PROTEIN_GENES = {
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
}


def first_existing(df, names):
    lower = {str(c).strip().lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def read_mitocarta(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in {".xls", ".xlsx"}:
        sheets = pd.read_excel(path, sheet_name=None, dtype=str)
        candidates = []
        for name, df in sheets.items():
            symbol = first_existing(df, ["Symbol", "Gene Symbol", "GeneSymbol"])
            pathway = first_existing(df, ["MitoPathways", "MitoPathway"])
            compartment = first_existing(df, ["Sub-compartment", "Subcompartment", "SubMitoLocalization"])
            if symbol is not None:
                score = int(pathway is not None) + int(compartment is not None)
                candidates.append((score, len(df), name, df))
        if not candidates:
            raise ValueError("No MitoCarta sheet with a Symbol column was found.")
        candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
        df = candidates[0][3].copy()
    else:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)

    symbol_col = first_existing(df, ["Symbol", "Gene Symbol", "GeneSymbol"])
    if symbol_col is None:
        raise ValueError("MitoCarta input needs a Symbol column.")

    pathway_col = first_existing(df, ["MitoPathways", "MitoPathway"])
    compartment_col = first_existing(df, ["Sub-compartment", "Subcompartment", "SubMitoLocalization"])
    maestro_col = first_existing(df, ["Maestro score", "Maestro Score", "MaestroScore"])
    evidence_col = first_existing(df, ["Evidence"])
    description_col = first_existing(df, ["Description"])

    keep = pd.DataFrame({
        "MITOCARTA_SYMBOL": df[symbol_col].fillna("").astype(str).str.upper().str.strip(),
        "MITOCARTA_DESCRIPTION": df[description_col].fillna(".").astype(str) if description_col else ".",
        "MITOCARTA_MAESTRO_SCORE": df[maestro_col].fillna(".").astype(str) if maestro_col else ".",
        "MITOCARTA_EVIDENCE": df[evidence_col].fillna(".").astype(str) if evidence_col else ".",
        "MITOCARTA_SUBCOMPARTMENT": df[compartment_col].fillna(".").astype(str) if compartment_col else ".",
        "MITOCARTA_MITOPATHWAYS": df[pathway_col].fillna(".").astype(str) if pathway_col else ".",
    })
    keep = keep[keep["MITOCARTA_SYMBOL"].ne("")].drop_duplicates("MITOCARTA_SYMBOL")
    return keep.set_index("MITOCARTA_SYMBOL")


def split_pathways(value):
    text = str(value or "").strip()
    if text in {"", ".", "nan", "None"}:
        return []
    return [x.strip() for x in text.replace("|", ";").split(";") if x.strip()]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--integrated", required=True)
    p.add_argument("--mitocarta", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    sv = pd.read_csv(args.integrated, sep="\t", dtype=str, low_memory=False)
    gene_col = first_existing(sv, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    pheno_col = first_existing(sv, ["PHENOTYPE_SCORE", "phenotype_score"])
    if gene_col is None:
        raise ValueError("Integrated table needs an overlapping-gene column.")

    mito = read_mitocarta(args.mitocarta)

    out = sv.copy()
    genes = out[gene_col].fillna(".").astype(str).str.upper().str.strip()
    phenotype = pd.to_numeric(out[pheno_col], errors="coerce").fillna(0) if pheno_col else pd.Series(0, index=out.index)

    defaults = {
        "MITOCARTA_STATUS": "NO",
        "MITOCARTA_ENCODING": "NOT_MITOCARTA",
        "MITOCARTA_DESCRIPTION": ".",
        "MITOCARTA_MAESTRO_SCORE": ".",
        "MITOCARTA_EVIDENCE": ".",
        "MITOCARTA_SUBCOMPARTMENT": ".",
        "MITOCARTA_MITOPATHWAYS": ".",
        "MITOCARTA_TOP_LEVEL_PATHWAYS": ".",
        "MITO_ON_CONTEXT": "NO",
        "MITO_ON_ASSOCIATION_CLASS": "NOT_MITOCARTA",
    }
    for col, value in defaults.items():
        out[col] = value

    matched = genes.isin(mito.index)
    for idx in out.index[matched]:
        gene = genes.loc[idx]
        rec = mito.loc[gene]
        out.at[idx, "MITOCARTA_STATUS"] = "YES"
        out.at[idx, "MITOCARTA_ENCODING"] = (
            "MTDNA_ENCODED_GENE" if gene in MTDNA_PROTEIN_GENES else "NUCLEAR_MITOCHONDRIAL_GENE"
        )
        for col in [
            "MITOCARTA_DESCRIPTION", "MITOCARTA_MAESTRO_SCORE", "MITOCARTA_EVIDENCE",
            "MITOCARTA_SUBCOMPARTMENT", "MITOCARTA_MITOPATHWAYS",
        ]:
            out.at[idx, col] = rec[col]

        pathways = split_pathways(rec["MITOCARTA_MITOPATHWAYS"])
        top = sorted({x.split(">", 1)[0].strip() for x in pathways if x})
        out.at[idx, "MITOCARTA_TOP_LEVEL_PATHWAYS"] = ";".join(top) if top else "."

        on_context = phenotype.loc[idx] > 0
        out.at[idx, "MITO_ON_CONTEXT"] = "YES" if on_context else "NO"
        if gene in MTDNA_PROTEIN_GENES:
            out.at[idx, "MITO_ON_ASSOCIATION_CLASS"] = (
                "MTDNA_GENE_WITH_ON_PHENOTYPE_CONTEXT"
                if on_context
                else "MTDNA_GENE_WITHOUT_ON_PHENOTYPE_CONTEXT"
            )
        else:
            out.at[idx, "MITO_ON_ASSOCIATION_CLASS"] = (
                "NUCLEAR_MITO_GENE_WITH_ON_PHENOTYPE_CONTEXT"
                if on_context
                else "NUCLEAR_MITO_GENE_WITHOUT_ON_PHENOTYPE_CONTEXT"
            )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    nuclear = out["MITOCARTA_ENCODING"].eq("NUCLEAR_MITOCHONDRIAL_GENE")
    on = out["MITO_ON_CONTEXT"].eq("YES")
    print(
        f"[OK] rows={len(out)} mitocarta_rows={int(matched.sum())} "
        f"nuclear_mito_rows={int(nuclear.sum())} nuclear_mito_on_rows={int((nuclear & on).sum())} "
        f"output={output}"
    )


if __name__ == "__main__":
    main()
