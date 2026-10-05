#!/usr/bin/env python3
"""Annotate SV-overlapping genes with MitoCarta3.0 context.

MitoCarta membership/localization and optic-neuropathy phenotype context remain
separate evidence layers. Optional Human.MitoPathways3.0.gmx is used as a
pathway fallback when the Excel inventory does not expose populated pathway
columns under the detected header.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd


MTDNA_PROTEIN_GENES = {
    "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
}
MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def first_existing(df, names):
    lower = {str(c).strip().lower(): c for c in df.columns}
    for name in names:
        if name in df.columns:
            return name
        hit = lower.get(name.lower())
        if hit is not None:
            return hit
    return None


def nonempty_count(series):
    text = series.fillna("").astype(str).str.strip().str.upper()
    return int((~text.isin(MISSING)).sum())


def select_excel_inventory(path):
    """Select the curated ~1,136-gene MitoCarta inventory sheet/header."""
    book = pd.ExcelFile(path)
    candidates = []

    for sheet_name in book.sheet_names:
        for header in range(0, 16):
            try:
                df = pd.read_excel(
                    book,
                    sheet_name=sheet_name,
                    header=header,
                    dtype=str,
                )
            except Exception:
                continue

            symbol_col = first_existing(df, ["Symbol", "Gene Symbol", "GeneSymbol"])
            if symbol_col is None:
                continue

            symbols = (
                df[symbol_col]
                .fillna("")
                .astype(str)
                .str.upper()
                .str.strip()
            )
            inventory_n = int(symbols.ne("").sum())
            pathway_col = first_existing(df, ["MitoPathways", "MitoPathway"])
            compartment_col = first_existing(
                df,
                ["Sub-compartment", "Subcompartment", "SubMitoLocalization"],
            )
            pathway_n = nonempty_count(df[pathway_col]) if pathway_col else 0
            compartment_n = nonempty_count(df[compartment_col]) if compartment_col else 0

            # Prefer: curated-size inventory, populated pathway/localization
            # columns, and then closeness to the official 1,136-gene inventory.
            curated_size = int(1000 <= inventory_n <= 1300)
            score = (
                curated_size,
                int(pathway_n > 0),
                int(compartment_n > 0),
                pathway_n + compartment_n,
                -abs(inventory_n - 1136),
            )
            candidates.append(
                (score, sheet_name, header, df, inventory_n, pathway_n, compartment_n)
            )

    if not candidates:
        raise ValueError("No MitoCarta Excel sheet with a Symbol column was found.")

    candidates.sort(key=lambda x: x[0], reverse=True)
    _, sheet_name, header, df, inventory_n, pathway_n, compartment_n = candidates[0]

    if not 1000 <= inventory_n <= 1300:
        raise ValueError(
            f"Selected MitoCarta sheet has {inventory_n} symbols; expected the "
            "curated human inventory (~1,136 genes), not an all-gene score sheet."
        )

    print(
        f"[MitoCarta] sheet={sheet_name!r} header={header} "
        f"inventory_symbols={inventory_n} pathway_nonempty={pathway_n} "
        f"subcompartment_nonempty={compartment_n}"
    )
    return df, sheet_name, header


def read_mitocarta(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    source_sheet = "."
    source_header = "."
    if path.suffix.lower() in {".xls", ".xlsx", ".xlsm"}:
        df, source_sheet, source_header = select_excel_inventory(path)
    else:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        df = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)

    symbol_col = first_existing(df, ["Symbol", "Gene Symbol", "GeneSymbol"])
    if symbol_col is None:
        raise ValueError("MitoCarta input needs a Symbol column.")

    pathway_col = first_existing(df, ["MitoPathways", "MitoPathway"])
    compartment_col = first_existing(
        df,
        ["Sub-compartment", "Subcompartment", "SubMitoLocalization"],
    )
    maestro_col = first_existing(df, ["Maestro score", "Maestro Score", "MaestroScore"])
    evidence_col = first_existing(df, ["Evidence"])
    description_col = first_existing(df, ["Description"])

    keep = pd.DataFrame(
        {
            "MITOCARTA_SYMBOL": df[symbol_col]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.strip(),
            "MITOCARTA_DESCRIPTION": (
                df[description_col].fillna(".").astype(str)
                if description_col
                else "."
            ),
            "MITOCARTA_MAESTRO_SCORE": (
                df[maestro_col].fillna(".").astype(str)
                if maestro_col
                else "."
            ),
            "MITOCARTA_EVIDENCE": (
                df[evidence_col].fillna(".").astype(str)
                if evidence_col
                else "."
            ),
            "MITOCARTA_SUBCOMPARTMENT": (
                df[compartment_col].fillna(".").astype(str)
                if compartment_col
                else "."
            ),
            "MITOCARTA_MITOPATHWAYS": (
                df[pathway_col].fillna(".").astype(str)
                if pathway_col
                else "."
            ),
        }
    )
    keep = (
        keep[keep["MITOCARTA_SYMBOL"].ne("")]
        .drop_duplicates("MITOCARTA_SYMBOL")
        .set_index("MITOCARTA_SYMBOL")
    )

    pathway_n = nonempty_count(keep["MITOCARTA_MITOPATHWAYS"])
    compartment_n = nonempty_count(keep["MITOCARTA_SUBCOMPARTMENT"])
    print(
        f"[MitoCarta] retained_genes={len(keep)} "
        f"pathway_genes={pathway_n} subcompartment_genes={compartment_n}"
    )
    return keep, str(source_sheet), str(source_header)


def read_gmx(path):
    """Return gene -> sorted pathway names from a standard GMX file."""
    if not path:
        return {}
    path = Path(path)
    if not path.exists():
        print(f"[MitoCarta] pathway GMX not found; workbook pathways only: {path}")
        return {}

    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))

    if len(rows) < 3:
        raise ValueError(f"MitoPathways GMX is unexpectedly short: {path}")

    names = [str(x).strip() for x in rows[0]]
    mapping = {}
    for col, pathway in enumerate(names):
        if not pathway:
            continue
        for row in rows[2:]:
            if col >= len(row):
                continue
            gene = str(row[col]).strip().upper()
            if gene in MISSING:
                continue
            mapping.setdefault(gene, set()).add(pathway)

    result = {gene: sorted(paths) for gene, paths in mapping.items()}
    print(
        f"[MitoCarta] GMX pathways={sum(bool(x) for x in names)} "
        f"genes_with_pathway={len(result)} source={path}"
    )
    return result


def split_pathways(value):
    text = str(value or "").strip()
    if text.upper() in MISSING:
        return []
    return [x.strip() for x in text.replace("|", ";").split(";") if x.strip()]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--integrated", required=True)
    p.add_argument("--mitocarta", required=True)
    p.add_argument("--pathways-gmx", default=None)
    p.add_argument(
        "--ranking",
        required=True,
        help="*_ranked_candidates.tsv with ON-anchor HPO evidence",
    )
    p.add_argument("--annotation-only", action="store_true", help="Write only MitoCarta annotation columns plus an internal row key.")
    p.add_argument("--output", required=True)
    p.add_argument("--summary-output", default=None)
    args = p.parse_args()

    sv = pd.read_csv(args.integrated, sep="\t", dtype=str, low_memory=False)
    gene_col = first_existing(sv, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    panel_col = first_existing(sv, ["PANEL_STATUS", "panel_gene"])
    if gene_col is None:
        raise ValueError("Integrated table needs an overlapping-gene column.")

    ranking = pd.read_csv(args.ranking, sep="\t", dtype=str, low_memory=False)
    ranking_gene_col = first_existing(ranking, ["gene", "Gene", "GENE", "SYMBOL"])
    anchor_col = first_existing(
        ranking,
        ["optic_neuropathy_anchor_HPO_count", "anchor_HPO_count"],
    )
    if ranking_gene_col is None:
        raise ValueError("Ranking table needs a gene column.")

    anchor_by_gene = {}
    if anchor_col:
        for _, row in ranking.iterrows():
            gene = str(row[ranking_gene_col]).upper().strip()
            try:
                count = int(float(row[anchor_col]))
            except Exception:
                count = 0
            anchor_by_gene[gene] = count

    mito, source_sheet, source_header = read_mitocarta(args.mitocarta)
    workbook_pathways_by_gene = {
        gene: split_pathways(value)
        for gene, value in mito["MITOCARTA_MITOPATHWAYS"].items()
    }
    gmx = read_gmx(args.pathways_gmx)

    # GMX is a fallback/augmentation. Preserve hierarchical workbook paths
    # when present; otherwise retain the curated pathway names from the GMX.
    if gmx:
        for gene, pathways in gmx.items():
            if gene not in mito.index:
                continue
            workbook_paths = split_pathways(mito.at[gene, "MITOCARTA_MITOPATHWAYS"])
            combined = sorted(set(workbook_paths) | set(pathways))
            mito.at[gene, "MITOCARTA_MITOPATHWAYS"] = (
                ";".join(combined) if combined else "."
            )

    base_columns = set(sv.columns)
    out = sv.copy()
    out["_INTEGRATED_ROW_INDEX"] = range(len(out))
    genes = out[gene_col].fillna(".").astype(str).str.upper().str.strip()
    if panel_col:
        panel_text = out[panel_col].fillna("").astype(str).str.upper().str.strip()
        panel_gene = panel_text.isin(["PANEL_GENE", "YES", "TRUE", "1"])
    else:
        panel_gene = pd.Series(False, index=out.index)

    defaults = {
        "MITOCARTA_STATUS": "NO",
        "MITOCARTA_ENCODING": "NOT_MITOCARTA",
        "MITOCARTA_DESCRIPTION": ".",
        "MITOCARTA_MAESTRO_SCORE": ".",
        "MITOCARTA_EVIDENCE": ".",
        "MITOCARTA_SUBCOMPARTMENT": ".",
        "MITOCARTA_MITOPATHWAYS": ".",
        "MITOCARTA_TOP_LEVEL_PATHWAYS": ".",
        "MITOCARTA_SOURCE_SHEET": source_sheet,
        "MITOCARTA_SOURCE_HEADER": source_header,
        "MITOCARTA_PATHWAY_SOURCE": "NONE",
        "MITO_ON_ANCHOR_HPO_COUNT": 0,
        "MITO_ON_CONTEXT": "NO",
        "MITO_ON_CONTEXT_SOURCE": "NONE",
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
            "MTDNA_ENCODED_GENE"
            if gene in MTDNA_PROTEIN_GENES
            else "NUCLEAR_MITOCHONDRIAL_GENE"
        )
        for col in [
            "MITOCARTA_DESCRIPTION",
            "MITOCARTA_MAESTRO_SCORE",
            "MITOCARTA_EVIDENCE",
            "MITOCARTA_SUBCOMPARTMENT",
            "MITOCARTA_MITOPATHWAYS",
        ]:
            out.at[idx, col] = rec[col]

        pathways = split_pathways(rec["MITOCARTA_MITOPATHWAYS"])
        top = sorted(
            {
                x.split(">", 1)[0].strip()
                for x in pathways
                if ">" in x and x.split(">", 1)[0].strip()
            }
        )
        out.at[idx, "MITOCARTA_TOP_LEVEL_PATHWAYS"] = (
            ";".join(top) if top else "."
        )

        workbook_original = workbook_pathways_by_gene.get(gene, [])
        if gene in gmx and workbook_original:
            out.at[idx, "MITOCARTA_PATHWAY_SOURCE"] = "WORKBOOK_AND_GMX"
        elif gene in gmx:
            out.at[idx, "MITOCARTA_PATHWAY_SOURCE"] = "GMX"
        elif workbook_original:
            out.at[idx, "MITOCARTA_PATHWAY_SOURCE"] = "WORKBOOK"

        anchor_count = int(anchor_by_gene.get(gene, 0))
        is_panel = bool(panel_gene.loc[idx])
        on_context = is_panel or anchor_count > 0

        out.at[idx, "MITO_ON_ANCHOR_HPO_COUNT"] = anchor_count
        out.at[idx, "MITO_ON_CONTEXT"] = "YES" if on_context else "NO"

        if is_panel and anchor_count > 0:
            context_source = "ON_PANEL_AND_ANCHOR_HPO"
        elif is_panel:
            context_source = "ON_PANEL_GENE"
        elif anchor_count > 0:
            context_source = "ON_ANCHOR_HPO"
        else:
            context_source = "NONE"
        out.at[idx, "MITO_ON_CONTEXT_SOURCE"] = context_source

        encoding = "MTDNA" if gene in MTDNA_PROTEIN_GENES else "NUCLEAR_MITO"
        out.at[idx, "MITO_ON_ASSOCIATION_CLASS"] = (
            f"{encoding}_{context_source}"
            if on_context
            else f"{encoding}_NO_ON_SPECIFIC_CONTEXT"
        )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if args.annotation_only:
        annotation_columns = [
            c for c in out.columns
            if c not in base_columns and c != "_INTEGRATED_ROW_INDEX"
        ]
        output_df = out[["_INTEGRATED_ROW_INDEX"] + annotation_columns]
    else:
        output_df = out.drop(columns=["_INTEGRATED_ROW_INDEX"], errors="ignore")
    output_df.to_csv(output, sep="\t", index=False)

    nuclear = out["MITOCARTA_ENCODING"].eq("NUCLEAR_MITOCHONDRIAL_GENE")
    on = out["MITO_ON_CONTEXT"].eq("YES")
    path_reported = ~out["MITOCARTA_MITOPATHWAYS"].fillna(".").isin([".", ""])
    summary = pd.DataFrame(
        [
            {
                "input_rows": len(out),
                "unique_input_genes": int(genes[~genes.isin(MISSING)].nunique()),
                "mitocarta_rows": int(matched.sum()),
                "unique_mitocarta_genes": int(
                    out.loc[out["MITOCARTA_STATUS"].eq("YES"), gene_col].nunique()
                ),
                "nuclear_mito_rows": int(nuclear.sum()),
                "unique_nuclear_mito_genes": int(
                    out.loc[nuclear, gene_col].nunique()
                ),
                "nuclear_mito_on_rows": int((nuclear & on).sum()),
                "rows_with_pathways": int(path_reported.sum()),
                "rows_with_subcompartment": int(
                    (
                        ~out["MITOCARTA_SUBCOMPARTMENT"]
                        .fillna(".")
                        .isin([".", ""])
                    ).sum()
                ),
                "workbook_sheet": source_sheet,
                "workbook_header": source_header,
                "gmx_gene_count": len(gmx),
                "pathway_source_counts": ";".join(
                    f"{key}={value}"
                    for key, value in out["MITOCARTA_PATHWAY_SOURCE"]
                    .value_counts()
                    .to_dict()
                    .items()
                ),
            }
        ]
    )

    summary_output = (
        Path(args.summary_output)
        if args.summary_output
        else output.with_name(output.stem + "_summary.tsv")
    )
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_output, sep="\t", index=False)

    print(
        f"[OK] rows={len(out)} mitocarta_rows={int(matched.sum())} "
        f"nuclear_mito_rows={int(nuclear.sum())} "
        f"nuclear_mito_on_rows={int((nuclear & on).sum())} "
        f"rows_with_pathways={int(path_reported.sum())} "
        f"summary={summary_output} output={output}"
    )


if __name__ == "__main__":
    main()
