#!/usr/bin/env python3
"""Patient-to-gene HPO semantic similarity using directional Resnik coverage.

The primary ranking score asks how well each patient HPO term is explained by
the gene's reference HPO annotations. It is normalized to the patient query's
self-similarity, so extra unrelated annotations on pleiotropic genes do not
penalize the candidate. Symmetric BMA is retained only as a diagnostic field.

This is gene/disease prioritization only. It does not classify an SV and it
does not replace SV mechanism, population frequency, segregation or validation.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import pandas as pd

from ranking_common import semantic_engine


MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}


def load_patient_terms(path):
    if not path:
        return []
    terms = []
    with open(path, encoding="utf-8") as handle:
        first = handle.readline()
        handle.seek(0)
        if "\t" in first and "HPO" in first.upper():
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                value = (
                    row.get("hpo_id")
                    or row.get("HPO_ID")
                    or row.get("term")
                    or row.get("HPO")
                    or ""
                )
                value = str(value).strip()
                if value.startswith("HP:"):
                    terms.append(value)
        else:
            for line in handle:
                value = line.strip().split("\t")[0]
                if value.startswith("HP:"):
                    terms.append(value)
    return sorted(set(terms))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True, help="Gene candidate table")
    p.add_argument("--reference", required=True, help="HON gene-disease-HPO reference")
    p.add_argument(
        "--gene-phenotypes",
        required=True,
        help="Full Monarch gene-HPO table for the SV-overlapping genes",
    )
    p.add_argument("--edges", required=True, help="Monarch KG edges for HPO hierarchy")
    p.add_argument("--patient-hpo", default=None)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    genes_df = pd.read_csv(args.genes, sep="\t", dtype=str, low_memory=False)
    reference = pd.read_csv(args.reference, sep="\t", dtype=str, low_memory=False)
    gene_pheno = pd.read_csv(
        args.gene_phenotypes,
        sep="\t",
        dtype=str,
        low_memory=False,
    )
    genes = sorted(set(genes_df["GENE"].dropna().astype(str)))
    patient_terms = load_patient_terms(args.patient_hpo)

    fields = [
        "GENE", "HPO_SEMANTIC_STATUS", "PATIENT_HPO_COUNT",
        "GENE_REFERENCE_HPO_COUNT", "HPO_EXACT_MATCH_COUNT",
        "HPO_QUERY_RESNIK_MEAN", "HPO_QUERY_RESNIK_NORMALIZED",
        "HPO_BMA_RESNIK", "HPO_BMA_RESNIK_NORMALIZED",
        "HPO_SEMANTIC_METHOD", "BEST_MATCHED_PATIENT_HPO",
        "INTERPRETATION",
    ]

    if not patient_terms:
        out = pd.DataFrame([
            {
                "GENE": gene,
                "HPO_SEMANTIC_STATUS": "PATIENT_HPO_NOT_AVAILABLE",
                "PATIENT_HPO_COUNT": 0,
                "GENE_REFERENCE_HPO_COUNT": int(
                    len(
                        set(
                            gene_pheno.loc[
                                gene_pheno["gene_symbol"].eq(gene),
                                "hpo_id",
                            ]
                            .dropna()
                            .astype(str)
                        )
                        | set(
                            reference.loc[
                                reference["ASSOCIATED_GENE"].eq(gene),
                                "HPO_ID",
                            ]
                            .dropna()
                            .astype(str)
                        )
                    )
                ),
                "HPO_EXACT_MATCH_COUNT": ".",
                "HPO_QUERY_RESNIK_MEAN": ".",
                "HPO_QUERY_RESNIK_NORMALIZED": ".",
                "HPO_BMA_RESNIK": ".",
                "HPO_BMA_RESNIK_NORMALIZED": ".",
                "HPO_SEMANTIC_METHOD": (
                    "ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF"
                ),
                "BEST_MATCHED_PATIENT_HPO": ".",
                "INTERPRETATION": "No patient-specific HPO terms supplied; no semantic similarity was inferred.",
            }
            for gene in genes
        ], columns=fields)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output, sep="\t", index=False)
        print(f"[OK] genes={len(out)} patient_hpo=0 output={args.output}")
        return

    engine = semantic_engine(args.edges)
    all_gene_terms = engine["all_gene_terms"]
    pair_similarity = engine["pair_similarity"]
    bma = engine["bma"]
    normalized_bma = engine["normalized_bma"]
    directional_best_match = engine["directional_best_match"]
    normalized_query_coverage = engine["normalized_query_coverage"]

    gene_terms = defaultdict(set)

    # Primary phenotype reference: all human Monarch gene-HPO associations for
    # genes actually intersected by the sample's SVs.
    for _, row in gene_pheno.iterrows():
        gene = str(row.get("gene_symbol", "")).strip()
        hp = str(row.get("hpo_id", "")).strip()
        if gene and hp.startswith("HP:"):
            gene_terms[gene].add(hp)

    # Supplement with the richer HON panel disease-HPO reference.
    for _, row in reference.iterrows():
        gene = str(row.get("ASSOCIATED_GENE", "")).strip()
        hp = str(row.get("HPO_ID", "")).strip()
        if gene and hp.startswith("HP:"):
            gene_terms[gene].add(hp)

    rows = []
    for gene in genes:
        ref_terms = sorted(gene_terms.get(gene, set()))
        if not ref_terms:
            rows.append({
                "GENE": gene,
                "HPO_SEMANTIC_STATUS": "NO_REFERENCE_HPO_FOR_GENE",
                "PATIENT_HPO_COUNT": len(patient_terms),
                "GENE_REFERENCE_HPO_COUNT": 0,
                "HPO_EXACT_MATCH_COUNT": 0,
                "HPO_QUERY_RESNIK_MEAN": 0.0,
                "HPO_QUERY_RESNIK_NORMALIZED": 0.0,
                "HPO_BMA_RESNIK": 0.0,
                "HPO_BMA_RESNIK_NORMALIZED": 0.0,
                "HPO_SEMANTIC_METHOD": (
                    "ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF"
                ),
                "BEST_MATCHED_PATIENT_HPO": ".",
                "INTERPRETATION": "No reference HPO annotations were available for this gene.",
            })
            continue

        query_score = directional_best_match(
            patient_terms,
            ref_terms,
        )
        query_normalized = normalized_query_coverage(
            patient_terms,
            ref_terms,
        )
        bma_score = bma(patient_terms, ref_terms)
        bma_normalized = normalized_bma(patient_terms, ref_terms)

        best = []
        for patient in patient_terms:
            best_score = max(pair_similarity(patient, ref) for ref in ref_terms)
            best_refs = sorted(
                ref for ref in ref_terms
                if pair_similarity(patient, ref) == best_score
            )
            best.append(f"{patient}->{','.join(best_refs[:3])}")

        rows.append({
            "GENE": gene,
            "HPO_SEMANTIC_STATUS": "EVALUATED",
            "PATIENT_HPO_COUNT": len(patient_terms),
            "GENE_REFERENCE_HPO_COUNT": len(ref_terms),
            "HPO_EXACT_MATCH_COUNT": len(set(patient_terms) & set(ref_terms)),
            "HPO_QUERY_RESNIK_MEAN": round(query_score, 6),
            "HPO_QUERY_RESNIK_NORMALIZED": round(query_normalized, 6),
            "HPO_BMA_RESNIK": round(bma_score, 6),
            "HPO_BMA_RESNIK_NORMALIZED": round(bma_normalized, 6),
            "HPO_SEMANTIC_METHOD": (
                "ASYMMETRIC_RESNIK_QUERY_COVERAGE_NORMALIZED_TO_QUERY_SELF"
            ),
            "BEST_MATCHED_PATIENT_HPO": ";".join(best),
            "INTERPRETATION": (
                "Primary ranking uses patient-query-to-gene Resnik coverage "
                "normalized to the patient query self-match. Symmetric BMA is "
                "reported only for comparison. Semantic similarity does not "
                "establish that a specific SV is causal or mechanistically compatible."
            ),
        })

    out = pd.DataFrame(rows, columns=fields).sort_values(
        ["HPO_QUERY_RESNIK_NORMALIZED", "GENE"],
        ascending=[False, True],
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)
    print(
        f"[OK] genes={len(out)} patient_hpo={len(patient_terms)} "
        f"background_genes={len(all_gene_terms)} "
        f"evaluated={(out['HPO_SEMANTIC_STATUS'] == 'EVALUATED').sum()} output={output}"
    )


if __name__ == "__main__":
    main()
