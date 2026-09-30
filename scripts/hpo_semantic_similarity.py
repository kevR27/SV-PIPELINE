#!/usr/bin/env python3
"""Patient-to-gene HPO semantic similarity using Resnik best-match average.

This is gene/disease prioritization only. It does not classify an SV and it
does not replace SV mechanism, population frequency, segregation or validation.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import pandas as pd


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


def load_parent_graph(edges_path):
    parents = defaultdict(set)
    with open(edges_path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for edge in reader:
            subject = edge.get("subject", "")
            obj = edge.get("object", "")
            predicate = edge.get("predicate", "")
            if (
                subject.startswith("HP:")
                and obj.startswith("HP:")
                and "subclass" in predicate.lower()
            ):
                parents[subject].add(obj)
    return parents


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genes", required=True, help="Gene candidate table")
    p.add_argument("--reference", required=True, help="HON gene-disease-HPO reference")
    p.add_argument("--edges", required=True, help="Monarch KG edges for HPO hierarchy")
    p.add_argument("--patient-hpo", default=None)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    genes_df = pd.read_csv(args.genes, sep="\t", dtype=str, low_memory=False)
    reference = pd.read_csv(args.reference, sep="\t", dtype=str, low_memory=False)
    genes = sorted(set(genes_df["GENE"].dropna().astype(str)))
    patient_terms = load_patient_terms(args.patient_hpo)

    fields = [
        "GENE", "HPO_SEMANTIC_STATUS", "PATIENT_HPO_COUNT",
        "GENE_REFERENCE_HPO_COUNT", "HPO_EXACT_MATCH_COUNT",
        "HPO_BMA_RESNIK", "HPO_BMA_RESNIK_NORMALIZED",
        "BEST_MATCHED_PATIENT_HPO",
        "INTERPRETATION",
    ]

    if not patient_terms:
        out = pd.DataFrame([
            {
                "GENE": gene,
                "HPO_SEMANTIC_STATUS": "PATIENT_HPO_NOT_AVAILABLE",
                "PATIENT_HPO_COUNT": 0,
                "GENE_REFERENCE_HPO_COUNT": int(
                    reference.loc[reference["ASSOCIATED_GENE"].eq(gene), "HPO_ID"].nunique()
                ),
                "HPO_EXACT_MATCH_COUNT": ".",
                "HPO_BMA_RESNIK": ".",
                "HPO_BMA_RESNIK_NORMALIZED": ".",
                "BEST_MATCHED_PATIENT_HPO": ".",
                "INTERPRETATION": "No patient-specific HPO terms supplied; no semantic similarity was inferred.",
            }
            for gene in genes
        ], columns=fields)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output, sep="\t", index=False)
        print(f"[OK] genes={len(out)} patient_hpo=0 output={args.output}")
        return

    parents = load_parent_graph(args.edges)

    @lru_cache(maxsize=None)
    def ancestors(term):
        result = {term}
        stack = list(parents.get(term, ()))
        while stack:
            current = stack.pop()
            if current in result:
                continue
            result.add(current)
            stack.extend(parents.get(current, ()))
        return frozenset(result)

    gene_terms = defaultdict(set)
    for _, row in reference.iterrows():
        gene = str(row.get("ASSOCIATED_GENE", "")).strip()
        hp = str(row.get("HPO_ID", "")).strip()
        if gene and hp.startswith("HP:"):
            gene_terms[gene].add(hp)

    # Information content based on the proportion of reference genes annotated
    # to a term or one of its descendants after ancestor propagation.
    reference_genes = sorted(gene_terms)
    propagated_count = defaultdict(int)
    for gene in reference_genes:
        propagated = set()
        for term in gene_terms[gene]:
            propagated.update(ancestors(term))
        for term in propagated:
            propagated_count[term] += 1

    n_genes = max(len(reference_genes), 1)

    def ic(term):
        # add-one smoothing avoids infinite values for patient-only terms
        return -math.log((propagated_count.get(term, 0) + 1) / (n_genes + 1))

    max_ic = max([ic(term) for term in propagated_count] or [1.0])

    @lru_cache(maxsize=None)
    def pair_similarity(a, b):
        common = ancestors(a) & ancestors(b)
        if not common:
            return 0.0
        return max(ic(term) for term in common)

    def bma(left, right):
        if not left or not right:
            return 0.0
        left_best = [max(pair_similarity(a, b) for b in right) for a in left]
        right_best = [max(pair_similarity(b, a) for a in left) for b in right]
        return (sum(left_best) / len(left_best) + sum(right_best) / len(right_best)) / 2

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
                "HPO_BMA_RESNIK": 0.0,
                "HPO_BMA_RESNIK_NORMALIZED": 0.0,
                "BEST_MATCHED_PATIENT_HPO": ".",
                "INTERPRETATION": "No reference HPO annotations were available for this gene.",
            })
            continue

        score = bma(patient_terms, ref_terms)
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
            "HPO_BMA_RESNIK": round(score, 6),
            "HPO_BMA_RESNIK_NORMALIZED": round(score / max_ic if max_ic else 0.0, 6),
            "BEST_MATCHED_PATIENT_HPO": ";".join(best),
            "INTERPRETATION": (
                "Phenotype-to-gene semantic similarity only; it does not establish "
                "that a specific SV is causal or mechanistically compatible."
            ),
        })

    out = pd.DataFrame(rows, columns=fields).sort_values(
        ["HPO_BMA_RESNIK_NORMALIZED", "GENE"],
        ascending=[False, True],
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)
    print(
        f"[OK] genes={len(out)} patient_hpo={len(patient_terms)} "
        f"evaluated={(out['HPO_SEMANTIC_STATUS'] == 'EVALUATED').sum()} output={output}"
    )


if __name__ == "__main__":
    main()
