#!/usr/bin/env python3
"""Build a provenance-preserving HON gene-disease-HPO reference from Monarch KG.

The resulting table is a phenotype/gene reference, not variant-specific
pathogenicity evidence. It can therefore be used to prioritize genes affected
by SNVs, CNVs or other SVs, while variant mechanism is assessed separately.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))

VISUAL_ROOT = "HP:0000478"  # Abnormality of the eye


def load_simple_list(path):
    values = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                values.append(value)
    return values


def load_seeds(path):
    out = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["hpo_id"]] = row.get("hpo_label", "")
    return out


def split_tokens(value):
    if value is None:
        return []
    text = str(value).strip()
    if text in {"", ".", "NA", "N/A", "nan", "None"}:
        return []
    return [x.strip() for x in re.split(r"[|;,]", text) if x.strip()]


def publication_support(edge):
    values = []
    for token in split_tokens(edge.get("publications", "")):
        upper = token.upper()
        if (
            upper.startswith("PMID:")
            or upper.startswith("PMCID:")
            or upper.startswith("DOI:")
            or token.lower().startswith("http://doi.org/")
            or token.lower().startswith("https://doi.org/")
        ):
            values.append(token)
    return sorted(set(values))


def knowledge_sources(edge):
    values = []
    for key in (
        "primary_knowledge_source",
        "aggregator_knowledge_source",
        "provided_by",
    ):
        values.extend(split_tokens(edge.get(key, "")))
    return sorted(set(values))


def read_gmx(path):
    by_gene = defaultdict(set)
    if not path:
        return by_gene
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            pathway = fields[0].strip()
            for gene in fields[2:]:
                gene = gene.strip().upper()
                if gene:
                    by_gene[gene].add(pathway)
    return by_gene


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--panel", required=True)
    p.add_argument("--nodes", required=True)
    p.add_argument("--edges", required=True)
    p.add_argument("--hon-seeds", required=True)
    p.add_argument("--mitopathways-gmx", default=None)
    p.add_argument("--output", required=True)
    args = p.parse_args()

    panel = {x.upper() for x in load_simple_list(args.panel)}
    seeds = load_seeds(args.hon_seeds)
    mito = read_gmx(args.mitopathways_gmx)

    symbol_to_hgnc = {}
    node_name = {}
    hpo_name = {}
    with open(args.nodes, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ident = row.get("id", "")
            name = row.get("name") or row.get("symbol") or ""
            if ident:
                node_name[ident] = name
            if ident.startswith("HP:"):
                hpo_name[ident] = name
            if ident.startswith("HGNC:") and name:
                symbol_to_hgnc.setdefault(name.upper(), ident)
                for syn in split_tokens(row.get("synonym", "")):
                    symbol_to_hgnc.setdefault(syn.upper(), ident)

    hgnc_to_gene = {
        hgnc: gene
        for gene in panel
        if (hgnc := symbol_to_hgnc.get(gene)) is not None
    }

    gene_disease = defaultdict(list)
    direct_gene_hpo = defaultdict(list)
    parents = defaultdict(set)

    with open(args.edges, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for edge in reader:
            subject = edge.get("subject", "")
            obj = edge.get("object", "")
            category = edge.get("category", "")
            predicate = edge.get("predicate", "")

            if (
                subject.startswith("HP:")
                and obj.startswith("HP:")
                and "subclass" in predicate.lower()
            ):
                parents[subject].add(obj)

            if subject not in hgnc_to_gene:
                continue

            if "GeneToDiseaseAssociation" in category:
                gene_disease[subject].append(
                    {
                        "disease_id": obj,
                        "publications": publication_support(edge),
                        "knowledge_sources": knowledge_sources(edge),
                        "predicate": predicate,
                    }
                )
            elif (
                "GeneToPhenotypicFeatureAssociation" in category
                and obj.startswith("HP:")
            ):
                direct_gene_hpo[subject].append(
                    {
                        "hpo_id": obj,
                        "publications": publication_support(edge),
                        "knowledge_sources": knowledge_sources(edge),
                        "predicate": predicate,
                    }
                )

    selected_diseases = {
        item["disease_id"]
        for values in gene_disease.values()
        for item in values
        if item["disease_id"]
    }

    disease_hpo = defaultdict(list)
    if selected_diseases:
        with open(args.edges, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for edge in reader:
                subject = edge.get("subject", "")
                obj = edge.get("object", "")
                category = edge.get("category", "")
                if (
                    subject in selected_diseases
                    and obj.startswith("HP:")
                    and "DiseaseToPhenotypicFeatureAssociation" in category
                ):
                    disease_hpo[subject].append(
                        {
                            "hpo_id": obj,
                            "publications": publication_support(edge),
                        "knowledge_sources": knowledge_sources(edge),
                            "predicate": edge.get("predicate", ""),
                        }
                    )

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

    def visual(term):
        return VISUAL_ROOT in ancestors(term)

    disease_class = {}
    for disease, annotations in disease_hpo.items():
        terms = {x["hpo_id"] for x in annotations}
        visual_terms = {x for x in terms if visual(x)}
        nonvisual_terms = terms - visual_terms
        if visual_terms and nonvisual_terms:
            disease_class[disease] = "LIKELY_SYNDROMIC_BY_HPO_BREADTH"
        elif visual_terms:
            disease_class[disease] = "PREDOMINANTLY_OCULAR_BY_HPO_BREADTH"
        else:
            disease_class[disease] = "UNRESOLVED"

    rows = []
    seen = set()
    for hgnc, gene in sorted(hgnc_to_gene.items(), key=lambda x: x[1]):
        diseases = gene_disease.get(hgnc, [])
        if not diseases:
            diseases = [{
                "disease_id": ".",
                "publications": [],
                "knowledge_sources": [],
                "predicate": ".",
            }]

        for gd in diseases:
            disease = gd["disease_id"]
            annotations = disease_hpo.get(disease, [])
            phenotype_scope = "DISEASE_SPECIFIC"
            if not annotations:
                annotations = direct_gene_hpo.get(hgnc, [])
                phenotype_scope = (
                    "GENE_LEVEL_FALLBACK"
                    if annotations
                    else "NO_PHENOTYPE_ASSOCIATION"
                )

            if not annotations:
                annotations = [{
                    "hpo_id": ".",
                    "publications": [],
                    "knowledge_sources": [],
                    "predicate": ".",
                }]

            for hp in annotations:
                hpo_id = hp["hpo_id"]
                if hpo_id == ".":
                    hon_relevance = "NO_HPO_ASSOCIATION_AVAILABLE"
                elif hpo_id in seeds:
                    hon_relevance = "CORE_HON"
                elif visual(hpo_id):
                    hon_relevance = "VISUAL_SYSTEM_RELATED"
                elif disease != "." and any(
                    x["hpo_id"] in seeds or visual(x["hpo_id"])
                    for x in disease_hpo.get(disease, [])
                ):
                    hon_relevance = "SYSTEMIC_FEATURE_OF_HON_RELATED_DISEASE"
                else:
                    hon_relevance = "PANEL_GENE_DISEASE_CONTEXT"

                literature = sorted(
                    set(gd.get("publications", []))
                    | set(hp.get("publications", []))
                )
                sources = sorted(
                    set(gd.get("knowledge_sources", []))
                    | set(hp.get("knowledge_sources", []))
                )
                key = (gene, disease, hpo_id)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "HPO_ID": hpo_id,
                        "HPO_NAME": hpo_name.get(hpo_id, "") if hpo_id != "." else ".",
                        "HON_RELEVANCE": hon_relevance,
                        "ASSOCIATED_DISEASE_ID": disease,
                        "ASSOCIATED_DISEASE": node_name.get(disease, "") if disease != "." else ".",
                        "ASSOCIATED_GENE": gene,
                        "PANEL_GENE": "YES",
                        "ISOLATED_SYNDROMIC": disease_class.get(disease, "UNRESOLVED"),
                        "ISOLATED_SYNDROMIC_BASIS": (
                            "AUTOMATED_HPO_PHENOTYPE_BREADTH"
                            if disease != "."
                            else "NO_DISEASE_ASSOCIATION_AVAILABLE"
                        ),
                        "MITOCHONDRIAL_PATHWAY": ";".join(sorted(mito.get(gene, set()))) or ".",
                        "LITERATURE_SUPPORT": ";".join(literature) or ".",
                        "KNOWLEDGE_SOURCE": ";".join(sources) or ".",
                        "GENE_MAPPING_STATUS": "HGNC_MAPPED",
                        "GENE_DISEASE_PREDICATE": gd.get("predicate", "."),
                        "PHENOTYPE_PREDICATE": hp.get("predicate", "."),
                        "PHENOTYPE_ASSOCIATION_SCOPE": phenotype_scope,
                    }
                )

    mapped_panel_genes = set(hgnc_to_gene.values())
    for gene in sorted(panel - mapped_panel_genes):
        rows.append({
            "HPO_ID": ".",
            "HPO_NAME": ".",
            "HON_RELEVANCE": "NO_HGNC_MAPPING",
            "ASSOCIATED_DISEASE_ID": ".",
            "ASSOCIATED_DISEASE": ".",
            "ASSOCIATED_GENE": gene,
            "PANEL_GENE": "YES",
            "ISOLATED_SYNDROMIC": "UNRESOLVED",
            "ISOLATED_SYNDROMIC_BASIS": "NO_HGNC_MAPPING",
            "MITOCHONDRIAL_PATHWAY": ";".join(sorted(mito.get(gene, set()))) or ".",
            "LITERATURE_SUPPORT": ".",
            "KNOWLEDGE_SOURCE": ".",
            "GENE_MAPPING_STATUS": "NO_HGNC_MAPPING",
            "GENE_DISEASE_PREDICATE": ".",
            "PHENOTYPE_PREDICATE": ".",
            "PHENOTYPE_ASSOCIATION_SCOPE": "NO_PHENOTYPE_ASSOCIATION",
        })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "HPO_ID", "HPO_NAME", "HON_RELEVANCE", "ASSOCIATED_DISEASE_ID",
        "ASSOCIATED_DISEASE", "ASSOCIATED_GENE", "PANEL_GENE",
        "ISOLATED_SYNDROMIC", "ISOLATED_SYNDROMIC_BASIS",
        "MITOCHONDRIAL_PATHWAY", "LITERATURE_SUPPORT", "KNOWLEDGE_SOURCE",
        "GENE_MAPPING_STATUS", "GENE_DISEASE_PREDICATE", "PHENOTYPE_PREDICATE",
        "PHENOTYPE_ASSOCIATION_SCOPE",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[OK] panel_genes={len(panel)} mapped_genes={len(hgnc_to_gene)} "
        f"reference_rows={len(rows)} output={output}"
    )


if __name__ == "__main__":
    main()
