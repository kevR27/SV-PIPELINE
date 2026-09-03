#!/usr/bin/env python3
"""
Retrieve human gene -> phenotype associations from Monarch KG.

Design:
- input is a list of human gene symbols;
- resolve each symbol through Monarch search to an HGNC entity;
- query that HGNC entity for GeneToPhenotypicFeatureAssociation;
- do NOT traverse orthologs/model organisms;
- retain HPO phenotype associations only;
- record whether an association matches one of the optic-neuropathy
  anchor HPO terms supplied in the input file.

This is a gene-associated phenotype layer, NOT a prediction of the
patient's actual phenotype.
"""

import argparse
import csv
import time
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--genes", required=True)
parser.add_argument("--hpo-anchors", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--api-base", required=True)
parser.add_argument("--timeout", type=int, default=30)
parser.add_argument("--sleep", type=float, default=0.05)
args = parser.parse_args()

session = requests.Session()
session.headers.update({"Accept": "application/json", "User-Agent": "SV-Pipeline/1.0"})

def load_genes(path):
    genes = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            x = line.strip()
            if x and not x.startswith("#"):
                genes.append(x)
    return sorted(set(genes))

def load_anchors(path):
    anchors = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            anchors[row["hpo_id"]] = row.get("hpo_label", "")
    return anchors

def get_json(url, params=None):
    r = session.get(url, params=params, timeout=args.timeout)
    r.raise_for_status()
    return r.json()

def resolve_hgnc(symbol):
    # Monarch search is used rather than a model-organism mapping.
    data = get_json(
        f"{args.api_base}/search",
        params={"q": symbol, "limit": 20}
    )
    items = data.get("items", data.get("results", []))

    exact = []
    for item in items:
        entity = item.get("entity", item)
        ident = entity.get("id", "")
        label = entity.get("label", entity.get("name", ""))
        category = entity.get("category", "")

        if ident.startswith("HGNC:") and label.upper() == symbol.upper():
            exact.append(ident)

    if exact:
        return exact[0]

    # Fallback: first HGNC result.
    for item in items:
        entity = item.get("entity", item)
        ident = entity.get("id", "")
        if ident.startswith("HGNC:"):
            return ident

    return None

def extract_items(data):
    if isinstance(data, dict):
        for key in ("items", "results", "associations"):
            if isinstance(data.get(key), list):
                return data[key]
    return []

def field_id(obj):
    if isinstance(obj, str):
        return obj
    if isinstance(obj, dict):
        return obj.get("id", "")
    return ""

def field_label(obj):
    if isinstance(obj, str):
        return obj
    if isinstance(obj, dict):
        return obj.get("label", obj.get("name", ""))
    return ""

anchors = load_anchors(args.hpo_anchors)
genes = load_genes(args.genes)

rows = []

for symbol in genes:
    try:
        hgnc = resolve_hgnc(symbol)
        if not hgnc:
            rows.append({
                "gene_symbol": symbol,
                "hgnc_id": "",
                "hpo_id": "",
                "hpo_label": "",
                "association_category": "",
                "association_predicate": "",
                "source": "",
                "human_association": "NO_HGNC_MAPPING",
                "optic_neuropathy_anchor": "0",
            })
            continue

        time.sleep(args.sleep)

        # Direct human gene -> phenotype associations.
        url = f"{args.api_base}/entity/{hgnc}/biolink:GeneToPhenotypicFeatureAssociation"
        data = get_json(url, params={"limit": 200})
        items = extract_items(data)

        found = False
        for assoc in items:
            subject = field_id(assoc.get("subject", ""))
            obj = assoc.get("object", {})
            hpo_id = field_id(obj)
            hpo_label = field_label(obj)

            # Enforce human gene subject and HPO object.
            if subject and not subject.startswith("HGNC:"):
                continue
            if not hpo_id.startswith("HP:"):
                continue

            found = True
            rows.append({
                "gene_symbol": symbol,
                "hgnc_id": hgnc,
                "hpo_id": hpo_id,
                "hpo_label": hpo_label,
                "association_category": assoc.get("category", "GeneToPhenotypicFeatureAssociation"),
                "association_predicate": assoc.get("predicate", ""),
                "source": assoc.get("primary_knowledge_source", assoc.get("source", "")),
                "human_association": "YES",
                "optic_neuropathy_anchor": "1" if hpo_id in anchors else "0",
            })

        if not found:
            rows.append({
                "gene_symbol": symbol,
                "hgnc_id": hgnc,
                "hpo_id": "",
                "hpo_label": "",
                "association_category": "",
                "association_predicate": "",
                "source": "",
                "human_association": "NO_HPO_ASSOCIATION",
                "optic_neuropathy_anchor": "0",
            })

    except Exception as exc:
        rows.append({
            "gene_symbol": symbol,
            "hgnc_id": "",
            "hpo_id": "",
            "hpo_label": "",
            "association_category": "",
            "association_predicate": "",
            "source": f"ERROR:{type(exc).__name__}:{exc}",
            "human_association": "ERROR",
            "optic_neuropathy_anchor": "0",
        })

fields = [
    "gene_symbol", "hgnc_id", "hpo_id", "hpo_label",
    "association_category", "association_predicate", "source",
    "human_association", "optic_neuropathy_anchor"
]

with open(args.output, "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
