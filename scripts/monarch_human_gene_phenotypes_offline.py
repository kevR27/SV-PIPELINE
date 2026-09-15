#!/usr/bin/env python3
"""
Retrieve human gene -> phenotype associations from a local Monarch KG
(KGX TSV) dump instead of the live API.

Download the KG once:
    curl -L -o monarch-kg.tar.gz \
        https://data.monarchinitiative.org/monarch-kg/latest/monarch-kg.tar.gz
    tar -xzf monarch-kg.tar.gz
This produces monarch-kg_nodes.tsv and monarch-kg_edges.tsv.

Design mirrors the API version:
- input is a list of human gene symbols;
- resolve each symbol to an HGNC id via the nodes file (id or symbol match);
- keep only biolink:GeneToPhenotypicFeatureAssociation edges whose
  subject is HGNC:... and whose object is HP:...;
- do NOT traverse orthologs/model organisms;
- record whether an association matches one of the optic-neuropathy
  anchor HPO terms supplied in the input file.
"""

import argparse
import csv
import sys

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser()
parser.add_argument("--genes", required=True)
parser.add_argument("--hpo-anchors", required=True)
parser.add_argument("--nodes", required=True, help="monarch-kg_nodes.tsv")
parser.add_argument("--edges", required=True, help="monarch-kg_edges.tsv")
parser.add_argument("--output", required=True)
args = parser.parse_args()


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


def build_hgnc_lookup(nodes_path):
    """Map gene symbol (upper) -> HGNC id, and HGNC id -> label."""
    symbol_to_hgnc = {}
    hgnc_label = {}
    with open(nodes_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ident = row.get("id", "")
            if not ident.startswith("HGNC:"):
                continue
            label = row.get("name") or row.get("symbol") or ""
            hgnc_label[ident] = label
            if label:
                symbol_to_hgnc.setdefault(label.upper(), ident)
            # Also index any pipe-separated synonyms, if present.
            synonym_field = row.get("synonym", "")
            if synonym_field:
                for syn in synonym_field.split("|"):
                    syn = syn.strip()
                    if syn:
                        symbol_to_hgnc.setdefault(syn.upper(), ident)
    return symbol_to_hgnc, hgnc_label


def build_hpo_label_lookup(nodes_path):
    hpo_label = {}
    with open(nodes_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ident = row.get("id", "")
            if ident.startswith("HP:"):
                hpo_label[ident] = row.get("name", "")
    return hpo_label


print("Loading nodes...", file=sys.stderr)
symbol_to_hgnc, hgnc_label = build_hgnc_lookup(args.nodes)
hpo_label_lookup = build_hpo_label_lookup(args.nodes)

anchors = load_anchors(args.hpo_anchors)
genes = load_genes(args.genes)

hgnc_to_symbol = {}
unmapped = set(genes)
for symbol in genes:
    hgnc = symbol_to_hgnc.get(symbol.upper())
    if hgnc:
        hgnc_to_symbol[hgnc] = symbol
        unmapped.discard(symbol)

print(f"Resolved {len(hgnc_to_symbol)}/{len(genes)} genes to HGNC ids.", file=sys.stderr)

rows = []
found_for_hgnc = set()

print("Scanning edges...", file=sys.stderr)
with open(args.edges, encoding="utf-8") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for edge in reader:
        category = edge.get("category", "")
        if "GeneToPhenotypicFeatureAssociation" not in category:
            continue

        subject = edge.get("subject", "")
        if subject not in hgnc_to_symbol:
            continue

        obj = edge.get("object", "")
        if not obj.startswith("HP:"):
            continue

        symbol = hgnc_to_symbol[subject]
        found_for_hgnc.add(subject)
        rows.append({
            "gene_symbol": symbol,
            "hgnc_id": subject,
            "hpo_id": obj,
            "hpo_label": hpo_label_lookup.get(obj, ""),
            "association_category": category,
            "association_predicate": edge.get("predicate", ""),
            "source": edge.get("primary_knowledge_source", ""),
            "human_association": "YES",
            "optic_neuropathy_anchor": "1" if obj in anchors else "0",
        })

# Genes with no HGNC mapping at all.
for symbol in sorted(unmapped):
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

# Genes that mapped to HGNC but had zero HPO edges.
for hgnc, symbol in hgnc_to_symbol.items():
    if hgnc not in found_for_hgnc:
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

fields = [
    "gene_symbol", "hgnc_id", "hpo_id", "hpo_label",
    "association_category", "association_predicate", "source",
    "human_association", "optic_neuropathy_anchor"
]

with open(args.output, "w", newline="", encoding="utf-8") as out:
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {args.output}", file=sys.stderr)
