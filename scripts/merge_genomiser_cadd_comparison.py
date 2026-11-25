#!/usr/bin/env python3
"""Append optional Genomiser+CADD columns to the existing clinical SV table.

Do not recompute or overwrite the baseline clinical ranking.  The CADD run is
an independent Genomiser sensitivity analysis, not an isolated CADD score.
"""
import argparse
import csv
import json
from pathlib import Path

import pandas as pd
from annotate_genomiser_comparison import annotate_genomiser

PREFIX = "GENOMISER_"
CADD_PREFIX = "GENOMISER_CADD_"


def check_same_inputs(baseline_manifest, cadd_manifest):
    base = json.loads(Path(baseline_manifest).read_text(encoding="utf-8"))
    cadd = json.loads(Path(cadd_manifest).read_text(encoding="utf-8"))
    for field in ("sample_id", "small_vcf", "sv_vcf", "phenopacket_sha256"):
        if not base.get(field) or base[field] != cadd.get(field):
            raise ValueError(f"Baseline and CADD runs differ in {field}")
    if "CADD" in base.get("pathogenicity_sources", []):
        raise ValueError("Baseline Genomiser already contains CADD")
    if "CADD" not in cadd.get("pathogenicity_sources", []):
        raise ValueError("CADD Genomiser run does not include CADD")
    if base.get("analysis_yaml_sha256") == cadd.get("analysis_yaml_sha256"):
        raise ValueError("CADD and baseline analysis YAML hashes are identical")


def clean_rank(value):
    try:
        number = int(float(str(value)))
        return number if number > 0 else None
    except (TypeError, ValueError):
        return None


def merge_patient(clinical, genes, variants, baseline_manifest, cadd_manifest, output):
    check_same_inputs(baseline_manifest, cadd_manifest)
    baseline = pd.read_csv(clinical, sep="\t", dtype=str, low_memory=False).fillna(".")
    if "GENOMISER_GENE_RANK" not in baseline.columns:
        raise ValueError("Expected baseline Genomiser columns in clinical candidate TSV")
    cadd_annotation = annotate_genomiser(baseline, genes, variants)
    comparison_cols = [c for c in cadd_annotation.columns if c.startswith(PREFIX)]
    comparison = baseline.copy()
    for col in comparison_cols:
        new_col = CADD_PREFIX + col[len(PREFIX):]
        if new_col in comparison.columns:
            raise ValueError(f"CADD column already exists: {new_col}")
        comparison[new_col] = cadd_annotation[col]
    delta = []
    for base_rank, cadd_rank in zip(
        comparison["GENOMISER_GENE_RANK"], comparison["GENOMISER_CADD_GENE_RANK"]
    ):
        base_number = clean_rank(base_rank)
        cadd_number = clean_rank(cadd_rank)
        delta.append(str(base_number - cadd_number)
                     if base_number is not None and cadd_number is not None
                     else ".")
    # Positive change means gene rose toward rank 1 in the CADD run.
    comparison["GENOMISER_CADD_RANK_IMPROVEMENT"] = delta
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(path, sep="\t", index=False)


def merge_cohort(table_specs, output):
    tables = []
    seen = set()
    columns = []
    for spec in table_specs:
        if "=" not in spec:
            raise ValueError("Use --table SAMPLE=PATH")
        sample, filename = spec.split("=", 1)
        sample = sample.strip()
        if not sample or sample in seen:
            raise ValueError("Every --table needs a distinct nonempty sample ID")
        seen.add(sample)
        path = Path(filename)
        if not path.is_file():
            raise FileNotFoundError(path)
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = reader.fieldnames or []
            if not {"GENOMISER_GENE_RANK", "GENOMISER_CADD_GENE_RANK"}.issubset(fields):
                raise ValueError(f"Missing comparison columns in {path}")
            for field in fields:
                if field != "SAMPLE_ID" and field not in columns:
                    columns.append(field)
            tables.append((sample, list(reader)))
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SAMPLE_ID"] + columns,
                                delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for sample, rows in tables:
            for row in rows:
                writer.writerow({"SAMPLE_ID": sample, **{
                    k: v for k, v in row.items() if k != "SAMPLE_ID"
                }})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical")
    parser.add_argument("--cadd-genes")
    parser.add_argument("--cadd-variants")
    parser.add_argument("--baseline-manifest")
    parser.add_argument("--cadd-manifest")
    parser.add_argument("--table", action="append", default=[], help="SAMPLE=table.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.table:
        if any([args.clinical, args.cadd_genes, args.cadd_variants,
                args.baseline_manifest, args.cadd_manifest]):
            parser.error("--table cannot be mixed with per-patient inputs")
        merge_cohort(args.table, args.output)
    else:
        required = (args.clinical, args.cadd_genes, args.cadd_variants,
                    args.baseline_manifest, args.cadd_manifest)
        if not all(required):
            parser.error("Individual comparison requires clinical, both CADD TSVs and two manifests")
        merge_patient(*required, args.output)


if __name__ == "__main__":
    main()
