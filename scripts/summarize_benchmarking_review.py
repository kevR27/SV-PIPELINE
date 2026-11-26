#!/usr/bin/env python3
"""One per-patient index of benchmarking and ranking-sensitivity results.

Truth-based benchmarking, agreement between methods, and filter sensitivity
remain three different analyses. This report does not re-rank variants.
"""
import argparse
import csv
import json
from pathlib import Path


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def report(sample, comparison, stability=None, truth=None):
    rows = []
    sources = {"rank_agreement": str(comparison)}
    agreement = read_json(comparison)
    if agreement.get("sample_id") not in (None, sample):
        raise ValueError("Rank-agreement report belongs to another sample")
    for k, values in sorted(agreement.get("top_k", {}).items(), key=lambda item: int(item[0])):
        rows.append({
            "SAMPLE": sample, "CATEGORY": "METHOD_AGREEMENT",
            "ANALYSIS": "CUSTOM_VS_GENOMISER", "SCENARIO": "top_" + k,
            "METRIC": "jaccard", "VALUE": values.get("jaccard"),
            "STATUS": "DESCRIPTIVE_ONLY", "SOURCE": str(comparison),
        })
    if stability:
        sources["ranking_stability"] = str(stability)
        with open(stability, newline="", encoding="utf-8") as handle:
            for item in csv.DictReader(handle, delimiter="\t"):
                rows.append({
                    "SAMPLE": sample, "CATEGORY": "SENSITIVITY_ANALYSIS",
                    "ANALYSIS": "FILTER_ELIGIBILITY_ONLY",
                    "SCENARIO": item.get("scenario", "."),
                    "METRIC": "baseline_top_n_retained",
                    "VALUE": item.get("baseline_top_n_retained", "."),
                    "STATUS": item.get("status", "UNKNOWN"),
                    "SOURCE": str(stability),
                })
        stability_status = "EVALUATED"
    else:
        rows.append({
            "SAMPLE": sample, "CATEGORY": "SENSITIVITY_ANALYSIS",
            "ANALYSIS": "FILTER_ELIGIBILITY_ONLY", "SCENARIO": ".",
            "METRIC": "baseline_top_n_retained", "VALUE": ".",
            "STATUS": "DISABLED", "SOURCE": ".",
        })
        stability_status = "DISABLED"

    if truth:
        sources["confirmed_truth"] = str(truth)
        validated = read_json(truth)
        if validated.get("sample_id") not in (None, sample):
            raise ValueError("Confirmed-truth benchmark belongs to another sample")
        records = validated.get("ground_truth_records", 0)
        rows.append({
            "SAMPLE": sample, "CATEGORY": "CONFIRMED_TRUTH_BENCHMARK",
            "ANALYSIS": "INDEPENDENT_CONFIRMED_VARIANTS",
            "SCENARIO": "curated_truth", "METRIC": "ground_truth_records",
            "VALUE": records, "STATUS": "EVALUATED" if records else "NO_CONFIRMED_RECORDS",
            "SOURCE": str(truth),
        })
        truth_status = "EVALUATED" if records else "NO_CONFIRMED_RECORDS"
    else:
        rows.append({
            "SAMPLE": sample, "CATEGORY": "CONFIRMED_TRUTH_BENCHMARK",
            "ANALYSIS": "INDEPENDENT_CONFIRMED_VARIANTS",
            "SCENARIO": "curated_truth", "METRIC": "ground_truth_records",
            "VALUE": ".", "STATUS": "NOT_CONFIGURED", "SOURCE": ".",
        })
        truth_status = "NOT_CONFIGURED"
    manifest = {
        "sample_id": sample,
        "ranking_stability": stability_status,
        "confirmed_truth_benchmark": truth_status,
        "rank_agreement": "DESCRIPTIVE_ONLY",
        "sources": sources,
        "interpretation": (
            "Rank agreement is not diagnostic accuracy. Filter eligibility sensitivity "
            "reselects genes from fixed baseline ranks; it is not full algorithm reranking. "
            "A truth benchmark requires independently curated confirmed variants."
        ),
    }
    return rows, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--comparison", required=True)
    parser.add_argument("--stability")
    parser.add_argument("--truth")
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()
    rows, manifest = report(args.sample, args.comparison, args.stability, args.truth)
    fields = ["SAMPLE", "CATEGORY", "ANALYSIS", "SCENARIO", "METRIC", "VALUE", "STATUS", "SOURCE"]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    m = Path(args.manifest)
    m.parent.mkdir(parents=True, exist_ok=True)
    m.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
