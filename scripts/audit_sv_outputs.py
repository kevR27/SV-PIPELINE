#!/usr/bin/env python3
"""Account for each master SV, including records AnnotSV did not annotate.

The output is a report. It does not remove variants or change clinical scores.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from annotsv_reconciliation import read_skip_log, record_evidence, summarize
from build_integrated_sv_gene_tsv import read_vcf, read_tsv, build_annotsv_indexes, annotsv_matches_for_sv


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--annotsv", required=True)
    p.add_argument("--unannotated")
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--output-json", required=True)
    a = p.parse_args()
    svs, _ = read_vcf(a.vcf)
    annotations = read_tsv(a.annotsv)
    indexes = build_annotsv_indexes(annotations)
    notes = read_skip_log(a.unannotated, svs)
    records = []
    for sv in svs:
        matches, method = annotsv_matches_for_sv(sv, *indexes)
        records.append({**record_evidence(sv, matches, notes, bool(a.unannotated)),
                        "CHROM": sv["CHROM"], "START": sv["START"], "SVTYPE": sv["SVTYPE"], "ANNOTSV_MATCH": method})
    Path(a.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    columns = list(records[0]) if records else ["SV_ID", "ANNOTSV_RECORD_STATUS", "ANNOTSV_MATCH"]
    with open(a.output_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, columns, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(records)
    summary = summarize(records)
    summary["input_sha256"] = {k: digest(v) for k, v in vars(a).items() if k in {"vcf", "annotsv", "unannotated"} and v}
    Path(a.output_json).parent.mkdir(parents=True, exist_ok=True)
    Path(a.output_json).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "input_sha256"}))


if __name__ == "__main__":
    main()
