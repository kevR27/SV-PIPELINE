#!/usr/bin/env python3
"""Check SRS files and references before expensive analyses begin."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def record_file_check(
    path: str,
    label: str,
    checks: list[dict[str, str]],
) -> Path:
    """Record whether a required file or directory is available."""
    checked_path = Path(path)

    if checked_path.is_dir():
        available = True
    elif checked_path.is_file():
        available = checked_path.stat().st_size > 0
    else:
        available = False

    checks.append(
        {
            "check": label,
            "path": str(checked_path),
            "status": "PASS" if available else "FAIL",
        }
    )
    return checked_path


def write_results(
    checks: list[dict[str, str]],
    output_tsv: str,
    output_json: str,
) -> list[dict[str, str]]:
    """Write both a readable table and a machine-readable report."""
    failed = [check for check in checks if check["status"] == "FAIL"]

    tsv_path = Path(output_tsv)
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["check", "path", "status"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(checks)

    report = {
        "status": "FAIL" if failed else "PASS",
        "failed": failed,
        "checks": checks,
    }
    Path(output_json).write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )

    return failed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--panel-bed", required=True)
    parser.add_argument("--gene-list", required=True)
    parser.add_argument("--nuclear-mito-bed", required=True)
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--exclude-bed", required=True)
    parser.add_argument("--repeat-catalog", required=True)
    parser.add_argument("--hpo-seeds", required=True)
    parser.add_argument("--monarch-nodes", required=True)
    parser.add_argument("--monarch-edges", required=True)
    parser.add_argument("--annotsv-dir", required=True)
    parser.add_argument("--vep-cache-dir", required=True)
    parser.add_argument("--gridss-enabled", choices=["true", "false"])
    parser.add_argument("--mt-enabled", choices=["true", "false"])
    parser.add_argument("--mt-contig", default="chrM")
    parser.add_argument("--mitocarta-enabled", choices=["true", "false"])
    parser.add_argument("--mitocarta")
    parser.add_argument("--pathways-gmx")
    parser.add_argument("--bam", nargs="+", default=[])
    parser.add_argument("--bai", nargs="+", default=[])
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    args = parser.parse_args()

    checks: list[dict[str, str]] = []
    reference = record_file_check(args.reference, "reference FASTA", checks)
    reference_index = record_file_check(
        str(reference) + ".fai",
        "reference FASTA index",
        checks,
    )

    required_resources = [
        (args.panel_bed, "panel BED"),
        (args.gene_list, "candidate gene list"),
        (args.nuclear_mito_bed, "nuclear mitochondrial candidate BED"),
        (args.gene_bed, "genome-wide gene BED"),
        (args.exclude_bed, "DELLY exclusion BED"),
        (args.repeat_catalog, "ExpansionHunter catalog"),
        (args.hpo_seeds, "HPO seed file"),
        (args.monarch_nodes, "Monarch nodes"),
        (args.monarch_edges, "Monarch edges"),
        (args.annotsv_dir, "AnnotSV directory"),
        (args.vep_cache_dir, "VEP cache directory"),
    ]
    for path, label in required_resources:
        record_file_check(path, label, checks)

    for bam in args.bam:
        record_file_check(bam, "sample BAM", checks)
    for index in args.bai:
        record_file_check(index, "sample BAM index", checks)

    if args.gridss_enabled == "true":
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            record_file_check(
                str(reference) + suffix,
                "GRIDSS/BWA reference index",
                checks,
            )

    if args.mitocarta_enabled == "true":
        if args.mitocarta:
            record_file_check(args.mitocarta, "MitoCarta inventory", checks)
        else:
            checks.append(
                {
                    "check": "MitoCarta inventory",
                    "path": ".",
                    "status": "FAIL",
                }
            )

        if args.pathways_gmx:
            record_file_check(args.pathways_gmx, "MitoCarta pathways", checks)

    if args.mt_enabled == "true" and reference_index.exists():
        contigs = {
            line.split("\t", 1)[0]
            for line in reference_index.read_text(encoding="utf-8").splitlines()
            if line
        }
        checks.append(
            {
                "check": "mtDNA contig in reference",
                "path": args.mt_contig,
                "status": "PASS" if args.mt_contig in contigs else "FAIL",
            }
        )

    failed = write_results(checks, args.output_tsv, args.output_json)
    if failed:
        details = "; ".join(
            f"{check['check']}: {check['path']}" for check in failed
        )
        raise SystemExit(f"Preflight failed: {details}")

    print(f"[OK] preflight_checks={len(checks)} output={args.output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
