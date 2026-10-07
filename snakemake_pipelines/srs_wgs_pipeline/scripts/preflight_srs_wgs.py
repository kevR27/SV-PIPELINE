#!/usr/bin/env python3
"""Validate SRS-WGS inputs before any expensive caller starts."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def check(path: str, label: str, rows: list[dict], nonempty: bool = True) -> Path:
    p = Path(path)
    ok = p.exists() and (not nonempty or (p.is_file() and p.stat().st_size > 0) or p.is_dir())
    rows.append({"check": label, "path": str(p), "status": "PASS" if ok else "FAIL"})
    return p


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True)
    p.add_argument("--panel-bed", required=True)
    p.add_argument("--gene-list", required=True)
    p.add_argument("--nuclear-mito-bed", required=True)
    p.add_argument("--exclude-bed", required=True)
    p.add_argument("--repeat-catalog", required=True)
    p.add_argument("--hpo-seeds", required=True)
    p.add_argument("--monarch-nodes", required=True)
    p.add_argument("--monarch-edges", required=True)
    p.add_argument("--annotsv-dir", required=True)
    p.add_argument("--vep-cache-dir", required=True)
    p.add_argument("--gridss-enabled", choices=["true", "false"], default="false")
    p.add_argument("--mt-enabled", choices=["true", "false"], default="false")
    p.add_argument("--mt-contig", default="chrM")
    p.add_argument("--mitocarta-enabled", choices=["true", "false"], default="false")
    p.add_argument("--mitocarta")
    p.add_argument("--pathways-gmx")
    p.add_argument("--bam", nargs="+", default=[])
    p.add_argument("--bai", nargs="+", default=[])
    p.add_argument("--output-tsv", required=True)
    p.add_argument("--output-json", required=True)
    a = p.parse_args()

    rows: list[dict] = []
    ref = check(a.reference, "reference FASTA", rows)
    fai = check(str(ref) + ".fai", "reference FASTA index", rows)
    for path, label in [
        (a.panel_bed, "panel BED"), (a.gene_list, "candidate gene list"),
        (a.nuclear_mito_bed, "nuclear mitochondrial candidate BED"),
        (a.exclude_bed, "DELLY exclusion BED"), (a.repeat_catalog, "ExpansionHunter catalog"),
        (a.hpo_seeds, "HPO seed file"), (a.monarch_nodes, "Monarch nodes"),
        (a.monarch_edges, "Monarch edges"), (a.annotsv_dir, "AnnotSV directory"),
        (a.vep_cache_dir, "VEP cache directory"),
    ]:
        check(path, label, rows)
    for path in a.bam:
        check(path, "sample BAM", rows)
    for path in a.bai:
        check(path, "sample BAM index", rows)

    if a.gridss_enabled == "true":
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
            check(str(ref) + suffix, "GRIDSS/BWA reference index", rows)
    if a.mitocarta_enabled == "true":
        if not a.mitocarta:
            rows.append({"check": "MitoCarta inventory", "path": ".", "status": "FAIL"})
        else:
            check(a.mitocarta, "MitoCarta inventory", rows)
        if a.pathways_gmx:
            check(a.pathways_gmx, "MitoCarta pathways", rows)

    if a.mt_enabled == "true" and fai.exists():
        contigs = {line.split("\t", 1)[0] for line in fai.read_text(encoding="utf-8").splitlines() if line}
        rows.append({
            "check": "mtDNA contig in reference", "path": a.mt_contig,
            "status": "PASS" if a.mt_contig in contigs else "FAIL",
        })

    failed = [row for row in rows if row["status"] == "FAIL"]
    out_tsv = Path(a.output_tsv); out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "path", "status"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    Path(a.output_json).write_text(json.dumps({"status": "FAIL" if failed else "PASS", "failed": failed, "checks": rows}, indent=2) + "\n", encoding="utf-8")
    if failed:
        raise SystemExit("Preflight failed: " + "; ".join(f"{x['check']}: {x['path']}" for x in failed))
    print(f"[OK] preflight_checks={len(rows)} output={out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
