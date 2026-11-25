#!/usr/bin/env python3
"""Audit archived LRS results before rerunning BAM-free downstream analysis.

The full downstream path consumes files produced during the original BAM-based
run. Checking them explicitly avoids Snakemake silently scheduling read-based
jobs or accidentally accepting partial archive folders.
"""
import argparse
import csv
from pathlib import Path


def discover_samples(root):
    if not root.is_dir():
        raise ValueError(f"Input root does not exist: {root}")
    return sorted(
        (folder.name for folder in root.iterdir()
         if folder.is_dir() and folder.name.isdecimal()),
        key=lambda value: (int(value), value),
    )


def required_sources(sample):
    prefix = f"{sample}/"
    return {
        "small_variants": prefix + f"snp_clair3/{sample}.vcf.gz",
        "merged_structural_variants": prefix + f"sv/merged/{sample}_merged_SV.vcf.gz",
        "integrated_evidence": prefix + f"gene_discovery/{sample}_integrated_SV_gene_analysis.tsv",
        "initial_gene_ranking": prefix + f"gene_discovery/{sample}_ranked_candidates.tsv",
        "gene_phenotype_evidence": prefix + f"gene_discovery/{sample}_human_gene_phenotypes.tsv",
        "straglr_annotation": prefix + f"sv/straglr/{sample}_straglr.annotated.tsv",
        "tldr_mobile_elements": prefix + f"mei/tldr/{sample}.tldr.table.txt",
        "longphase_small_variants": prefix + f"phasing_longphase/{sample}.longphase.vcf.gz",
        "longphase_sv": prefix + f"phasing_longphase/{sample}.longphase_SV.vcf.gz",
        "mosdepth_10kb_bins": prefix + f"coverage/binned/{sample}.10kb.regions.bed.gz",
    }


def optional_sources(sample):
    prefix = f"{sample}/"
    return {
        "whatshap_phased_variants": prefix + f"phasing/{sample}.phased.vcf.gz",
        "modkit_methylation": prefix + f"methylation/{sample}.cpg.bedmethyl.gz",
        "caller_support_summary": prefix + f"sv/merged/{sample}_caller_support_summary.tsv",
        "dorado_qc": prefix + f"qc/{sample}.dorado_qc.tsv",
        "read_length_qc": prefix + f"qc/{sample}.n50.tsv",
        "read_level_samplot_manifest": prefix + "plots/13_breakpoint_evidence/top_20/samplot_manifest.tsv",
    }


def audit(root, sample):
    rows = []
    for group, requirements in (
        ("REQUIRED", required_sources(sample)),
        ("OPTIONAL", optional_sources(sample)),
    ):
        for name, relative in requirements.items():
            file = root / relative
            exists = file.is_file() and file.stat().st_size > 0
            rows.append({
                "sample": sample,
                "category": group,
                "component": name,
                "status": "OK" if exists else "MISSING",
                "path": str(file),
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="Parent of numbered patient result folders")
    parser.add_argument("--sample", action="append", dest="samples",
                        help="Patient folder name, repeatable. Default: detect all numeric folders.")
    parser.add_argument("--report", help="Write a TSV inventory to this path")
    parser.add_argument("--check", action="store_true",
                        help="Exit nonzero if an archived input required by the full path is missing.")
    args = parser.parse_args()

    root = Path(args.root).expanduser().resolve()
    samples = args.samples if args.samples else discover_samples(root)
    if not samples:
        parser.error(f"No numeric sample folders found in {root}. Use --sample 11446.")

    rows = [row for sample in samples for row in audit(root, sample)]
    missing = [row for row in rows
               if row["category"] == "REQUIRED" and row["status"] == "MISSING"]

    for sample in samples:
        sample_rows = [row for row in rows if row["sample"] == sample]
        count = sum(row["category"] == "REQUIRED" and row["status"] == "OK"
                    for row in sample_rows)
        total = sum(row["category"] == "REQUIRED" for row in sample_rows)
        print(f"{sample}: {count}/{total} required archived inputs present")
        for row in sample_rows:
            if row["status"] == "MISSING":
                print(f"  {row['category']}: {row['component']}: {row['path']}")

    if args.report:
        target = Path(args.report).expanduser().resolve()
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=["sample", "category", "component", "status", "path"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved inventory: {target}")

    if args.check and missing:
        raise SystemExit(
            f"{len(missing)} required file(s) missing. Restore these archived "
            "outputs, or select only ready patient folders using --sample."
        )


if __name__ == "__main__":
    main()
