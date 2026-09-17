#!/usr/bin/env python3
"""Annotate Straglr repeat loci with overlapping human genes."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

A_COLUMNS = 4  # temporary Straglr BED: chrom, start, end, locus_id


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate Straglr repeat loci with overlapping human genes."
    )
    parser.add_argument("--straglr-tsv", required=True)
    parser.add_argument("--straglr-bed", required=True)
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-overlap", type=int, default=1)
    return parser.parse_args()


def read_tsv(path):
    loci = defaultdict(list)
    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            sys.exit("ERROR: Straglr TSV has no header.")

        required = {
            "chrom",
            "start",
            "end",
            "target_repeat",
            "locus",
            "coverage",
            "genotype",
            "actual_repeat",
            "read_name",
            "copy_number",
            "size",
            "read_status",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            sys.exit(
                "ERROR: Missing Straglr TSV columns: " + ", ".join(sorted(missing))
            )

        for row in reader:
            key = row["locus"] or f'{row["chrom"]}:{row["start"]}-{row["end"]}'
            loci[key].append(row)
    return loci


def read_bed(path):
    records = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            try:
                records.append(
                    {
                        "chrom": fields[0],
                        "start": int(fields[1]),
                        "end": int(fields[2]),
                    }
                )
            except ValueError:
                continue
    return records


def make_locus_bed(bed_records, output):
    with open(output, "w", encoding="utf-8") as fh:
        for i, record in enumerate(bed_records, 1):
            fh.write(
                f'{record["chrom"]}\t{record["start"]}\t{record["end"]}\tSTRAGLR_{i}\n'
            )


def intersect_genes(locus_bed, gene_bed, output):
    if shutil.which("bedtools") is None:
        sys.exit("ERROR: bedtools not found in PATH.")

    cmd = ["bedtools", "intersect", "-a", locus_bed, "-b", gene_bed, "-wa", "-wb"]
    with open(output, "w", encoding="utf-8") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit("ERROR: bedtools intersect failed:\n" + result.stderr)


def infer_gene_symbol(gene_fields: list[str]) -> str:
    """Return the gene symbol from common gene-BED layouts.

    The repository gene BED is expected to have chromosome/start/end followed
    by one or more annotation columns.  Prefer the first non-empty annotation
    field that resembles a gene symbol and is not a strand/score/Ensembl ID.
    """
    if len(gene_fields) < 4:
        return ""

    annotations = gene_fields[3:]
    # First annotation column is the conventional BED4 gene-name field.
    preferred = annotations[0].strip() if annotations else ""
    if preferred and preferred not in {".", "NA", "N/A"}:
        # If this field is a compound label, keep the first useful symbol-like
        # component rather than silently discarding the overlap.
        for token in re.split(r"[;,|]", preferred):
            token = token.strip()
            if token and token not in {".", "NA", "N/A", "+", "-"}:
                return token

    for value in annotations[1:]:
        value = value.strip()
        if not value or value in {".", "NA", "N/A", "+", "-"}:
            continue
        if value.startswith(("ENSG", "ENST")):
            continue
        try:
            float(value)
            continue
        except ValueError:
            pass
        return value
    return ""


def read_gene_annotations(path, min_overlap):
    genes = defaultdict(set)
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.rstrip().split("\t")
            if len(fields) < A_COLUMNS + 4:
                continue

            locus_id = fields[3]
            gene_fields = fields[A_COLUMNS:]
            try:
                a_start = int(fields[1])
                a_end = int(fields[2])
                b_start = int(gene_fields[1])
                b_end = int(gene_fields[2])
            except (ValueError, IndexError):
                continue

            overlap = min(a_end, b_end) - max(a_start, b_start)
            if overlap < min_overlap:
                continue

            gene = infer_gene_symbol(gene_fields)
            if gene:
                genes[locus_id].add(gene)
    return genes


def unique(rows, column):
    return sorted(
        {
            row[column]
            for row in rows
            if row.get(column) not in {"", ".", "NA", None}
        }
    )


def summarize(rows):
    first = rows[0]
    statuses = defaultdict(int)
    reads = set()
    for row in rows:
        if row.get("read_name"):
            reads.add(row["read_name"])
        if row.get("read_status"):
            statuses[row["read_status"]] += 1

    return {
        "target_repeat": first["target_repeat"],
        "locus": first["locus"],
        "coverage": first["coverage"],
        "genotype": ";".join(unique(rows, "genotype")),
        "actual_repeat": ";".join(unique(rows, "actual_repeat")),
        "copy_number": ";".join(unique(rows, "copy_number")),
        "size": ";".join(unique(rows, "size")),
        "supporting_reads": len(reads),
        "full_reads": statuses["full"],
        "partial_reads": statuses["partial"],
        "skipped_reads": statuses["skipped"],
        "failed_reads": statuses["failed"],
    }


def main():
    args = parse_args()
    tsv_loci = read_tsv(args.straglr_tsv)
    bed_records = read_bed(args.straglr_bed)
    if not bed_records:
        sys.exit("ERROR: No Straglr loci found in BED.")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="straglr_") as tmp:
        locus_bed = os.path.join(tmp, "straglr_loci.bed")
        intersections = os.path.join(tmp, "gene_intersections.tsv")
        make_locus_bed(bed_records, locus_bed)
        intersect_genes(locus_bed, args.gene_bed, intersections)
        genes = read_gene_annotations(intersections, args.min_overlap)

        columns = [
            "chrom",
            "start",
            "end",
            "target_repeat",
            "locus",
            "coverage",
            "genotype",
            "actual_repeat",
            "copy_number",
            "size",
            "supporting_reads",
            "full_reads",
            "partial_reads",
            "skipped_reads",
            "failed_reads",
            "overlapping_genes",
        ]

        with open(args.output, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
            writer.writeheader()

            for i, bed_record in enumerate(bed_records, 1):
                chrom = bed_record["chrom"]
                start = bed_record["start"]
                end = bed_record["end"]
                locus_rows = []

                for rows in tsv_loci.values():
                    for row in rows:
                        try:
                            row_start = int(row["start"])
                            row_end = int(row["end"])
                        except ValueError:
                            continue
                        if row["chrom"] == chrom and row_start < end and row_end > start:
                            locus_rows.append(row)

                if locus_rows:
                    summary = summarize(locus_rows)
                else:
                    summary = {
                        "target_repeat": "",
                        "locus": f"{chrom}:{start}-{end}",
                        "coverage": "",
                        "genotype": "",
                        "actual_repeat": "",
                        "copy_number": "",
                        "size": "",
                        "supporting_reads": 0,
                        "full_reads": 0,
                        "partial_reads": 0,
                        "skipped_reads": 0,
                        "failed_reads": 0,
                    }

                summary.update(
                    {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "overlapping_genes": ";".join(
                            sorted(genes.get(f"STRAGLR_{i}", set()))
                        ),
                    }
                )
                writer.writerow(summary)

    print(f"[OK] Straglr annotation written to: {args.output}")


if __name__ == "__main__":
    main()
