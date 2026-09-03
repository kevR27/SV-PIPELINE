#!/usr/bin/env python3

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def args():
    p = argparse.ArgumentParser(
        description="Annotate Straglr repeat loci with overlapping genes."
    )
    p.add_argument("--straglr-tsv", required=True)
    p.add_argument("--straglr-bed", required=True)
    p.add_argument("--gene-bed", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-overlap", type=int, default=1)
    return p.parse_args()


def read_tsv(path):
    loci = defaultdict(lambda: {"rows": [], "reads": set(), "status": defaultdict(int)})

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "chrom", "start", "end", "target_repeat", "locus",
            "coverage", "genotype", "actual_repeat", "read_name",
            "copy_number", "size", "read_status"
        }

        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"Missing TSV columns: {', '.join(sorted(missing))}")

        for row in reader:
            key = (
                row["chrom"],
                row["start"],
                row["end"],
                row["locus"],
                row["target_repeat"]
            )
            loci[key]["rows"].append(row)
            if row["read_name"]:
                loci[key]["reads"].add(row["read_name"])
            if row["read_status"]:
                loci[key]["status"][row["read_status"]] += 1

    return loci


def read_bed(path):
    records = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                fields = line.rstrip().split()

            if len(fields) < 3:
                continue

            try:
                records.append((fields[0], int(fields[1]), int(fields[2])))
            except ValueError:
                continue

    return records


def create_locus_bed(loci, straglr_bed, output):
    bed_lookup = set(straglr_bed)

    with open(output, "w") as f:
        for i, key in enumerate(loci, 1):
            chrom, start, end, locus, repeat = key

            try:
                start, end = int(start), int(end)
            except ValueError:
                continue

            if (chrom, start, end) not in bed_lookup:
                continue

            f.write(
                f"{chrom}\t{start}\t{end}\t"
                f"STRAGLR_{i}\t{locus}\t{repeat}\n"
            )


def intersect_genes(locus_bed, gene_bed, output):
    if shutil.which("bedtools") is None:
        sys.exit("ERROR: bedtools was not found in PATH.")

    command = [
        "bedtools", "intersect",
        "-a", locus_bed,
        "-b", gene_bed,
        "-wa", "-wb"
    ]

    with open(output, "w") as f:
        result = subprocess.run(
            command,
            stdout=f,
            stderr=subprocess.PIPE,
            text=True
        )

    if result.returncode:
        sys.exit(f"bedtools failed:\n{result.stderr}")


def read_intersections(path, min_overlap):
    genes = defaultdict(set)

    with open(path) as f:
        for line in f:
            fields = line.rstrip().split("\t")
            if len(fields) < 10:
                continue

            a_start, a_end = int(fields[1]), int(fields[2])
            b_start, b_end = int(fields[7]), int(fields[8])
            overlap = min(a_end, b_end) - max(a_start, b_start)

            if overlap >= min_overlap and fields[9] not in {"", ".", "NA"}:
                genes[fields[3]].add(fields[9])

    return genes


def unique_values(rows, column):
    return sorted({
        row[column]
        for row in rows
        if row[column] not in {"", ".", "NA"}
    })


def summarize(rows, reads, status):
    row = rows[0]

    return {
        "chrom": row["chrom"],
        "start": row["start"],
        "end": row["end"],
        "target_repeat": row["target_repeat"],
        "locus": row["locus"],
        "coverage": row["coverage"],
        "genotype": ";".join(unique_values(rows, "genotype")),
        "actual_repeat": ";".join(unique_values(rows, "actual_repeat")),
        "copy_number": ";".join(unique_values(rows, "copy_number")),
        "size": ";".join(unique_values(rows, "size")),
        "supporting_reads": len(reads),
        "full_reads": status["full"],
        "partial_reads": status["partial"],
        "skipped_reads": status["skipped"],
        "failed_reads": status["failed"]
    }


def main():
    a = args()
    loci = read_tsv(a.straglr_tsv)
    bed = read_bed(a.straglr_bed)

    if not loci:
        sys.exit("ERROR: No loci found in Straglr TSV.")

    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="straglr_") as tmp:
        locus_bed = os.path.join(tmp, "loci.bed")
        intersections = os.path.join(tmp, "intersections.tsv")

        create_locus_bed(loci, bed, locus_bed)
        intersect_genes(
            locus_bed,
            a.gene_bed,
            intersections
        )

        genes = read_intersections(
            intersections,
            a.min_overlap
        )

        columns = [
            "chrom", "start", "end", "target_repeat", "locus",
            "coverage", "genotype", "actual_repeat", "copy_number",
            "size", "supporting_reads", "full_reads", "partial_reads",
            "skipped_reads", "failed_reads", "overlapping_genes"
        ]

        with open(a.output, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=columns,
                delimiter="\t"
            )
            writer.writeheader()

            for i, key in enumerate(loci, 1):
                data = loci[key]
                result = summarize(
                    data["rows"],
                    data["reads"],
                    data["status"]
                )
                result["overlapping_genes"] = ";".join(
                    sorted(genes.get(f"STRAGLR_{i}", set()))
                )
                writer.writerow(result)

    print(f"Output written to: {a.output}")


if __name__ == "__main__":
    main()
