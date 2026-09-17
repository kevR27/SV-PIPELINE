#!/usr/bin/env python3

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate Straglr repeat loci with overlapping human genes."
    )
    p.add_argument("--straglr-tsv", required=True)
    p.add_argument("--straglr-bed", required=True)
    p.add_argument("--gene-bed", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--min-overlap", type=int, default=1)
    return p.parse_args()


def read_tsv(path):
    loci = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            sys.exit("ERROR: Straglr TSV has no header.")
        required = {
            "chrom", "start", "end", "target_repeat", "locus",
            "coverage", "genotype", "actual_repeat", "read_name",
            "copy_number", "size", "read_status",
        }
        missing = required - set(reader.fieldnames)
        if missing:
            sys.exit("ERROR: Missing Straglr TSV columns: " + ", ".join(sorted(missing)))
        for row in reader:
            key = row["locus"] or f'{row["chrom"]}:{row["start"]}-{row["end"]}'
            loci[key].append(row)
    return loci


def read_bed(path):
    records = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1]); end = int(fields[2])
            except ValueError:
                continue
            records.append({"chrom": fields[0], "start": start, "end": end, "fields": fields})
    return records


def make_locus_bed(bed_records, output):
    # Deliberately write exactly four A columns. read_gene_annotations relies on
    # this stable boundary to separate A from the variable-width gene BED (B).
    with open(output, "w", encoding="utf-8") as f:
        for i, record in enumerate(bed_records, 1):
            f.write(f'{record["chrom"]}\t{record["start"]}\t{record["end"]}\tSTRAGLR_{i}\n')


def intersect_genes(locus_bed, gene_bed, output):
    if shutil.which("bedtools") is None:
        sys.exit("ERROR: bedtools not found in PATH.")
    cmd = ["bedtools", "intersect", "-a", locus_bed, "-b", gene_bed, "-wa", "-wb"]
    with open(output, "w", encoding="utf-8") as f:
        result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        sys.exit("ERROR: bedtools intersect failed:\n" + result.stderr)


def infer_gene_symbol(b_fields):
    """Infer a gene label from BED columns after chrom/start/end.

    The repository gene BED is an annotation BED, but BED4/BED5/BED6 layouts
    are all common. Prefer a symbol-like value over Ensembl IDs, numeric score
    fields and strand fields. If no symbol-like value exists, keep the first
    non-empty annotation rather than silently dropping the overlap.
    """
    annotations = [x.strip() for x in b_fields[3:] if x.strip() not in {"", ".", "NA", "+", "-"}]
    if not annotations:
        return ""

    for value in annotations:
        if re.fullmatch(r"ENSG\d+(?:\.\d+)?", value, flags=re.IGNORECASE):
            continue
        if re.fullmatch(r"[-+]?\d+(?:\.\d+)?", value):
            continue
        # Typical HGNC symbols contain letters/numbers/hyphen/dot only.
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]*", value):
            return value

    return annotations[0]


def read_gene_annotations(path, min_overlap):
    genes = defaultdict(set)
    with open(path, encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip().split("\t")
            # A has exactly 4 fields. B must contain at least chrom/start/end/name.
            if len(fields) < 8:
                continue
            a_fields = fields[:4]
            b_fields = fields[4:]
            try:
                a_start = int(a_fields[1]); a_end = int(a_fields[2])
                b_start = int(b_fields[1]); b_end = int(b_fields[2])
            except (ValueError, IndexError):
                continue

            overlap = min(a_end, b_end) - max(a_start, b_start)
            if overlap < min_overlap:
                continue

            gene = infer_gene_symbol(b_fields)
            if gene:
                genes[a_fields[3]].add(gene)
    return genes


def unique(rows, column):
    return sorted({r[column] for r in rows if r.get(column) not in {"", ".", "NA", None}})


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
    a = parse_args()
    tsv_loci = read_tsv(a.straglr_tsv)
    bed_records = read_bed(a.straglr_bed)
    if not bed_records:
        sys.exit("ERROR: No Straglr loci found in BED.")

    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="straglr_") as tmp:
        locus_bed = os.path.join(tmp, "straglr_loci.bed")
        intersections = os.path.join(tmp, "gene_intersections.tsv")
        make_locus_bed(bed_records, locus_bed)
        intersect_genes(locus_bed, a.gene_bed, intersections)
        genes = read_gene_annotations(intersections, a.min_overlap)

        columns = [
            "chrom", "start", "end", "target_repeat", "locus", "coverage",
            "genotype", "actual_repeat", "copy_number", "size",
            "supporting_reads", "full_reads", "partial_reads", "skipped_reads",
            "failed_reads", "overlapping_genes",
        ]
        with open(a.output, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            for i, bed_record in enumerate(bed_records, 1):
                chrom = bed_record["chrom"]; start = bed_record["start"]; end = bed_record["end"]
                locus_rows = []
                for rows in tsv_loci.values():
                    for row in rows:
                        try:
                            rstart = int(row["start"]); rend = int(row["end"])
                        except ValueError:
                            continue
                        if row["chrom"] == chrom and rstart < end and rend > start:
                            locus_rows.append(row)

                if locus_rows:
                    summary = summarize(locus_rows)
                else:
                    summary = {
                        "target_repeat": "", "locus": f"{chrom}:{start}-{end}",
                        "coverage": "", "genotype": "", "actual_repeat": "",
                        "copy_number": "", "size": "", "supporting_reads": 0,
                        "full_reads": 0, "partial_reads": 0, "skipped_reads": 0,
                        "failed_reads": 0,
                    }
                summary.update(
                    chrom=chrom,
                    start=start,
                    end=end,
                    overlapping_genes=";".join(sorted(genes.get(f"STRAGLR_{i}", set()))),
                )
                writer.writerow(summary)

    print(f"Straglr annotation written to: {a.output}")


if __name__ == "__main__":
    main()
