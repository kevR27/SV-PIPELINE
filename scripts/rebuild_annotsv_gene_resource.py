#!/usr/bin/env python3
"""Rebuild AnnotSV's human GRCh38 Ensembl genes from a local GTF.

UCSC gtfToGenePred supplies the transcript geometry. Its tab-separated info
file supplies gene names and gene IDs, including IDs for unnamed genes.
The existing reference and the installed annotation directory are never edited.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path


CANONICAL = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}
BAD_NAMES = {"CMPL", "INCMPL"}
INFO_FIELDS = ["transId", "geneId", "source", "chrom", "start", "end",
               "strand", "proteinId", "geneName", "transcriptName",
               "geneType", "transcriptType"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_id(value, prefix):
    match = re.fullmatch(rf"({prefix}\d+)(?:\.(\d+))?", value)
    if not match:
        raise ValueError(f"Expected a human Ensembl {prefix} identifier: {value!r}")
    return match.group(1), match.group(2) or ""


def chromosome(value):
    value = value.removeprefix("chr")
    return "M" if value == "MT" else value


def validate_bed(fields, location):
    if len(fields) != 10:
        raise ValueError(f"{location}: expected ten tab-separated BED fields")
    chrom, start, end, strand, name, tx, cds_start, cds_end, left, right = fields
    start, end, cds_start, cds_end = map(int, (start, end, cds_start, cds_end))
    starts = [int(x) for x in left.rstrip(",").split(",")]
    ends = [int(x) for x in right.rstrip(",").split(",")]
    if chrom not in CANONICAL or strand not in {"+", "-"}:
        raise ValueError(f"{location}: unexpected chromosome or strand")
    stable_id(tx, "ENST")
    if not name or any(c.isspace() for c in name):
        raise ValueError(f"{location}: empty gene identifier or whitespace in it")
    if not 0 <= start < end or not start <= cds_start <= cds_end <= end:
        raise ValueError(f"{location}: invalid transcript/CDS bounds")
    if len(starts) != len(ends) or starts[0] != start or ends[-1] != end:
        raise ValueError(f"{location}: exon count or terminal exon mismatch")
    if any(not start <= a < b <= end for a, b in zip(starts, ends)):
        raise ValueError(f"{location}: an exon lies outside the transcript")
    if any(b > a for b, a in zip(ends, starts[1:])):
        raise ValueError(f"{location}: exons overlap or are out of order")


def read_bed(path):
    rows = {}
    with Path(path).open() as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            validate_bed(fields, f"{path}:{n}")
            if fields[5] in rows:
                raise ValueError(f"{path}:{n}: duplicate transcript {fields[5]}")
            rows[fields[5]] = fields
    if not rows:
        raise ValueError(f"No transcripts in {path}")
    return rows


def read_info(path):
    rows = {}
    with Path(path).open() as handle:
        header = handle.readline().rstrip("\r\n").lstrip("#").split("\t")
        if header != INFO_FIELDS:
            raise ValueError("Unexpected gtfToGenePred -infoOut header")
        for n, fields in enumerate(csv.reader(handle, delimiter="\t"), 2):
            if len(fields) != len(INFO_FIELDS):
                raise ValueError(f"{path}:{n}: malformed transcript information")
            row = dict(zip(INFO_FIELDS, fields))
            if row["transId"] in rows:
                raise ValueError(f"Duplicate transcript information: {row['transId']}")
            rows[row["transId"]] = row
    return rows


def convert_tables(gene_pred, info_path):
    info = read_info(info_path)
    rows, identities, skipped, seen = {}, {}, Counter(), set()
    with Path(gene_pred).open() as handle:
        for n, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 15:
                raise ValueError(f"{gene_pred}:{n}: expected 15 genePredExt fields")
            full_tx, chrom, strand, start, end, cs, ce, count, starts, ends = fields[:10]
            if full_tx in seen:
                raise ValueError(f"Duplicate genePred transcript: {full_tx}")
            seen.add(full_tx)
            source = info.get(full_tx)
            if source is None:
                raise ValueError(f"Missing source information for {full_tx}")
            if source["geneId"] != fields[11]:
                raise ValueError(f"Conflicting source gene IDs for {full_tx}")
            if chromosome(source["chrom"]) != chromosome(chrom) or source["strand"] != strand:
                raise ValueError(f"Conflicting chromosome or strand for {full_tx}")
            chrom = chromosome(chrom)
            if chrom not in CANONICAL:
                skipped[chrom] += 1
                continue
            tx, version = stable_id(full_tx, "ENST")
            gene_id, gene_version = stable_id(source["geneId"], "ENSG")
            symbol = source["geneName"]
            if symbol.upper() in BAD_NAMES:
                raise ValueError(f"The source GTF itself has an invalid gene name: {full_tx}")
            if symbol == ".":
                symbol = ""
            name = symbol or gene_id
            bed = [chrom, start, end, strand, name, tx, cs, ce, starts, ends]
            validate_bed(bed, full_tx)
            if int(count) != len(starts.rstrip(",").split(",")):
                raise ValueError(f"Incorrect exon count for {full_tx}")
            if tx in rows:
                raise ValueError(f"Multiple versions of transcript {tx}")
            rows[tx] = bed
            identities[tx] = [tx, version, gene_id, gene_version, symbol, name,
                              "GENE_SYMBOL" if symbol else "ENSEMBL_GENE_ID_NO_SYMBOL",
                              source["geneType"], source["transcriptType"]]
    if not rows:
        raise ValueError("No canonical human transcripts were produced")
    if seen != set(info):
        raise ValueError("genePred and transcript information contain different transcript sets")
    return rows, identities, dict(skipped)


def geometry(row):
    return tuple(row[:4] + row[6:8]) + tuple(
        tuple(int(v) for v in field.rstrip(",").split(",")) for field in row[8:10])


def compare(old, new):
    counts, details = Counter(), []
    for tx in sorted(set(old) | set(new)):
        before, after = old.get(tx), new.get(tx)
        statuses = []
        if before is None:
            statuses.append("ADDED_TRANSCRIPT")
        elif after is None:
            statuses.append("REMOVED_TRANSCRIPT")
        else:
            if geometry(before) != geometry(after):
                statuses.append("TRANSCRIPT_GEOMETRY_CHANGED")
            if before[4] != after[4]:
                statuses.append("INVALID_NAME_REPLACED" if before[4].upper() in BAD_NAMES
                                else "EXISTING_GENE_NAME_CHANGED")
        if statuses:
            counts.update(statuses)
            details.append([tx, ";".join(statuses), before[4] if before else ".",
                            after[4] if after else "."])
        else:
            counts["UNCHANGED_TRANSCRIPT"] += 1
    return dict(counts), details


def write_table(path, header, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def copy_gtf(source, destination):
    opener = gzip.open if source.suffix == ".gz" else open
    headers = []
    with opener(source, "rt") as handle, destination.open("w") as output:
        for line in handle:
            if line.startswith("#!"):
                headers.append(line.rstrip("\r\n"))
            output.write(line)
    if not any(re.match(r"#!genome-build\s+GRCh38(?:\.p\d+)?\s*$", h) for h in headers):
        raise ValueError("Expected an Ensembl GTF with a '#!genome-build GRCh38' header")
    return headers


def rebuild(args):
    gtf, existing, out = args.gtf.resolve(), args.existing_bed.resolve(), args.outdir.resolve()
    if out.exists():
        raise ValueError(f"Choose a new output folder; {out} already exists")
    if not gtf.is_file() or not existing.is_file():
        raise ValueError("The GTF and existing BED must both be readable local files")
    if args.source_release < 1:
        raise ValueError("Use the positive Ensembl release number for --source-release")
    named_release = re.search(r"Homo_sapiens\.GRCh38\.(\d+)(?:\.|$)", gtf.name)
    if named_release and int(named_release.group(1)) != args.source_release:
        raise ValueError("--source-release does not match the release in the GTF filename")
    executable = shutil.which(args.gtf_to_genepred)
    if executable is None:
        raise ValueError("gtfToGenePred was not found; use --gtf-to-genepred /path/to/gtfToGenePred")
    executable = str(Path(executable).resolve())
    print("[1/4] Checking the existing reference", flush=True)
    old = read_bed(existing)
    out.mkdir(parents=True)
    manifest = {"status": "BUILDING", "genome_build": "GRCh38", "source_release": args.source_release,
                "gtf": str(gtf), "gtf_sha256": sha256(gtf),
                "existing_bed": str(existing), "existing_bed_sha256": sha256(existing),
                "script_sha256": sha256(__file__), "gtfToGenePred": executable,
                "gtfToGenePred_sha256": sha256(executable),
                "allow_reference_update": args.allow_reference_update}
    try:
        with tempfile.TemporaryDirectory(prefix="conversion_", dir=out) as tmp:
            tmp = Path(tmp)
            source = tmp / "source.gtf"
            print("[2/4] Reading the GTF and converting its transcripts", flush=True)
            manifest["gtf_headers"] = copy_gtf(gtf, source)
            gp, info = tmp / "genes.genePred", tmp / "transcript_info.tsv"
            command = [executable, "-genePredExt", "-includeVersion", f"-infoOut={info}",
                       str(source), str(gp)]
            manifest["conversion_options"] = ["-genePredExt", "-includeVersion", "-infoOut"]
            with (out / "gtfToGenePred.log").open("w") as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            new, identities, skipped = convert_tables(gp, info)
        print("[3/4] Comparing transcript identities and coordinates", flush=True)
        counts, differences = compare(old, new)
        manifest.update({"old_transcripts": len(old), "new_transcripts": len(new),
                         "old_invalid_names": sum(r[4].upper() in BAD_NAMES for r in old.values()),
                         "names_from_gene_ids": sum(not r[4] for r in identities.values()),
                         "skipped_noncanonical_transcripts": skipped, "comparison": counts})
        write_table(out / "reference_changes.tsv", ["TRANSCRIPT_ID", "CHANGE", "OLD_GENE", "NEW_GENE"], differences)
        blocked = {k: counts[k] for k in ("REMOVED_TRANSCRIPT", "TRANSCRIPT_GEOMETRY_CHANGED",
                                         "EXISTING_GENE_NAME_CHANGED") if counts.get(k)}
        if blocked and not args.allow_reference_update:
            raise ValueError(f"The source differs beyond repairing missing names: {blocked}. "
                             "Read reference_changes.tsv and use the matching GTF. "
                             "--allow-reference-update is only for an intentional annotation release change.")
        ordered = sorted(new, key=lambda tx: (new[tx][0], int(new[tx][1]), int(new[tx][2]), tx))
        genes_dir = out / "Annotations_Human" / "Genes" / "GRCh38"
        genes_dir.mkdir(parents=True)
        bed_path = genes_dir / "genes.ENSEMBL.sorted.bed"
        version_path = genes_dir / "transcript_version.ENSEMBL.tsv"
        write_table(bed_path, None, (new[tx] for tx in ordered))
        write_table(version_path, None, ([tx, identities[tx][1]] for tx in ordered))
        write_table(out / "gene_identities.tsv",
                    ["TRANSCRIPT_ID", "TRANSCRIPT_VERSION", "GENE_ID", "GENE_VERSION", "GENE_SYMBOL",
                     "ANNOTSV_GENE_LABEL", "LABEL_SOURCE", "GENE_TYPE", "TRANSCRIPT_TYPE"],
                    (identities[tx] for tx in ordered))
        manifest.update({"status": "READY", "bed_sha256": sha256(bed_path),
                         "transcript_versions_sha256": sha256(version_path)})
        print(f"[4/4] Created {len(new):,} transcript records in {genes_dir}", flush=True)
        print(f"Gene IDs retained where no symbol was supplied: {manifest['names_from_gene_ids']:,}")
        print("The installed reference has not been changed. See the installation guide before use.")
    except Exception as error:
        manifest.update({"status": "FAILED", "error": str(error)})
        raise
    finally:
        (out / "reference_build.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True, help="Local Ensembl GRCh38 .gtf or .gtf.gz")
    parser.add_argument("--existing-bed", type=Path, required=True)
    parser.add_argument("--source-release", type=int, required=True, help="Recorded Ensembl release number")
    parser.add_argument("--outdir", type=Path, required=True, help="A new, separate staging folder")
    parser.add_argument("--gtf-to-genepred", default="gtfToGenePred")
    parser.add_argument("--allow-reference-update", action="store_true")
    args = parser.parse_args()
    try:
        rebuild(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"[ERROR] {error}\n")


if __name__ == "__main__":
    main()
