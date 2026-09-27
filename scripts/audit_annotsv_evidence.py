#!/usr/bin/env python3
"""Audit AnnotSV source-field reporting and processed resource availability.

Reads the v3.5.10 prepared human SV BED layout, never infers a negative database
search from an empty output cell. Eight source channels and six database labels
are distinct dimensions. Cache only the expensive bundle scan, not sample output.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from pathlib import Path

from sv_evidence_common import MISSING

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
VERSION = "annotsv-evidence-audit-v1"
DATABASES = ("CLINVAR", "DBVAR", "GNOMAD", "DGV", "1000G", "CLINGEN")
TYPES = {"DEL": "loss", "DUP": "gain", "INS": "ins", "INV": "inv"}
FIELDS = tuple(f"{side}_{kind}_source" for side in ("P", "B") for kind in TYPES.values())
PATTERNS = {
    "CLINVAR": r"(?:^|[^A-Z0-9])(?:CLINVAR|CLN)(?:$|[^A-Z0-9])",
    "DBVAR": r"(?:^|[^A-Z0-9])DBVAR(?:$|[^A-Z0-9])",
    "GNOMAD": r"(?:^|[^A-Z0-9])GNOMAD(?:$|[^A-Z0-9])",
    "DGV": r"(?:^|[^A-Z0-9])DGV(?:$|[^A-Z0-9])",
    "1000G": r"(?:^|[^A-Z0-9])(?:1000G|1000 ?GENOMES)(?:$|[^A-Z0-9])",
    "CLINGEN": r"(?:^|[^A-Z0-9])(?:CLINGEN|HI3|TS3|HI40|TS40)(?:$|[^A-Z0-9])",
}


def present(value):
    return value is not None and str(value).strip().upper() not in MISSING


def mentions(value, database):
    return bool(re.search(PATTERNS[database], str(value).upper())) if present(value) else False


def resource_paths(root, build="GRCh38"):
    root = Path(root) / "Annotations_Human"
    paths = {}
    for field in FIELDS:
        side, kind, _ = field.split("_")
        directory = "FtIncludedInSV/PathogenicSV" if side == "P" else "SVincludedInFt/BenignSV"
        prefix = "pathogenic" if side == "P" else "benign"
        paths[field] = root / directory / build / f"{prefix}_{kind.title()}_SV_{build}.sorted.bed"
    return paths


def fingerprint(paths):
    result = {}
    for field, path in paths.items():
        try:
            stat = path.stat()
            result[field] = [str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]
        except FileNotFoundError:
            result[field] = [str(path.resolve()), "NOT_FOUND"]
        except OSError as error:
            result[field] = [str(path.absolute()), type(error).__name__]
    return result


def scan_resource(field, path):
    record = {"path": str(path), "status": "PROCESSED_RESOURCE_NOT_FOUND", "records": 0,
              "nonempty_source_records": 0, "malformed_records": 0,
              "database_records": {db: 0 for db in DATABASES}, "sha256": "."}
    digest = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for raw in fh:
                digest.update(raw)
                line = raw.decode("utf-8").rstrip("\r\n")
                if not line or line.startswith(("#", "track ", "browser ")):
                    continue
                cells = line.split("\t")
                expected = 7 if field.startswith("P_") else 6
                if len(cells) != expected or not cells[1].isdigit() or not cells[2].isdigit():
                    record["malformed_records"] += 1
                    continue
                record["records"] += 1
                source = cells[5] if field.startswith("P_") else cells[3]
                if present(source):
                    record["nonempty_source_records"] += 1
                    for db in DATABASES:
                        record["database_records"][db] += int(mentions(source, db))
        record["sha256"] = digest.hexdigest()
        if record["malformed_records"]:
            record["status"] = "RESOURCE_SCHEMA_UNRECOGNIZED"
        elif not record["records"]:
            record["status"] = "PROCESSED_RESOURCE_EMPTY"
        elif not record["nonempty_source_records"]:
            record["status"] = "RESOURCE_SOURCE_VALUES_EMPTY"
        else:
            record["status"] = "PROCESSED_RESOURCE_AVAILABLE"
    except FileNotFoundError:
        pass
    except (OSError, UnicodeError) as error:
        record["status"] = "RESOURCE_UNREADABLE"
        record["error"] = type(error).__name__
    return record


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False, encoding="utf-8") as fh:
        temp = fh.name
        json.dump(data, fh, indent=2)
        fh.write("\n")
    os.replace(temp, path)


def bundle_inventory(root, build, cache=None):
    paths = resource_paths(root, build)
    key = {"version": VERSION, "build": build, "files": fingerprint(paths)}
    if cache and Path(cache).is_file():
        try:
            cached = json.loads(Path(cache).read_text())
            if cached.get("key") == key:
                return cached
        except (OSError, ValueError):
            pass
    result = {"key": key, "channels": {field: scan_resource(field, path) for field, path in paths.items()}}
    if fingerprint(paths) != key["files"]:
        raise RuntimeError("AnnotSV resources changed during audit; rerun after annotation preparation finishes")
    if cache:
        atomic_json(cache, result)
    return result


def sample_columns(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""):
            digest.update(chunk)
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        channels = {field: {"column_present": field in header, "nonempty_rows": 0,
                           "database_rows": {db: 0 for db in DATABASES}} for field in FIELDS}
        count = 0
        for row in reader:
            count += 1
            for field, info in channels.items():
                value = row.get(field)
                if present(value):
                    info["nonempty_rows"] += 1
                    for db in DATABASES:
                        info["database_rows"][db] += int(mentions(value, db))
    return {"path": str(Path(path).resolve()), "sha256": digest.hexdigest(), "rows": count, "channels": channels}


def channel_status(value, field, audit, has_annotation=True):
    if present(value):
        return "OVERLAP_REPORTED"
    if not has_annotation:
        return "NO_ANNOTSV_MATCH"
    if not audit:
        return "AVAILABILITY_NOT_AUDITED"
    output = audit["sample"]["channels"][field]
    if not output["column_present"]:
        return "SOURCE_COLUMN_MISSING"
    resource = audit["bundle"]["channels"][field]
    if resource["status"] == "PROCESSED_RESOURCE_AVAILABLE" or output["nonempty_rows"]:
        return "NO_OVERLAP_REPORTED"
    return resource["status"]


def evidence_availability(ann_row, svtype, audit, has_annotation=True):
    kind = TYPES.get(svtype)
    if kind is None:
        return {"SV_PATHOGENIC_DB_STATUS": "NOT_APPLICABLE", "SV_BENIGN_DB_STATUS": "NOT_APPLICABLE",
                **{f"SV_DB_{db}_AVAILABILITY": "NOT_APPLICABLE" for db in DATABASES}}
    fields = [f"P_{kind}_source", f"B_{kind}_source"]
    result = {label: channel_status(ann_row.get(field), field, audit, has_annotation)
              for field, label in zip(fields, ["SV_PATHOGENIC_DB_STATUS", "SV_BENIGN_DB_STATUS"])}
    for db in DATABASES:
        if any(mentions(ann_row.get(field), db) for field in fields):
            status = "OVERLAP_REPORTED"
        elif not has_annotation:
            status = "NO_ANNOTSV_MATCH"
        elif not audit:
            status = "AVAILABILITY_NOT_AUDITED"
        elif any(not audit["sample"]["channels"][field]["column_present"] for field in fields):
            status = "SOURCE_COLUMN_MISSING"
        else:
            resources = [audit["bundle"]["channels"][field] for field in fields]
            observed = any(audit["sample"]["channels"][field]["database_rows"][db] for field in fields)
            represented = observed or any(r["database_records"][db] for r in resources)
            complete = all(r["status"] in {"PROCESSED_RESOURCE_AVAILABLE", "PROCESSED_RESOURCE_EMPTY"} for r in resources)
            if represented and complete:
                status = "NO_OVERLAP_REPORTED"
            elif represented:
                status = "DATABASE_OBSERVED_AVAILABILITY_INCOMPLETE"
            elif complete:
                status = "DATABASE_NOT_REPRESENTED_IN_PROCESSED_RESOURCES"
            else:
                status = "AVAILABILITY_UNDETERMINED"
        result[f"SV_DB_{db}_AVAILABILITY"] = status
    return result


def load_audit(path, annotsv):
    if not path:
        return None
    audit = json.loads(Path(path).read_text())
    if audit.get("version") != VERSION:
        raise ValueError("Unsupported AnnotSV availability audit version")
    digest = hashlib.sha256()
    with open(annotsv, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024*1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != audit["sample"]["sha256"]:
        raise ValueError("Availability audit belongs to a different/stale AnnotSV TSV; regenerate the audit")
    return audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotations-dir", required=True)
    p.add_argument("--genome-build", choices=["GRCh37", "GRCh38", "CHM13"], default="GRCh38")
    p.add_argument("--annotsv", required=True)
    p.add_argument("--cache")
    p.add_argument("--output-json", required=True)
    p.add_argument("--output-tsv", required=True)
    a = p.parse_args()
    audit = {"version": VERSION, "resource_layout": "AnnotSV-v3.5.10-prepared-human-SV-BED",
             "bundle": bundle_inventory(a.annotations_dir, a.genome_build, a.cache),
             "sample": sample_columns(a.annotsv)}
    atomic_json(a.output_json, audit)
    output = Path(a.output_tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["source_field", "database", "column_present", "nonempty_output_rows", "database_output_rows",
               "processed_resource_status", "processed_resource_records", "database_resource_records", "resource_path", "resource_sha256"]
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for field in FIELDS:
            sample, resource = audit["sample"]["channels"][field], audit["bundle"]["channels"][field]
            for db in DATABASES:
                writer.writerow(dict(zip(columns, [field, db, "YES" if sample["column_present"] else "NO",
                    sample["nonempty_rows"], sample["database_rows"][db], resource["status"], resource["records"],
                    resource["database_records"][db], resource["path"], resource["sha256"]])))
    print(f"[OK] audited {len(FIELDS)} source channels and {len(DATABASES)} database labels; no clinical negatives inferred")


if __name__ == "__main__":
    main()
