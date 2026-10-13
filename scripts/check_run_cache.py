#!/usr/bin/env python3
"""Reuse costly tool outputs only when inputs/settings and saved outputs match.

Large BAM/reference/backend files are identified by path, size and nanosecond
mtime. Content inputs (e.g. the query VCF) and generated outputs are hashed.
This is run provenance, not validation of biological accuracy.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def checksum(path, uncompressed=False):
    digest = hashlib.sha256()
    opener = gzip.open if uncompressed and str(path).endswith(".gz") else open
    with opener(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def snapshot(inputs, content_inputs, tags):
    records = []
    for source in inputs:
        source = Path(source).resolve()
        files = sorted(p for p in source.rglob("*") if p.is_file()) if source.is_dir() else [source]
        if not files:
            raise ValueError(f"Empty cache input directory: {source}")
        for path in files:
            stat = path.stat()
            if not stat.st_size:
                continue  # Empty lock/marker files are not scientific resources.
            row = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
            if stat.st_size <= 1024 * 1024:
                row["sha256"] = checksum(path)
            records.append(row)
    for source in content_inputs:
        records.append({"path": str(Path(source).resolve()), "uncompressed_sha256": checksum(source, uncompressed=True)})
    return {"inputs": records, "settings": tags}


def output_records(paths):
    result = []
    for path in paths:
        path = Path(path).resolve()
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f"Missing/empty cached output: {path}")
        # Reading gzip through its trailer also detects incomplete compressed outputs.
        result.append({"path": str(path), "sha256": checksum(path, uncompressed=True)})
    return result


def matches(manifest, current):
    try:
        saved = json.loads(Path(manifest).read_text())
        return (saved["snapshot"] == current and bool(saved["outputs"])
                and saved["outputs"] == output_records([r["path"] for r in saved["outputs"]]))
    except (OSError, ValueError, KeyError, TypeError, EOFError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["snapshot", "check", "complete", "output"])
    parser.add_argument("--inputs", nargs="*", default=[])
    parser.add_argument("--content-inputs", nargs="*", default=[])
    parser.add_argument("--tags-json", default="{}")
    parser.add_argument("--tool-paths", nargs="*", default=[], help="Paths passed to a tool that does not quote its own shell variables")
    parser.add_argument("--snapshot")
    parser.add_argument("--manifest")
    parser.add_argument("--outputs", nargs="+", default=[])
    parser.add_argument("--index", type=int, default=0)
    args = parser.parse_args()
    if any(any(char.isspace() or char in "*?[]" for char in path) for path in args.tool_paths):
        parser.error("needLR 4.1 uses unquoted paths internally. Use tool input/output paths without whitespace or glob characters.")
    if args.command == "snapshot":
        data = snapshot(args.inputs, args.content_inputs, json.loads(args.tags_json))
        Path(args.snapshot).write_text(json.dumps(data, sort_keys=True) + "\n")
    elif args.command == "check":
        return 0 if matches(args.manifest, json.loads(Path(args.snapshot).read_text())) else 1
    elif args.command == "complete":
        data = {"snapshot": json.loads(Path(args.snapshot).read_text()), "outputs": output_records(args.outputs)}
        target = Path(args.manifest)
        temporary = target.with_suffix(target.suffix + ".tmp")
        temporary.write_text(json.dumps(data, sort_keys=True) + "\n")
        temporary.replace(target)
    else:
        print(json.loads(Path(args.manifest).read_text())["outputs"][args.index]["path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
