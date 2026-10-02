#!/usr/bin/env python3
"""Create a panel-only view from VEP tab-delimited output.

Supports both the legacy default VEP format with an Extra column and explicit
--tab output where SYMBOL and other annotations are already separate columns.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def load_genes(path: str) -> set[str]:
    with open(path, encoding="utf-8") as fh:
        return {
            line.strip().upper()
            for line in fh
            if line.strip() and not line.startswith("#")
        }


def parse_extra(extra_field: str) -> dict[str, str]:
    extra: dict[str, str] = {}
    if extra_field in {"-", "", "."}:
        return extra
    for pair in extra_field.split(";"):
        if "=" in pair:
            key, value = pair.split("=", 1)
            extra[key] = value
    return extra


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    genes = load_genes(args.genes)
    header: list[str] | None = None
    matched_rows: list[list[str]] = []

    with open(args.vep, encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fields = line.lstrip("#").rstrip("\n").split("\t")
                if "Uploaded_variation" in fields:
                    header = fields
                continue
            if header is None:
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields.extend([""] * (len(header) - len(fields)))
            row = dict(zip(header, fields))

            symbol = row.get("SYMBOL", "").strip().upper()
            if not symbol and "Extra" in row:
                symbol = parse_extra(row.get("Extra", "")).get("SYMBOL", "").upper()

            if symbol in genes:
                matched_rows.append(fields[: len(header)])

    if header is None:
        raise SystemExit("No VEP header line found in input file.")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fout:
        fout.write("\t".join(header) + "\n")
        for fields in matched_rows:
            fout.write("\t".join(fields) + "\n")

    print(f"[OK] panel_rows={len(matched_rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
