#!/usr/bin/env python3
"""Create a panel-only view from VEP tab-delimited output."""

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
    base_fields = None
    extra_idx = None
    extra_keys: list[str] = []
    extra_keys_seen: set[str] = set()
    matched_rows: list[tuple[list[str], dict[str, str]]] = []

    with open(args.vep, encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                base_fields = line.lstrip("#").rstrip("\n").split("\t")
                extra_idx = base_fields.index("Extra") if "Extra" in base_fields else None
                continue
            if base_fields is None or extra_idx is None:
                continue

            fields = line.rstrip("\n").split("\t")
            if extra_idx >= len(fields):
                continue
            extra_dict = parse_extra(fields[extra_idx])
            symbol = extra_dict.get("SYMBOL", "").upper()
            if not symbol or symbol not in genes:
                continue

            for key in extra_dict:
                if key not in extra_keys_seen:
                    extra_keys_seen.add(key)
                    extra_keys.append(key)
            matched_rows.append((fields, extra_dict))

    if base_fields is None:
        raise SystemExit("No VEP header line found in input file.")
    if extra_idx is None:
        raise SystemExit("VEP output has no Extra column; cannot filter by SYMBOL.")

    out_base_fields = [field for i, field in enumerate(base_fields) if i != extra_idx]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as fout:
        fout.write("\t".join(out_base_fields + extra_keys) + "\n")
        for fields, extra_dict in matched_rows:
            out_fields = [field for i, field in enumerate(fields) if i != extra_idx]
            out_fields += [extra_dict.get(key, "") for key in extra_keys]
            fout.write("\t".join(out_fields) + "\n")

    print(f"[OK] panel_rows={len(matched_rows)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
