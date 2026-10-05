#!/usr/bin/env python3
"""Merge compact annotation delta tables into one final SV-gene table.

All annotation scripts receive the same compact ranked SV-gene table and emit
only newly added columns plus _INTEGRATED_ROW_INDEX. This merger validates a
strict one-to-one row relationship before attaching those annotations. It
therefore avoids serial full-table copies and prevents many-to-many row
expansion during MitoCarta, gnomAD-SV, complementary-evidence, and multimodal
context integration.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ROW_KEY = "_INTEGRATED_ROW_INDEX"


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)


def validate_delta(delta: pd.DataFrame, base_rows: int, label: str) -> None:
    if ROW_KEY not in delta.columns:
        raise ValueError(f"{label}: missing required {ROW_KEY} column")
    if len(delta) != base_rows:
        raise ValueError(
            f"{label}: row count differs from base table "
            f"({len(delta)} != {base_rows})"
        )
    if delta[ROW_KEY].duplicated().any():
        raise ValueError(f"{label}: duplicate {ROW_KEY} values detected")

    expected = set(map(str, range(base_rows)))
    observed = set(delta[ROW_KEY].astype(str))
    if observed != expected:
        missing = len(expected - observed)
        extra = len(observed - expected)
        raise ValueError(
            f"{label}: row-key coverage mismatch "
            f"(missing={missing}, unexpected={extra})"
        )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True)
    p.add_argument("--delta", action="append", default=[], help="Compact annotation delta table; may be repeated.")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    base = read_table(args.base).reset_index(drop=True)
    base[ROW_KEY] = base.index.astype(str)

    for delta_path in args.delta:
        delta = read_table(delta_path)
        label = Path(delta_path).name
        validate_delta(delta, len(base), label)

        delta[ROW_KEY] = delta[ROW_KEY].astype(str)
        new_columns = [c for c in delta.columns if c != ROW_KEY]
        overlaps = sorted(set(new_columns) & set(base.columns))
        if overlaps:
            raise ValueError(
                f"{label}: annotation columns already exist in base/final table: "
                + ", ".join(overlaps)
            )

        base = base.merge(
            delta,
            on=ROW_KEY,
            how="left",
            validate="one_to_one",
            sort=False,
        )

    base[ROW_KEY] = pd.to_numeric(base[ROW_KEY], errors="raise")
    base = base.sort_values(ROW_KEY).drop(columns=[ROW_KEY]).reset_index(drop=True)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    base.to_csv(output, sep="\t", index=False)

    print(
        f"[OK] final_integrated_rows={len(base)} columns={len(base.columns)} "
        f"deltas={len(args.delta)} output={output}"
    )


if __name__ == "__main__":
    main()
