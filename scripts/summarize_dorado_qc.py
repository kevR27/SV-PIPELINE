#!/usr/bin/env python3
"""Summarize read-level Dorado sequencing-summary output."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def n50(lengths: pd.Series) -> int:
    values = pd.to_numeric(lengths, errors="coerce").dropna().astype(int)
    if values.empty:
        return 0
    values = values.sort_values(ascending=False)
    return int(values.iloc[np.searchsorted(values.cumsum().to_numpy(), values.sum() / 2)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    length_col = "sequence_length_template"
    q_col = "mean_qscore_template"

    if length_col not in df:
        raise ValueError(f"Dorado summary is missing {length_col}")

    lengths = pd.to_numeric(df[length_col], errors="coerce")
    qscores = pd.to_numeric(df[q_col], errors="coerce") if q_col in df else pd.Series(dtype=float)
    duration = pd.to_numeric(df["duration"], errors="coerce") if "duration" in df else pd.Series(dtype=float)

    row = {
        "SAMPLE": args.sample,
        "READ_COUNT": int(lengths.notna().sum()),
        "TOTAL_BASES": int(lengths.sum()),
        "YIELD_GBP": round(float(lengths.sum()) / 1e9, 3),
        "MEAN_READ_LENGTH": round(float(lengths.mean()), 1),
        "MEDIAN_READ_LENGTH": round(float(lengths.median()), 1),
        "READ_N50": n50(lengths),
        "MEAN_QSCORE": round(float(qscores.mean()), 3) if not qscores.empty else ".",
        "MEDIAN_QSCORE": round(float(qscores.median()), 3) if not qscores.empty else ".",
        "MEDIAN_DURATION_SECONDS": round(float(duration.median()), 3) if not duration.empty else ".",
        "CHANNELS_USED": int(df["channel"].nunique()) if "channel" in df else ".",
    }

    if "passes_filtering" in df:
        passed = df["passes_filtering"].astype(str).str.upper().eq("TRUE")
        row["PASSING_READ_PERCENT"] = round(100 * float(passed.mean()), 2)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(output, sep="\t", index=False)
    print(f"[OK] dorado_qc={output}")


if __name__ == "__main__":
    main()
