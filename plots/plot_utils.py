#!/usr/bin/env python3
"""Shared utilities for thesis-quality SV pipeline figures."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MISSING_VALUES = {"", ".", "NA", "N/A", "nan", "NaN", "NAN", "None", None}

# Okabe-Ito-inspired, color-blind-friendly palette kept consistent across figures.
SVTYPE_COLORS = {
    "DEL": "#0072B2",
    "INS": "#E69F00",
    "DUP": "#009E73",
    "INV": "#CC79A7",
    "BND": "#D55E00",
    "CNV": "#56B4E9",
    "OTHER": "#999999",
}

CALLER_COLORS = {
    "Sniffles2": "#0072B2",
    "cuteSV": "#E69F00",
    "Delly": "#009E73",
}

SUPPORT_COLORS = {
    "SINGLE_CALLER": "#BDBDBD",
    "TWO_CALLER": "#56B4E9",
    "MULTICALLER": "#0072B2",
    "UNKNOWN": "#999999",
}


def set_thesis_style() -> None:
    """Apply one reproducible publication/thesis style to all figures."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Liberation Sans"],
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.bbox": "tight",
        }
    )


def read_tsv(path: str | Path, required: Iterable[str] | None = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    if required:
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {', '.join(missing)}. "
                f"Available columns: {', '.join(df.columns)}"
            )
    return df


def first_existing(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    for name in candidates:
        if name in df.columns:
            return name
        found = lower.get(name.lower())
        if found is not None:
            return found
    return None


def clean_string(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.replace(list(MISSING_VALUES), np.nan), errors="coerce")


def normalize_svtype(series: pd.Series) -> pd.Series:
    raw = clean_string(series).str.upper()
    out = pd.Series("OTHER", index=raw.index, dtype=object)
    for svtype in ["DEL", "INS", "DUP", "INV", "BND", "CNV"]:
        out.loc[raw.str.contains(svtype, regex=False, na=False)] = svtype
    return out


def abs_svlen(df: pd.DataFrame) -> pd.Series:
    col = first_existing(df, ["SVLEN", "SV_Length", "SV_length", "Length"])
    if col is not None:
        values = numeric(df[col]).abs()
    else:
        values = pd.Series(np.nan, index=df.index, dtype=float)

    start_col = first_existing(df, ["START", "POS", "SV_start", "Start"])
    end_col = first_existing(df, ["END", "SV_end", "End"])
    if start_col and end_col:
        fallback = (numeric(df[end_col]) - numeric(df[start_col])).abs()
        values = values.fillna(fallback)
    return values


def unique_master_svs(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse integrated one-row-per-(SV,gene) tables to one row per master SV."""
    id_col = first_existing(df, ["SV_ID", "ID", "AnnotSV_ID"])
    if id_col is None:
        raise ValueError("Could not identify an SV ID column for deduplication.")
    return df.drop_duplicates(subset=[id_col], keep="first").copy()


def chromosome_sort_key(value: str):
    x = str(value).replace("chr", "", 1)
    if x.isdigit():
        return (0, int(x))
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, order.get(x.upper(), 1000)) if x.upper() in order else (1, x)


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="right",
    )


def style_axis(ax, grid_axis: str = "y") -> None:
    ax.grid(axis=grid_axis, color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines["left"].set_color("#BDBDBD")
    ax.spines["bottom"].set_color("#BDBDBD")


def save_figure(fig, out_prefix: str | Path, dpi: int = 600) -> list[Path]:
    """Save vector PDF/SVG plus a high-resolution PNG."""
    prefix = Path(out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for ext, kwargs in [
        ("pdf", {}),
        ("svg", {}),
        ("png", {"dpi": dpi}),
    ]:
        path = prefix.with_suffix(f".{ext}")
        fig.savefig(path, **kwargs)
        outputs.append(path)
    return outputs


def safe_log10(values: pd.Series) -> pd.Series:
    values = numeric(values)
    values = values.where(values > 0)
    return np.log10(values)


def split_tokens(value) -> list[str]:
    if value in MISSING_VALUES:
        return []
    text = str(value)
    for sep in [";", "|", ","]:
        text = text.replace(sep, " ")
    return [x for x in text.split() if x]


def parse_semicolon_numeric(value) -> list[float]:
    if value in MISSING_VALUES:
        return []
    values = []
    for token in str(value).replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            x = float(token)
        except ValueError:
            continue
        if math.isfinite(x):
            values.append(x)
    return values
