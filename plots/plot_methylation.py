#!/usr/bin/env python3
"""Plot candidate-region methylation from modkit bedMethyl or extract output."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import first_existing, save_figure, set_thesis_style, style_axis

BEDMETHYL_COLUMNS = [
    "chrom", "start", "end", "mod_code", "score", "strand", "thick_start", "thick_end",
    "item_rgb", "valid_coverage", "fraction_modified", "n_modified", "n_canonical",
    "n_other_mod", "n_delete", "n_fail", "n_diff", "n_nocall",
]


def parse_args():
    p = argparse.ArgumentParser(description="Plot a candidate-region modkit methylation track.")
    p.add_argument("--input", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--region", default=None, help="Optional chr:start-end restriction")
    p.add_argument("--format", choices=["auto", "bedmethyl", "extract"], default="auto")
    p.add_argument("--min-coverage", type=int, default=5)
    p.add_argument("--title", default="Candidate-region CpG methylation")
    return p.parse_args()


def parse_region(text: str | None):
    if not text:
        return None
    chrom, coords = text.split(":", 1)
    start, end = coords.replace(",", "").split("-", 1)
    return chrom, int(start), int(end)


def looks_like_bedmethyl(path: Path) -> bool:
    with path.open(encoding="utf-8", errors="replace") as fh:
        first = fh.readline().rstrip("\n").split("\t")
    if len(first) < 11:
        return False
    try:
        int(first[1]); int(first[2]); float(first[10])
        return True
    except Exception:
        return False


def load_bedmethyl(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str, low_memory=False)
    n = min(df.shape[1], len(BEDMETHYL_COLUMNS))
    df = df.iloc[:, :n]
    df.columns = BEDMETHYL_COLUMNS[:n]
    for col in ["start", "end", "valid_coverage", "fraction_modified"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # modkit bedMethyl reports fraction_modified as percentage in current output.
    if "fraction_modified" in df.columns and df["fraction_modified"].dropna().max() <= 1.0:
        df["fraction_modified"] = df["fraction_modified"] * 100.0
    return df


def load_extract(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    chrom = first_existing(df, ["chrom", "chromosome", "ref_name", "contig"])
    pos = first_existing(df, ["ref_position", "position", "pos", "reference_position"])
    if chrom is None or pos is None:
        raise ValueError("Could not infer chromosome/reference-position columns from modkit extract output.")

    mod_prob = first_existing(df, ["mod_qual", "mod_probability", "probability", "mod_prob"])
    mod_code = first_existing(df, ["mod_code", "modified_base", "modification"])
    if mod_prob is None:
        raise ValueError(
            "Could not infer a modification probability/quality column from modkit extract output. "
            "When data are available, prefer `modkit pileup` bedMethyl for thesis methylation figures."
        )

    out = pd.DataFrame({
        "chrom": df[chrom].astype(str),
        "start": pd.to_numeric(df[pos], errors="coerce"),
        "prob": pd.to_numeric(df[mod_prob], errors="coerce"),
    })
    if mod_code:
        out["mod_code"] = df[mod_code].astype(str)

    # Values above 1 are commonly Phred-like quality rather than probabilities; do not silently reinterpret them.
    if out["prob"].dropna().max() > 1.0:
        raise ValueError(
            f"Detected values >1 in {mod_prob}. This looks like a quality score rather than a probability. "
            "Use modkit pileup bedMethyl or revise the mapping after inspecting the real extract file."
        )

    grouped = out.groupby(["chrom", "start"], as_index=False).agg(
        fraction_modified=("prob", lambda x: float(np.nanmean(x)) * 100.0),
        valid_coverage=("prob", "count"),
    )
    grouped["end"] = grouped["start"] + 1
    return grouped


def main():
    args = parse_args()
    set_thesis_style()
    path = Path(args.input)
    if not path.exists():
        raise FileNotFoundError(path)

    fmt = args.format
    if fmt == "auto":
        fmt = "bedmethyl" if looks_like_bedmethyl(path) else "extract"
    df = load_bedmethyl(path) if fmt == "bedmethyl" else load_extract(path)

    required = {"chrom", "start", "fraction_modified", "valid_coverage"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Methylation input is missing required normalized fields: {sorted(missing)}")

    region = parse_region(args.region)
    if region:
        chrom, start, end = region
        df = df[(df["chrom"].astype(str) == chrom) & (df["start"] >= start) & (df["start"] <= end)]
    df = df[pd.to_numeric(df["valid_coverage"], errors="coerce") >= args.min_coverage].copy()
    if df.empty:
        raise ValueError("No methylation positions remain after region/coverage filtering.")

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(prefix.with_name(prefix.name + "_plot_data.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(12.0, 6.6), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.0]})
    ax1, ax2 = axes

    for chrom, sub in df.groupby("chrom", sort=False):
        sub = sub.sort_values("start")
        ax1.plot(sub["start"], sub["fraction_modified"], marker="o", markersize=2.8, linewidth=0.9, label=str(chrom))
        ax2.fill_between(sub["start"], 0, pd.to_numeric(sub["valid_coverage"], errors="coerce"), alpha=0.5)

    ax1.set_ylabel("Modified calls (%)")
    ax1.set_ylim(-2, 102)
    ax1.set_title(args.title, fontweight="bold")
    style_axis(ax1, "y")
    ax2.set_ylabel("Coverage")
    ax2.set_xlabel("Genomic position (bp)")
    style_axis(ax2, "y")
    if df["chrom"].nunique() > 1:
        ax1.legend(frameon=False, fontsize=8)

    fig.text(0.5, 0.01, f"Input mode: {fmt}; minimum coverage: {args.min_coverage}. For final candidate figures, phased modkit pileup (HP1/HP2) is preferred when available.", ha="center", fontsize=8.1)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
