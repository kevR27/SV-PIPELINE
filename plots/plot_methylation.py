#!/usr/bin/env python3
"""Plot modkit methylation from compressed/uncompressed bedMethyl files.

Without --region, generate a genome-wide summary using chunked reading.
With --region, generate a locus-level methylation/coverage track.
"""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import save_figure, set_thesis_style, style_axis

BEDMETHYL_COLUMNS = [
    "chrom", "start", "end", "mod_code", "score", "strand", "thick_start", "thick_end",
    "item_rgb", "valid_coverage", "fraction_modified", "n_modified", "n_canonical",
    "n_other_mod", "n_delete", "n_fail", "n_diff", "n_nocall",
]

CANONICAL_RE = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", re.IGNORECASE)


def parse_args():
    p = argparse.ArgumentParser(description="Plot modkit bedMethyl output.")
    p.add_argument("--input", required=True, help="modkit bedMethyl, .bed.gz, or .bedmethyl.gz")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--region", default=None, help="Optional chr:start-end candidate region")
    p.add_argument("--min-coverage", type=int, default=5)
    p.add_argument("--chunksize", type=int, default=500000)
    p.add_argument("--title", default="CpG methylation summary")
    return p.parse_args()


def parse_region(text):
    if not text:
        return None
    chrom, coords = text.split(":", 1)
    start, end = coords.replace(",", "").split("-", 1)
    return chrom, int(start), int(end)


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def looks_like_bedmethyl(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                return False
            try:
                int(fields[1])
                int(fields[2])
                float(fields[9])
                float(fields[10])
                return True
            except Exception:
                return False
    return False


def normalize_mod_code(value):
    text = str(value).strip().lower()
    if text in {"m", "5mc", "c+m"}:
        return "5mC"
    if text in {"h", "5hmc", "c+h"}:
        return "5hmC"
    return str(value).strip() or "other"


def detect_scale(path):
    vals = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                vals.append(float(fields[10]))
            except Exception:
                continue
            if len(vals) >= 1000:
                break
    return 100.0 if vals and max(vals) <= 1.0 else 1.0


def load_region(path, region, min_cov, scale):
    chrom, start, end = region
    try:
        import pysam
        idx = Path(str(path) + ".tbi")
        csi = Path(str(path) + ".csi")
        if idx.exists() or csi.exists():
            tbx = pysam.TabixFile(str(path))
            rows = []
            try:
                fetched = tbx.fetch(chrom, max(0, start), end + 1)
                for line in fetched:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 11:
                        continue
                    try:
                        cov = float(fields[9])
                        frac = float(fields[10]) * scale
                        pos = int(fields[1])
                        stop = int(fields[2])
                    except Exception:
                        continue
                    if cov < min_cov:
                        continue
                    rows.append(
                        {
                            "chrom": chrom,
                            "start": pos,
                            "end": stop,
                            "mod_code": normalize_mod_code(fields[3]),
                            "valid_coverage": cov,
                            "fraction_modified": frac,
                        }
                    )
            finally:
                tbx.close()
            return pd.DataFrame(rows)
    except Exception:
        pass

    rows = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        names=BEDMETHYL_COLUMNS,
        usecols=range(len(BEDMETHYL_COLUMNS)),
        dtype=str,
        compression="infer",
        chunksize=250000,
        low_memory=False,
    ):
        chunk["start"] = pd.to_numeric(chunk["start"], errors="coerce")
        chunk["end"] = pd.to_numeric(chunk["end"], errors="coerce")
        chunk["valid_coverage"] = pd.to_numeric(chunk["valid_coverage"], errors="coerce")
        chunk["fraction_modified"] = pd.to_numeric(chunk["fraction_modified"], errors="coerce") * scale
        sub = chunk[
            (chunk["chrom"].astype(str) == chrom)
            & (chunk["start"] <= end)
            & (chunk["end"] >= start)
            & (chunk["valid_coverage"] >= min_cov)
        ].copy()
        if not sub.empty:
            sub["mod_code"] = sub["mod_code"].map(normalize_mod_code)
            rows.append(sub[["chrom", "start", "end", "mod_code", "valid_coverage", "fraction_modified"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_region(df, prefix, region, title, min_cov):
    if df.empty:
        raise ValueError("No methylation positions remain in the requested region after coverage filtering.")

    df = df.sort_values(["mod_code", "start"])
    df.to_csv(prefix.with_name(prefix.name + "_plot_data.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(
        2, 1, figsize=(14.0, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    ax1, ax2 = axes

    for mod, sub in df.groupby("mod_code", sort=False):
        ax1.plot(
            sub["start"],
            sub["fraction_modified"],
            marker="o",
            markersize=3.0,
            linewidth=1.0,
            label=str(mod),
        )
    coverage = df.groupby("start", as_index=False)["valid_coverage"].max()
    ax2.fill_between(coverage["start"], 0, coverage["valid_coverage"], alpha=0.45)

    ax1.set_ylabel("Modified calls (%)")
    ax1.set_ylim(-2, 102)
    ax1.legend(frameon=False, ncol=3)
    style_axis(ax1, "y")
    ax2.set_ylabel("Coverage")
    ax2.set_xlabel(f"Genomic position in {region[0]}")
    style_axis(ax2, "y")

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5, 0.008,
        f"Regional modkit bedMethyl view; minimum coverage = {min_cov}. Methylation is contextual evidence, not an SV confirmation.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    return outputs


def genome_summary(path, min_cov, chunksize, scale):
    stats = defaultdict(lambda: {
        "n": 0,
        "sum_pct": 0.0,
        "sum_cov": 0.0,
        "hist": np.zeros(20, dtype=np.int64),
    })

    for chunk in pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        names=BEDMETHYL_COLUMNS,
        usecols=range(len(BEDMETHYL_COLUMNS)),
        dtype=str,
        compression="infer",
        chunksize=chunksize,
        low_memory=False,
    ):
        chunk["valid_coverage"] = pd.to_numeric(chunk["valid_coverage"], errors="coerce")
        chunk["fraction_modified"] = pd.to_numeric(chunk["fraction_modified"], errors="coerce") * scale
        chunk = chunk[
            chunk["chrom"].astype(str).str.match(CANONICAL_RE)
            & (chunk["valid_coverage"] >= min_cov)
            & chunk["fraction_modified"].notna()
        ].copy()
        if chunk.empty:
            continue
        chunk["mod_code"] = chunk["mod_code"].map(normalize_mod_code)

        for (chrom, mod), sub in chunk.groupby(["chrom", "mod_code"]):
            vals = sub["fraction_modified"].clip(lower=0, upper=100).to_numpy(float)
            cov = sub["valid_coverage"].to_numpy(float)
            key = (chrom, mod)
            stats[key]["n"] += len(sub)
            stats[key]["sum_pct"] += float(np.nansum(vals))
            stats[key]["sum_cov"] += float(np.nansum(cov))
            hist, _ = np.histogram(vals, bins=np.linspace(0, 100, 21))
            stats[key]["hist"] += hist

    rows = []
    for (chrom, mod), st in stats.items():
        if st["n"] == 0:
            continue
        rows.append(
            {
                "chrom": chrom,
                "mod_code": mod,
                "n_records": st["n"],
                "mean_modified_percent": st["sum_pct"] / st["n"],
                "mean_coverage": st["sum_cov"] / st["n"],
            }
        )
    return pd.DataFrame(rows), stats


def chrom_sort_key(chrom):
    x = str(chrom).replace("chr", "", 1).upper()
    if x.isdigit():
        return int(x)
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(x, 99)


def plot_genome(summary, stats, prefix, title, min_cov):
    if summary.empty:
        raise ValueError("No canonical-chromosome methylation records remain after coverage filtering.")

    summary = summary.sort_values("chrom", key=lambda s: s.map(chrom_sort_key))
    summary.to_csv(prefix.with_name(prefix.name + "_genome_summary.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(14.2, 9.0))
    ax1, ax2, ax3, ax4 = axes.flatten()

    chroms = sorted(summary["chrom"].unique(), key=chrom_sort_key)
    mods = summary["mod_code"].drop_duplicates().tolist()

    width = 0.8 / max(len(mods), 1)
    x = np.arange(len(chroms))
    for i, mod in enumerate(mods):
        sub = summary[summary["mod_code"] == mod].set_index("chrom").reindex(chroms)
        ax1.bar(x + (i - (len(mods)-1)/2) * width, sub["mean_modified_percent"], width=width, label=mod)
    ax1.set_xticks(x)
    ax1.set_xticklabels(chroms, rotation=45, ha="right")
    ax1.set_ylabel("Mean modified calls (%)")
    ax1.set_xlabel("Canonical chromosome")
    ax1.legend(frameon=False)
    style_axis(ax1, "y")

    bins = np.linspace(0, 100, 21)
    centers = (bins[:-1] + bins[1:]) / 2
    for mod in mods:
        total = np.zeros(20, dtype=np.int64)
        for (chrom, m), st in stats.items():
            if m == mod:
                total += st["hist"]
        ax2.plot(centers, total, marker="o", linewidth=1.5, markersize=3, label=mod)
    ax2.set_xlabel("Modified calls (%)")
    ax2.set_ylabel("CpG records")
    ax2.legend(frameon=False)
    style_axis(ax2, "both")

    coverage = summary.groupby("chrom", as_index=False)["mean_coverage"].mean()
    coverage = coverage.sort_values("chrom", key=lambda s: s.map(chrom_sort_key))
    ax3.bar(coverage["chrom"], coverage["mean_coverage"])
    ax3.tick_params(axis="x", rotation=45)
    ax3.set_xlabel("Canonical chromosome")
    ax3.set_ylabel("Mean CpG coverage")
    style_axis(ax3, "y")

    counts = summary.groupby("mod_code", as_index=False)["n_records"].sum().sort_values("n_records", ascending=False)
    ax4.bar(counts["mod_code"], counts["n_records"])
    ax4.set_xlabel("Modification class")
    ax4.set_ylabel("CpG records")
    style_axis(ax4, "y")

    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5, 0.008,
        f"Genome-wide modkit bedMethyl summary on canonical chromosomes; minimum coverage = {min_cov}.",
        ha="center", fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    return outputs


def main():
    args = parse_args()
    set_thesis_style()
    path = Path(args.input)
    if not path.exists():
        raise FileNotFoundError(path)
    if not looks_like_bedmethyl(path):
        raise ValueError(
            "Input does not look like bedMethyl. Compressed files may be named .bed.gz, "
            "but they still need standard bedMethyl columns."
        )

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    scale = detect_scale(path)
    region = parse_region(args.region)

    if region:
        df = load_region(path, region, args.min_coverage, scale)
        outputs = plot_region(df, prefix, region, args.title, args.min_coverage)
    else:
        summary, stats = genome_summary(path, args.min_coverage, args.chunksize, scale)
        outputs = plot_genome(summary, stats, prefix, "Genome-wide CpG methylation summary", args.min_coverage)

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
