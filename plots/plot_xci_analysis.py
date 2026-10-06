#!/usr/bin/env python3
"""Create separate static figures for X-chromosome inactivation analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import save_figure, set_thesis_style, style_axis


def read_bedmethyl(path: str) -> pd.DataFrame:
    cols = [
        "chrom", "start", "end", "name", "score", "strand",
        "thick_start", "thick_end", "color", "valid_cov",
        "percent_modified",
    ]
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        dtype=str,
        compression="infer",
        low_memory=False,
    )
    if df.shape[1] < 11:
        raise ValueError(f"bedMethyl has fewer than 11 columns: {path}")
    df = df.iloc[:, :11].copy()
    df.columns = cols
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["valid_cov"] = pd.to_numeric(df["valid_cov"], errors="coerce")
    df["percent_modified"] = pd.to_numeric(
        df["percent_modified"],
        errors="coerce",
    )
    return df.dropna(subset=["start", "valid_cov", "percent_modified"])


def binned_methylation(df: pd.DataFrame, bin_bp: int):
    work = df.copy()
    work["bin_start"] = (work["start"] // bin_bp) * bin_bp
    work["_weighted"] = work["percent_modified"] * work["valid_cov"]
    out = (
        work.groupby("bin_start", as_index=False)
        .agg(
            weighted_sum=("_weighted", "sum"),
            coverage_sum=("valid_cov", "sum"),
        )
    )
    out["methylation_percent"] = np.where(
        out["coverage_sum"] > 0,
        out["weighted_sum"] / out["coverage_sum"],
        np.nan,
    )
    out["bin_mid_mb"] = (out["bin_start"] + bin_bp / 2) / 1e6
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--block-skew", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--phase-blocks", required=True)
    p.add_argument("--phase-summary", required=True)
    p.add_argument("--hp1-bedmethyl", required=True)
    p.add_argument("--hp2-bedmethyl", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--methylation-bin-bp", type=int, default=5_000_000)
    args = p.parse_args()

    set_thesis_style()
    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    blocks = pd.read_csv(args.block_skew, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    phase = pd.read_csv(args.phase_blocks, sep="\t")
    phase_summary = pd.read_csv(args.phase_summary, sep="\t")

    global_p = np.nan
    ratio = "."
    status = "."
    if not summary.empty:
        global_p = pd.to_numeric(
            pd.Series([summary.iloc[0].get("GLOBAL_FOLDED_SKEW_P", np.nan)]),
            errors="coerce",
        ).iloc[0]
        ratio = str(summary.iloc[0].get("XCI_MAJOR_MINOR_RATIO", "."))
        status = str(summary.iloc[0].get("XCI_SKEW_STATUS", "."))

    outputs = []

    # 1. Weighted folded block-skew distribution.
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    if blocks.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No informative phased XCI blocks",
            ha="center", va="center",
            transform=ax.transAxes,
        )
    else:
        vals = pd.to_numeric(
            blocks["FOLDED_BLOCK_SKEW"], errors="coerce"
        )
        weights = pd.to_numeric(blocks["TRIALS"], errors="coerce").fillna(0)
        valid = vals.notna() & weights.gt(0)
        ax.hist(
            vals[valid],
            bins=np.arange(0, 0.525, 0.025),
            weights=weights[valid],
        )
        if np.isfinite(global_p):
            ax.axvline(global_p, linestyle="--", linewidth=2)
        ax.axvline(0.5, linestyle=":", linewidth=1.5)
        ax.set_xlim(0, 0.5)
        ax.set_xlabel("Folded block skew P")
        ax.set_ylabel("Supporting reads")
        ax.set_title(
            f"{args.sample}: X-inactivation block skew | global {ratio}"
        )
        style_axis(ax, "both")
    fig.tight_layout()
    outputs += save_figure(fig, outdir / f"{args.sample}_xci_block_skew_distribution")
    plt.close(fig)

    # 2. Non-folded block skew along chrX.
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    if blocks.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No informative phased XCI blocks",
            ha="center", va="center",
            transform=ax.transAxes,
        )
    else:
        start = pd.to_numeric(blocks["BLOCK_START"], errors="coerce")
        end = pd.to_numeric(blocks["BLOCK_END"], errors="coerce")
        x = (start + end) / 2 / 1e6
        y = pd.to_numeric(blocks["H1_Xa_SKEW"], errors="coerce")
        trials = pd.to_numeric(blocks["TRIALS"], errors="coerce").fillna(1)
        sizes = np.clip(np.sqrt(trials) * 9, 15, 180)
        ax.scatter(x, y, s=sizes, alpha=0.7)
        ax.axhline(0.5, linestyle="--", linewidth=1.4)
        ax.set_ylim(0, 1)
        ax.set_xlabel("chrX position (Mb)")
        ax.set_ylabel("H1 preferential-Xa probability")
        ax.set_title(
            f"{args.sample}: local haplotype-block XCI orientation"
        )
        style_axis(ax, "both")
    fig.tight_layout()
    outputs += save_figure(fig, outdir / f"{args.sample}_xci_chrX_block_skew")
    plt.close(fig)

    # 3. Haplotype-specific methylation from modkit --phased.
    hp1 = binned_methylation(
        read_bedmethyl(args.hp1_bedmethyl),
        args.methylation_bin_bp,
    )
    hp2 = binned_methylation(
        read_bedmethyl(args.hp2_bedmethyl),
        args.methylation_bin_bp,
    )
    fig, ax = plt.subplots(figsize=(12.5, 5.8))
    if hp1.empty and hp2.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No haplotype-specific chrX methylation records",
            ha="center", va="center",
            transform=ax.transAxes,
        )
    else:
        if not hp1.empty:
            ax.plot(
                hp1["bin_mid_mb"],
                hp1["methylation_percent"],
                marker="o",
                label="WhatsHap HP1",
            )
        if not hp2.empty:
            ax.plot(
                hp2["bin_mid_mb"],
                hp2["methylation_percent"],
                marker="o",
                label="WhatsHap HP2",
            )
        ax.set_xlabel("chrX position (Mb)")
        ax.set_ylabel("CpG 5mC (%)")
        ax.set_ylim(0, 100)
        ax.set_title(
            f"{args.sample}: phased chrX CpG methylation"
        )
        ax.legend(frameon=False)
        style_axis(ax, "both")
    fig.tight_layout()
    outputs += save_figure(fig, outdir / f"{args.sample}_xci_haplotype_methylation")
    plt.close(fig)

    # 4. LongPhase / WhatsHap phase concordance.
    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    if phase.empty:
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No shared phased chrX SNVs between WhatsHap and LongPhase",
            ha="center", va="center",
            transform=ax.transAxes,
        )
    else:
        pview = phase.copy()
        pview["ORIENTATION_CONCORDANCE"] = pd.to_numeric(
            pview["ORIENTATION_CONCORDANCE"], errors="coerce"
        )
        pview["SHARED_PHASED_SNVS"] = pd.to_numeric(
            pview["SHARED_PHASED_SNVS"], errors="coerce"
        ).fillna(0)
        pview = pview.sort_values(
            ["SHARED_PHASED_SNVS", "ORIENTATION_CONCORDANCE"],
            ascending=[False, False],
        ).head(30)
        labels = [
            f"W{w}/L{l}"
            for w, l in zip(
                pview["WHATSHAP_PS"],
                pview["LONGPHASE_PS"],
            )
        ]
        ax.bar(
            np.arange(len(pview)),
            pview["ORIENTATION_CONCORDANCE"],
        )
        ax.set_xticks(np.arange(len(pview)))
        ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=8)
        ax.set_ylim(0.5, 1.02)
        ax.set_ylabel("Flip-tolerant phase concordance")
        ax.set_xlabel("Overlapping phase-block pair")
        concordance = "."
        if not phase_summary.empty:
            concordance = phase_summary.iloc[0].get(
                "FLIP_TOLERANT_PHASE_CONCORDANCE", "."
            )
        ax.set_title(
            f"{args.sample}: WhatsHap vs LongPhase chrX phasing | overall {concordance}"
        )
        style_axis(ax, "y")
    fig.tight_layout()
    outputs += save_figure(fig, outdir / f"{args.sample}_xci_phase_concordance")
    plt.close(fig)

    manifest = pd.DataFrame(
        [
            {
                "sample": args.sample,
                "XCI_SKEW_STATUS": status,
                "XCI_MAJOR_MINOR_RATIO": ratio,
                "GLOBAL_FOLDED_SKEW_P": (
                    round(float(global_p), 6)
                    if np.isfinite(global_p)
                    else "."
                ),
                "plot_directory": str(outdir),
            }
        ]
    )
    manifest.to_csv(
        outdir / f"{args.sample}_xci_plot_manifest.tsv",
        sep="\t",
        index=False,
    )

    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
