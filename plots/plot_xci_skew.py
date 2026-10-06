#!/usr/bin/env python3
"""Create separate static thesis figures for X-chromosome inactivation."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_tsv(path):
    try:
        return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def save(fig, prefix):
    outputs = []
    for suffix in ("pdf", "svg", "png"):
        path = Path(str(prefix) + "." + suffix)
        fig.savefig(path, dpi=300, bbox_inches="tight")
        outputs.append(path)
    return outputs


def empty_panel(title, message, prefix):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.axis("off")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.text(
        0.5,
        0.5,
        message,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=10,
    )
    out = save(fig, prefix)
    plt.close(fig)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--block-skew", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--phase-detail", required=True)
    p.add_argument("--phase-summary", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    blocks = read_tsv(args.block_skew)
    summary = read_tsv(args.summary)
    phase = read_tsv(args.phase_detail)
    phase_summary = read_tsv(args.phase_summary)

    # 1. Folded XCI skew across chrX phase blocks.
    block_prefix = outdir / f"{args.sample}_xci_block_skew"
    if blocks.empty:
        empty_panel(
            "X-chromosome inactivation skew by phase block",
            "No informative methylation-phased chrX blocks were available.",
            block_prefix,
        )
    else:
        work = blocks.copy()
        work["BLOCK_START"] = numeric(work["BLOCK_START"])
        work["MINOR"] = numeric(work["FOLDED_MINOR_XI_FRACTION"])
        work["READS"] = numeric(work["TOTAL_READS"]).fillna(0)
        work = work[
            work["BLOCK_START"].notna() & work["MINOR"].notna()
        ].sort_values("BLOCK_START")

        fig, ax = plt.subplots(figsize=(12.5, 6.2))
        if work.empty:
            ax.axis("off")
            ax.text(
                0.5, 0.5,
                "No informative phase blocks after filtering.",
                transform=ax.transAxes,
                ha="center", va="center",
            )
        else:
            x = work["BLOCK_START"] / 1e6
            sizes = np.clip(18 + 2.2 * np.sqrt(work["READS"]), 22, 110)
            ax.scatter(x, work["MINOR"], s=sizes, alpha=0.75)
            ax.axhline(0.5, linestyle="--", linewidth=1)
            ax.axhline(0.2, linestyle=":", linewidth=1)
            ax.axhline(0.1, linestyle=":", linewidth=1)
            ax.set_ylim(-0.02, 0.52)
            ax.set_xlabel("chrX position of phase block (Mb)")
            ax.set_ylabel("Minor inactive-X fraction")
            ax.set_title(
                f"{args.sample}: XCI skew across phased methylation blocks"
            )
            ax.text(
                0.99, 0.205, "80:20",
                transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8,
            )
            ax.text(
                0.99, 0.105, "90:10",
                transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=8,
            )

            if not summary.empty:
                mle = pd.to_numeric(
                    summary.iloc[0].get(
                        "XCI_MINOR_INACTIVE_FRACTION_MLE", np.nan
                    ),
                    errors="coerce",
                )
                if pd.notna(mle):
                    ax.axhline(mle, linewidth=1.5)
                    ax.text(
                        0.01,
                        0.04,
                        f"Sample folded MLE = {mle:.3f}:{1-mle:.3f}",
                        transform=ax.transAxes,
                        fontsize=9,
                    )

            ax.grid(axis="y", alpha=0.2)
        save(fig, block_prefix)
        plt.close(fig)

    # 2. Xa/Xi read counts by WhatsHap haplotype.
    count_prefix = outdir / f"{args.sample}_xci_haplotype_xa_xi_counts"
    if blocks.empty:
        empty_panel(
            "Xa/Xi methylation-state counts by haplotype",
            "No informative phased methylation reads were available.",
            count_prefix,
        )
    else:
        work = blocks.copy()
        count_cols = ["H1_Xa", "H2_Xa", "H1_Xi", "H2_Xi"]
        for col in count_cols:
            work[col] = numeric(work[col]).fillna(0)
        work["BLOCK_START"] = numeric(work["BLOCK_START"])
        work = work.sort_values("BLOCK_START").reset_index(drop=True)

        # The figure remains readable for large numbers of blocks by using
        # genomic order rather than long phase-set labels on the x-axis.
        x = np.arange(len(work))
        fig, ax = plt.subplots(figsize=(max(11, len(work) * 0.34), 6.5))
        bottom = np.zeros(len(work))
        for col in count_cols:
            values = work[col].to_numpy(dtype=float)
            ax.bar(x, values, bottom=bottom, label=col.replace("_", " "))
            bottom += values

        ax.set_xlabel("Phase blocks ordered along chrX")
        ax.set_ylabel("Read count")
        ax.set_title(
            f"{args.sample}: phased reads classified as Xa-like or Xi-like"
        )
        ax.legend(frameon=False, ncol=4, fontsize=8)
        ax.grid(axis="y", alpha=0.2)
        save(fig, count_prefix)
        plt.close(fig)

    # 3. WhatsHap versus LongPhase phase consistency.
    phase_prefix = outdir / f"{args.sample}_xci_phase_concordance"
    if phase.empty:
        message = "No shared informative chrX phase blocks between WhatsHap and LongPhase."
        if not phase_summary.empty:
            message += (
                "\nStatus: "
                + str(
                    phase_summary.iloc[0].get(
                        "PHASE_CONCORDANCE_STATUS", "."
                    )
                )
            )
        empty_panel(
            "WhatsHap versus LongPhase chrX phase concordance",
            message,
            phase_prefix,
        )
    else:
        work = phase.copy()
        work["CONCORDANCE"] = numeric(work["PHASE_CONCORDANCE"])
        work["BLOCK_START"] = numeric(work["BLOCK_START"])
        work["N"] = numeric(work["SHARED_PHASED_HET_SNPS"]).fillna(0)
        work = work[
            work["CONCORDANCE"].notna()
        ].sort_values("BLOCK_START")

        fig, ax = plt.subplots(figsize=(12.5, 5.8))
        if work.empty:
            ax.axis("off")
            ax.text(
                0.5, 0.5,
                "No informative phase-concordance blocks.",
                transform=ax.transAxes,
                ha="center", va="center",
            )
        else:
            x = work["BLOCK_START"] / 1e6
            sizes = np.clip(18 + 3 * np.sqrt(work["N"]), 22, 100)
            ax.scatter(x, 100 * work["CONCORDANCE"], s=sizes, alpha=0.8)
            ax.axhline(90, linestyle="--", linewidth=1)
            ax.set_ylim(0, 102)
            ax.set_xlabel("chrX block position (Mb)")
            ax.set_ylabel("Phase concordance after block flip (%)")
            ax.set_title(
                f"{args.sample}: WhatsHap–LongPhase phase consistency"
            )
            ax.grid(axis="y", alpha=0.2)
        save(fig, phase_prefix)
        plt.close(fig)

    print(f"[OK] XCI plots written to {outdir}")


if __name__ == "__main__":
    main()
