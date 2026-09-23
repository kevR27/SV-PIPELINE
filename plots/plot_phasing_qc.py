#!/usr/bin/env python3
"""Summarize phased genotypes from WhatsHap and LongPhase VCFs."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import save_figure, set_thesis_style, style_axis


def parse_args():
    p = argparse.ArgumentParser(description="Plot phasing QC from phased VCFs.")
    p.add_argument("--clair3-vcf", default=None, help="Phased Clair3 small-variant VCF")
    p.add_argument("--whatshap-vcf", default=None, help="Legacy alias for a WhatsHap-phased small-variant VCF")
    p.add_argument("--longphase-vcf", default=None, help="LongPhase phased VCF")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Phasing quality-control summary")
    return p.parse_args()


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def summarize_vcf(path: str, label: str) -> tuple[dict, pd.DataFrame]:
    counts = {
        "source": label,
        "total_genotyped": 0,
        "heterozygous": 0,
        "phased_heterozygous": 0,
        "with_phase_set": 0,
    }
    ps_counts: dict[str, int] = {}

    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            data = dict(zip(fmt, sample))
            gt = data.get("GT", ".")
            if gt in {".", "./.", ".|."}:
                continue

            counts["total_genotyped"] += 1
            alleles = gt.replace("|", "/").split("/")
            is_het = len(alleles) == 2 and alleles[0] != alleles[1] and "." not in alleles
            if is_het:
                counts["heterozygous"] += 1
                if "|" in gt:
                    counts["phased_heterozygous"] += 1

            ps = data.get("PS", ".")
            if ps not in {".", "", None}:
                counts["with_phase_set"] += 1
                key = f"{fields[0]}:{ps}"
                ps_counts[key] = ps_counts.get(key, 0) + 1

    counts["fraction_het_phased"] = (
        counts["phased_heterozygous"] / counts["heterozygous"]
        if counts["heterozygous"]
        else np.nan
    )
    ps = pd.DataFrame({"phase_set": list(ps_counts), "variant_count": list(ps_counts.values())})
    ps["source"] = label
    return counts, ps


def main():
    args = parse_args()
    set_thesis_style()

    inputs = []
    if args.clair3_vcf:
        inputs.append((args.clair3_vcf, "Clair3 phased"))
    elif args.whatshap_vcf:
        inputs.append((args.whatshap_vcf, "WhatsHap phased"))
    if args.longphase_vcf:
        inputs.append((args.longphase_vcf, "LongPhase"))

    if not inputs:
        raise ValueError("Provide --whatshap-vcf and/or --longphase-vcf.")

    summaries = []
    phase_sets = []
    for path, label in inputs:
        summary, ps = summarize_vcf(path, label)
        summaries.append(summary)
        phase_sets.append(ps)

    summary_df = pd.DataFrame(summaries)
    ps_df = pd.concat(phase_sets, ignore_index=True) if phase_sets else pd.DataFrame()

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)
    ps_df.to_csv(prefix.with_name(prefix.name + "_phase_sets.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8))
    ax1, ax2 = axes

    ax1.bar(
        summary_df["source"],
        summary_df["fraction_het_phased"] * 100.0,
        color=["#0072B2", "#E69F00"][: len(summary_df)],
    )
    ax1.set_ylabel("Heterozygous genotypes phased (%)")
    ax1.set_ylim(0, 105)
    style_axis(ax1, "y")

    if not ps_df.empty:
        grouped = list(ps_df.groupby("source", sort=False))
        groups = [g["variant_count"].values for _, g in grouped]
        labels = [name for name, _ in grouped]
        ax2.boxplot(groups, tick_labels=labels, showfliers=False)
        ax2.set_ylabel("Variants per phase set")
        ax2.set_yscale("log")
    else:
        ax2.text(0.5, 0.5, "No PS field detected", transform=ax2.transAxes, ha="center", va="center")
        ax2.set_ylabel("Variants per phase set")
    style_axis(ax2, "y")

    fig.suptitle(args.title, fontsize=14, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.01,
        "This is phasing QC; candidate-locus haplotype interpretation should use the phased VCF/BAM evidence directly.",
        ha="center",
        fontsize=8.2,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
