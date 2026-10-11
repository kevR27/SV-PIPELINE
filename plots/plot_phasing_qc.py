#!/usr/bin/env python3
"""Summarize phased genotypes from WhatsHap and LongPhase VCFs."""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, read_tsv, save_figure, set_thesis_style, style_axis

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_snv_sv_candidates import pair_phase_relation, alt_haplotype


def parse_args():
    p = argparse.ArgumentParser(description="Plot phasing QC from phased VCFs.")
    p.add_argument("--whatshap-vcf", default=None, help="WhatsHap-phased small-variant VCF")
    p.add_argument("--longphase-vcf", default=None, help="LongPhase phased VCF")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--title", default="Phasing quality-control summary")
    p.add_argument("--pairs", help="Existing *_snv_sv_candidates.tsv")
    p.add_argument("--sv-candidates", help="SV table with coordinates and joint LongPhase GT/PS")
    p.add_argument("--top-n", type=int, default=10)
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


def draw_gene_phase(pairs, svs, out, top_n=10):
    """Draw only alleles whose phase can be verified from the two input tables."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    source_rows, manifest = [], []
    gene_col = "GENE" if "GENE" in svs else "GENES"
    for gene, group in list(pairs.groupby("GENE", sort=False))[:top_n]:
        drawable, unresolved = [], []
        for _, pair in group.iterrows():
            hit = svs[svs["SV_ID"].astype(str).eq(str(pair["SV_ID"])) & svs[gene_col].astype(str).eq(str(gene))]
            geometry = [col for col in ["CHROM", "START", "END", "LONGPHASE_GT", "LONGPHASE_PS", "LONGPHASE_PHASED"] if col in hit]
            consistent = not hit.empty and not hit[geometry].nunique(dropna=False).gt(1).any()
            sv = hit.iloc[0] if consistent else pd.Series(dtype=object)
            relation = pair_phase_relation(pair, sv)
            source_rows.append({**pair.to_dict(), "PLOT_PHASE_RELATION": relation, "SV_COORDINATES_CONSISTENT": consistent})
            if relation == "UNRESOLVED":
                unresolved.append(str(pair["SMALL_VARIANT"])+" + "+str(pair["SV_ID"]))
                continue
            drawable.append((pair, sv, relation))
        blocks = list(dict.fromkeys((str(pair["CHROM"]),str(pair["SMALL_PS"])) for pair, _sv, _relation in drawable))
        fig, axes = plt.subplots(max(1,len(blocks)), 1, figsize=(12, max(4, len(blocks)*3.2)), squeeze=False)
        for ax, block in zip(axes.flat, blocks):
            chrom, ps = block
            here = [(pair,sv,rel) for pair,sv,rel in drawable if (str(pair["CHROM"]),str(pair["SMALL_PS"])) == block]
            seen = set()
            for pair, sv, relation in here:
                for key, pos, end, gt, color, label in [
                    (str(pair["SV_ID"]), float(sv["START"]), float(sv["END"]) if sv.get("SVTYPE") in {"DEL","DUP","INV","CNV"} and pd.notna(sv.get("END")) else float(sv["START"]), str(sv["LONGPHASE_GT"]), SVTYPE_COLORS.get(str(sv["SVTYPE"]),"#999999"), str(sv["SVTYPE"])+" · "+str(pair["SV_ID"])),
                    (str(pair["SMALL_VARIANT"]), float(pair["POS"]), float(pair["POS"]), str(pair["SMALL_GT"]), "#222222", str(pair["SMALL_VARIANT"]))]:
                    if key in seen:
                        continue
                    seen.add(key)
                    hp = 1 if alt_haplotype(gt) == "H1" else 2
                    ax.plot([pos,end],[hp,hp],color=color,lw=6,solid_capstyle="butt")
                    ax.scatter([pos],[hp],color=color,s=40)
                    ax.annotate(label,(pos,hp),xytext=(0,9),textcoords="offset points",fontsize=8,rotation=15)
            ax.set_yticks([1,2],["HP1","HP2"])
            ax.set_ylim(.5,2.7)
            ax.set_title(f"{gene} · chromosome {chrom} · phase set {ps} · "+", ".join(sorted({rel for _,_,rel in here})),loc="left",fontsize=11)
            ax.set_xlabel("Reference position (bp)")
            ax.ticklabel_format(axis="x",style="plain",useOffset=False)
            ax.grid(axis="x",alpha=.15)
        if not blocks:
            axes.flat[0].axis("off")
            axes.flat[0].text(.5,.5,"Cis/trans unresolved\nNo verified same-chromosome, same-block phased pair",transform=axes.flat[0].transAxes,ha="center",fontsize=13)
        fig.suptitle(f"Allele positions for {gene}",fontweight="bold")
        fig.text(.5,.01,f"{len(unresolved)} unresolved pair(s) listed in source TSV. HP1/HP2 labels are local to each block, not maternal/paternal.\nThis shows called phase, not pathogenicity or proof of recessive disease. BNDs are shown at their local anchor only.",ha="center",fontsize=9)
        fig.tight_layout(rect=[0,.11,1,.94])
        import re
        stem = re.sub(r"[^A-Za-z0-9_.-]", "_", str(gene))
        save_figure(fig,out/stem)
        plt.close(fig)
        manifest.append({"GENE":gene,"FIGURE_PREFIX":str(out/stem),"RESOLVED_PAIRS":len(drawable),"UNRESOLVED_PAIRS":len(unresolved)})
    pd.DataFrame(source_rows).to_csv(out/"phase_pairs_source.tsv",sep="\t",index=False)
    pd.DataFrame(manifest,columns=["GENE","FIGURE_PREFIX","RESOLVED_PAIRS","UNRESOLVED_PAIRS"]).to_csv(out/"phase_gene_manifest.tsv",sep="\t",index=False)


def main():
    args = parse_args()
    set_thesis_style()
    if bool(args.pairs) != bool(args.sv_candidates):
        raise ValueError("Gene phase views need both --pairs and --sv-candidates")
    if args.pairs:
        draw_gene_phase(read_tsv(args.pairs), read_tsv(args.sv_candidates), Path(args.out_prefix).parent / "gene_phase", args.top_n)

    inputs = []
    if args.whatshap_vcf:
        inputs.append((args.whatshap_vcf, "WhatsHap"))
    if args.longphase_vcf:
        inputs.append((args.longphase_vcf, "LongPhase"))

    if not inputs and args.pairs:
        return
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

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))
    ax1, ax2 = axes

    ax1.bar(
        summary_df["source"],
        summary_df["fraction_het_phased"] * 100.0,
        color=["#0072B2", "#E69F00"][: len(summary_df)],
    )
    ax1.set_ylabel("Heterozygous genotypes phased (%)")
    ax1.set_ylim(0, 105)
    for i, value in enumerate(summary_df["fraction_het_phased"] * 100.0):
        if pd.notna(value):
            ax1.text(i, value + 1.5, f"{value:.1f}%", ha="center", va="bottom", fontsize=10)
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

    fig.suptitle(args.title, fontsize=16, fontweight="bold", y=0.995)
    fig.text(
        0.5,
        0.01,
        "This is phasing QC; candidate-locus haplotype interpretation should use the phased VCF/BAM evidence directly.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    outputs = save_figure(fig, prefix)
    plt.close(fig)
    print("[OK]", *outputs, sep="\n")


if __name__ == "__main__":
    main()
