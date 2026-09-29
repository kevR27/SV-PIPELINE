#!/usr/bin/env python3
"""Create one detailed locus figure per prioritized SV-gene pair.

The locus view links a master SV to its overlapping gene, nearby genes, optional
indexed modkit methylation calls and the principal evidence fields carried by
the integrated table.  It is an interpretation aid, not read-level validation.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, first_existing, normalize_svtype, numeric, read_tsv, save_figure, set_thesis_style


def parse_args():
    p = argparse.ArgumentParser(description="Generate detailed candidate-specific SV locus figures.")
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--gene-bed", required=True, help="BED with chrom, start, end, gene label")
    p.add_argument("--methylation-bed", default=None, help="Optional tabix-indexed modkit bedMethyl")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=8)
    p.add_argument("--flank", type=int, default=50000)
    p.add_argument("--title-prefix", default="Candidate locus")
    return p.parse_args()


def normalize_chrom(value):
    value = str(value)
    return value if value.startswith("chr") else "chr" + value


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "candidate"


def load_gene_bed(path):
    bed = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str, low_memory=False)
    if bed.shape[1] < 3:
        raise ValueError("Gene BED needs at least chromosome, start and end columns.")
    bed = bed.iloc[:, : max(4, min(6, bed.shape[1]))].copy()
    bed["chrom"] = bed.iloc[:, 0].map(normalize_chrom)
    bed["start"] = pd.to_numeric(bed.iloc[:, 1], errors="coerce")
    bed["end"] = pd.to_numeric(bed.iloc[:, 2], errors="coerce")
    if bed.shape[1] >= 4:
        bed["gene"] = bed.iloc[:, 3].fillna(".").astype(str)
    else:
        bed["gene"] = "."
    return bed.dropna(subset=["start", "end"])


def methylation_records(path, chrom, start, end):
    if not path:
        return []
    try:
        import pysam
    except ImportError:
        return []
    p = Path(path)
    if not p.exists() or not (Path(str(p) + ".tbi").exists() or Path(str(p) + ".csi").exists()):
        return []
    try:
        tbx = pysam.TabixFile(str(p))
    except Exception:
        return []
    rows = []
    try:
        for line in tbx.fetch(chrom, max(0, int(start)), int(end) + 1):
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            try:
                pos = (int(f[1]) + int(f[2])) / 2
                cov = float(f[9])
                pct = float(f[10])
            except ValueError:
                continue
            if np.isfinite(pct) and np.isfinite(cov):
                rows.append((pos, pct, cov, f[3]))
    except Exception:
        pass
    finally:
        tbx.close()
    return rows


def evidence_value(row, candidates):
    for col in candidates:
        if col in row.index:
            value = row[col]
            if pd.notna(value) and str(value).strip() not in {"", ".", "NA", "N/A", "nan", "None"}:
                return str(value)
    return "."


def compact_text(value, max_chars=72, max_tokens=4):
    text = str(value)
    if text in {"", ".", "NA", "N/A", "nan", "None"}:
        return "."
    tokens = [x.strip() for x in text.split(";") if x.strip()]
    if len(tokens) > max_tokens:
        text = ";".join(tokens[:max_tokens]) + f";(+{len(tokens) - max_tokens} more)"
    if len(text) > max_chars:
        text = text[: max_chars - 3] + "..."
    return text


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    genes = load_gene_bed(args.gene_bed)

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    end_col = first_existing(df, ["END"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    score_col = first_existing(df, ["INTEGRATED_DISCOVERY_SCORE", "integrated_discovery_score", "PHENOTYPE_SCORE", "ALLELE_RESEARCH_SCORE"])
    caller_col = first_existing(df, ["CALLER_COUNT", "SUPP"])
    af_col = first_existing(df, ["NEEDLR_AF"])

    if None in (id_col, gene_col, chrom_col, start_col):
        raise ValueError("Input needs SV_ID, gene, chromosome and start columns.")

    work = df.copy()
    work["_gene"] = work[gene_col].fillna(".").astype(str)
    work["_start"] = numeric(work[start_col])
    work["_end"] = numeric(work[end_col]) if end_col else work["_start"]
    work["_end"] = work["_end"].fillna(work["_start"])
    work["_priority"] = numeric(work[score_col]).fillna(0) if score_col else 0
    work["_caller_count"] = numeric(work[caller_col]).fillna(0) if caller_col else 0
    work["_svtype"] = normalize_svtype(work[type_col]) if type_col else "OTHER"
    work["_af"] = numeric(work[af_col]) if af_col else np.nan
    work["_population_rank"] = np.select(
        [work["_af"].notna() & work["_af"].le(0.01), work["_af"].isna()],
        [2, 1],
        default=0,
    )

    work = work[~work["_gene"].isin(["", ".", "NA", "N/A", "nan", "None"])].copy()
    work = (
        work.sort_values(
            ["_priority", "_population_rank", "_caller_count", id_col, "_gene"],
            ascending=[False, False, False, True, True],
        )
        .drop_duplicates([id_col, gene_col], keep="first")
        .head(args.top_n)
        .copy()
    )

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for rank, (_, row) in enumerate(work.iterrows(), 1):
        sv_id = str(row[id_col])
        gene = str(row["_gene"])
        chrom = normalize_chrom(row[chrom_col])
        start = int(row["_start"])
        end = int(row["_end"])
        if end < start:
            start, end = end, start
        svtype = str(row["_svtype"])
        is_breakpoint = svtype in {"INS", "BND"} or end == start

        sv_span = abs(end - start)
        target_gene = genes[
            genes["chrom"].eq(chrom)
            & genes["gene"].astype(str).eq(gene)
        ].copy()
        large_gene_centered = sv_span >= 1_000_000 and not target_gene.empty

        if large_gene_centered:
            gene_start = int(target_gene["start"].min())
            gene_end = int(target_gene["end"].max())
            locus_start = max(0, gene_start - args.flank)
            locus_end = gene_end + args.flank
        else:
            locus_start = max(0, start - args.flank)
            locus_end = max(end, start) + args.flank

        nearby = genes[
            genes["chrom"].eq(chrom)
            & (genes["end"] >= locus_start)
            & (genes["start"] <= locus_end)
        ].copy()
        if len(nearby) > 24:
            nearby["_distance"] = np.minimum(
                (nearby["start"] - start).abs(),
                (nearby["end"] - start).abs(),
            )
            keep_target = nearby["gene"].astype(str).eq(gene)
            nearby = pd.concat(
                [nearby[keep_target], nearby[~keep_target].nsmallest(23, "_distance")]
            ).drop_duplicates().sort_values("start")

        methyl = methylation_records(args.methylation_bed, chrom, locus_start, locus_end)

        fig, axes = plt.subplots(
            3,
            1,
            figsize=(14.5, 9.4),
            gridspec_kw={"height_ratios": [1.1, 1.5, 1.25]},
        )
        ax_sv, ax_gene, ax_ev = axes

        color = SVTYPE_COLORS.get(svtype, "#999999")
        if large_gene_centered:
            ax_sv.plot(
                [locus_start, locus_end],
                [0.5, 0.5],
                color=color,
                linewidth=9,
                alpha=0.45,
                solid_capstyle="butt",
            )
            if locus_start <= start <= locus_end:
                ax_sv.axvline(start, color=color, linewidth=3)
            if locus_start <= end <= locus_end:
                ax_sv.axvline(end, color=color, linewidth=3)
            scope_note = f"Gene-centered view within {sv_span / 1e6:.1f} Mb {svtype}"
            if svtype == "INV":
                scope_note += "; interval overlap does not imply breakpoint disruption"
            ax_sv.text(
                0.5, 0.78, scope_note,
                transform=ax_sv.transAxes,
                ha="center", va="center", fontsize=9,
            )
        elif is_breakpoint:
            ax_sv.axvline(start, color=color, linewidth=4)
            ax_sv.scatter([start], [0.5], s=90, color=color, zorder=3)
        else:
            ax_sv.plot([start, end], [0.5, 0.5], color=color, linewidth=9, solid_capstyle="butt")
            ax_sv.scatter([start, end], [0.5, 0.5], s=35, color=color, zorder=3)
        ax_sv.set_xlim(locus_start, locus_end)
        ax_sv.set_ylim(0, 1)
        ax_sv.set_yticks([])
        ax_sv.set_ylabel("SV")
        ax_sv.set_title(f"{svtype} {chrom}:{start:,}-{end:,} | {sv_id}", loc="left", fontsize=11)
        ax_sv.grid(axis="x", color="#E6E6E6", linewidth=0.7)

        if not nearby.empty:
            y_positions = np.arange(len(nearby))
            for y, (_, g) in zip(y_positions, nearby.iterrows()):
                target = str(g["gene"]) == gene
                face = "#0072B2" if target else "#C7CDD3"
                ax_gene.add_patch(
                    Rectangle(
                        (float(g["start"]), y - 0.28),
                        max(float(g["end"] - g["start"]), 1.0),
                        0.56,
                        facecolor=face,
                        edgecolor="none",
                    )
                )
                ax_gene.text(float(g["start"]), y + 0.34, str(g["gene"]), fontsize=8.5, va="bottom")
            ax_gene.set_ylim(-0.7, len(nearby) - 0.2)
        else:
            ax_gene.text(0.5, 0.5, "No gene BED entries in displayed interval", transform=ax_gene.transAxes, ha="center")
            ax_gene.set_ylim(0, 1)
        ax_gene.set_xlim(locus_start, locus_end)
        ax_gene.set_yticks([])
        ax_gene.set_ylabel("Genes")
        ax_gene.grid(axis="x", color="#E6E6E6", linewidth=0.7)

        if methyl:
            inset = ax_gene.inset_axes([0.0, -0.55, 1.0, 0.36], transform=ax_gene.transAxes)
            pos = [x[0] for x in methyl]
            pct = [x[1] for x in methyl]
            cov = [x[2] for x in methyl]
            inset.scatter(pos, pct, s=np.clip(np.sqrt(cov) * 4, 6, 36), alpha=0.65)
            inset.set_xlim(locus_start, locus_end)
            inset.set_ylim(0, 100)
            inset.set_ylabel("5mC %", fontsize=8)
            inset.tick_params(labelsize=7)
            inset.grid(axis="y", color="#EEEEEE", linewidth=0.6)

        evidence = [
            ("Overlapping gene", gene),
            ("Gene-discovery score", f"{row['_priority']:.3g}" if pd.notna(row["_priority"]) else "."),
            ("Callers", evidence_value(row, ["CALLERS"])),
            ("Caller count", evidence_value(row, ["CALLER_COUNT", "SUPP"])),
            ("Read support", evidence_value(row, ["CALLER_READ_SUPPORT"])),
            ("needLR AF", evidence_value(row, ["NEEDLR_AF"])),
            ("needLR status", evidence_value(row, ["NEEDLR_STATUS"])),
            ("Panel status", evidence_value(row, ["PANEL_STATUS"])),
            ("Phenotype score", evidence_value(row, ["PHENOTYPE_SCORE"])),
            ("Gene-disease evidence", compact_text(evidence_value(row, ["GENE_DISEASE_EVIDENCE_LEVEL", "GENCC"]))),
            ("AnnotSV class", evidence_value(row, ["ANNOTSV_GENERAL_CLASSIFICATION", "ANNotsv_Classification"])),
            ("Pathogenic SV DB", compact_text(evidence_value(row, ["SV_PATHOGENIC_DB_SOURCE"]))),
            ("Benign SV DB / AF", compact_text(evidence_value(row, ["SV_BENIGN_DB_SOURCE", "SV_BENIGN_DB_AFMAX"]))),
            ("LongPhase", evidence_value(row, ["LONGPHASE_PHASED", "LONGPHASE_MATCH"])),
            ("Straglr / TLDR", evidence_value(row, ["STRAGLR_MATCH"]) + " / " + evidence_value(row, ["TLDR_MATCH"])),
            ("Nearby phased SNVs", evidence_value(row, ["WHATSHAP_PHASED_HET_COUNT"])),
            ("Local mean 5mC", evidence_value(row, ["METHYLATION_5MC_MEAN_PERCENT"])),
        ]
        ax_ev.axis("off")
        left = evidence[:9]
        right = evidence[9:]
        for col_x, block in [(0.01, left), (0.52, right)]:
            y = 0.96
            for label, value in block:
                ax_ev.text(col_x, y, f"{label}:", fontweight="bold", fontsize=9.2, transform=ax_ev.transAxes, va="top")
                ax_ev.text(col_x + 0.19, y, str(value), fontsize=9.2, transform=ax_ev.transAxes, va="top", wrap=True)
                y -= 0.105

        formatter = FuncFormatter(lambda x, _: f"{x / 1e6:.3f}")
        ax_sv.xaxis.set_major_formatter(formatter)
        ax_gene.xaxis.set_major_formatter(formatter)
        ax_sv.tick_params(axis="x", labelbottom=False)
        if methyl:
            inset.xaxis.set_major_formatter(formatter)
            ax_gene.tick_params(axis="x", labelbottom=False)
            inset.set_xlabel(f"{chrom} position (Mb)", fontsize=8)
        else:
            ax_gene.set_xlabel(f"{chrom} position (Mb)")

        fig.suptitle(f"{args.title_prefix} {rank}: {gene}", fontsize=16, fontweight="bold", y=0.995)
        fig.text(
            0.5,
            0.008,
            "The locus plot links the master SV to its annotated gene and available contextual evidence. Methylation and nearby phased small variants are context, not independent SV validation.",
            ha="center",
            fontsize=8.5,
        )
        fig.tight_layout(rect=[0, 0.035, 1, 0.965])

        stem = f"{rank:02d}_{safe_name(gene)}_{safe_name(sv_id)}"
        outputs = save_figure(fig, outdir / stem)
        plt.close(fig)

        manifest.append(
            {
                "rank": rank,
                "SV_ID": sv_id,
                "gene": gene,
                "chrom": chrom,
                "start": start,
                "end": end,
                "svtype": svtype,
                "priority_score": row["_priority"],
                "selection_score_field": score_col or ".",
                "needLR_AF": row["_af"] if pd.notna(row["_af"]) else ".",
                "large_gene_centered_view": "YES" if large_gene_centered else "NO",
                "sv_span_bp": sv_span,
                "figure_prefix": str(outdir / stem),
                "methylation_records_plotted": len(methyl),
                "nearby_genes_plotted": len(nearby),
            }
        )

    pd.DataFrame(manifest).to_csv(outdir / "candidate_locus_manifest.tsv", sep="\t", index=False)
    print(f"[OK] candidate_loci={len(manifest)} output={outdir}")


if __name__ == "__main__":
    main()
