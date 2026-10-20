#!/usr/bin/env python3
"""Create one detailed locus figure per prioritized SV-gene pair.

The locus view links a master SV to its overlapping gene, nearby genes, optional
indexed modkit methylation calls and the principal evidence fields carried by
the integrated table.  It is an interpretation aid, not read-level validation.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

from plot_utils import SVTYPE_COLORS, first_existing, normalize_svtype, numeric, read_tsv, save_figure, set_thesis_style
from candidate_review import enrich_events, event_evidence


def parse_args():
    p = argparse.ArgumentParser(description="Generate detailed candidate-specific SV locus figures.")
    p.add_argument("--input", required=True, help="Integrated/multimodal SV-gene TSV")
    p.add_argument("--gene-bed", required=True, help="BED with chrom, start, end, gene label")
    p.add_argument("--gene-gtf", help="Optional matching GRCh38 GTF for exon-resolved gene tracks")
    p.add_argument("--evidence-table", help="Detailed table used to fill missing selected-row evidence")
    p.add_argument("--depth-bins", help="Optional existing large-SV depth-bin TSV")
    p.add_argument("--depth-bed", help="Optional actual mosdepth regions BED for local measured depth")
    p.add_argument("--platform", choices=["lrs","srs"], default="lrs")
    p.add_argument("--selection-table", help="Optional shortlist containing exact SV_ID and GENE keys")
    p.add_argument("--methylation-bed", default=None, help="Optional tabix-indexed modkit bedMethyl")
    p.add_argument("--methylation-units", choices=["percent", "fraction"], default="percent")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--flank", type=int, default=50000)
    p.add_argument("--title-prefix", default="Candidate locus")
    p.add_argument("--samplot-manifest", action="append", default=[], help="Existing Samplot manifest(s); joined by exact master ID")
    p.add_argument("--phase-pairs", help="Existing same-gene SNV/SV pairing TSV")
    p.add_argument("--assembly-manifest", help="Optional TSV with SV_ID and PAF_PATH for existing local assemblies")
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


def methylation_records(path, chrom, start, end, units="percent"):
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
            if f[3].lower() not in {"m", "5mc", "c+m"}:
                continue
            try:
                pos = (int(f[1]) + int(f[2])) / 2
                cov = float(f[9])
                pct = float(f[10])
            except ValueError:
                continue
            if units == "fraction":
                pct *= 100
            if np.isfinite(pct) and 0 <= pct <= 100 and np.isfinite(cov) and cov >= 5:
                rows.append((pos, pct, cov, f[3]))
    except Exception:
        pass
    finally:
        tbx.close()
    return rows


def load_selected_exons(path, selected_genes):
    """Keep true exon intervals for selected genes; do not infer exons from BED4."""
    records = []
    if not path:
        return pd.DataFrame(columns=["chrom", "start", "end", "gene", "transcript"])
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', fields[8]))
            gene = attrs.get("gene_name", attrs.get("gene_id", "."))
            if gene not in selected_genes:
                continue
            transcript = attrs.get("transcript_id", ".")
            if transcript == ".":
                continue
            records.append({"chrom": normalize_chrom(fields[0]), "start": int(fields[3])-1, "end": int(fields[4]), "gene": gene, "transcript": transcript})
    return pd.DataFrame(records, columns=["chrom", "start", "end", "gene", "transcript"]).drop_duplicates()


def representative_exons(exons, gene, chrom, preferred="."):
    sub = exons[exons["gene"].eq(gene) & exons["chrom"].eq(chrom)].copy()
    if sub.empty:
        return sub, ".", "NO_EXON_ANNOTATION"
    available = set(sub["transcript"])
    choices = [t.strip() for t in re.split(r"[;,|]", preferred) if t.strip() in available]
    if choices:
        chosen = sorted(choices)[0]
        basis = "VEP_CANONICAL_TRANSCRIPT_AVAILABLE_IN_GTF"
    else:
        sub["_length"] = sub["end"] - sub["start"]
        totals = sub.groupby("transcript")["_length"].sum()
        chosen = sorted(totals.index, key=lambda t: (-totals[t], t))[0]
        basis = "LARGEST_SUMMED_EXON_LENGTH_FOR_DISPLAY_ONLY"
    return sub[sub["transcript"].eq(chosen)].sort_values("start"), chosen, basis


def evidence_value(row, candidates):
    for col in candidates:
        if col in row.index:
            value = row[col]
            if pd.notna(value) and str(value).strip() not in {"", ".", "NA", "N/A", "nan", "None"}:
                return str(value)
    return "."


def compact_gene_label(value, max_chars=18):
    text = str(value)
    if len(text) <= max_chars:
        return text
    return text[: max_chars - 3] + "..."


def assign_gene_lanes(frame, locus_start, locus_end, max_lanes=6):
    """Pack nearby genes into a small number of non-overlapping visual lanes."""
    if frame.empty:
        return frame.copy()

    work = frame.sort_values(["start", "end", "gene"]).copy()
    lane_ends = []
    lane_values = []
    gap = max(250.0, (float(locus_end) - float(locus_start)) * 0.004)

    for _, row in work.iterrows():
        start = float(row["start"])
        end = float(row["end"])
        lane = None

        for idx, lane_end in enumerate(lane_ends):
            if start >= lane_end + gap:
                lane = idx
                lane_ends[idx] = end
                break

        if lane is None:
            if len(lane_ends) < max_lanes:
                lane = len(lane_ends)
                lane_ends.append(end)
            else:
                lane = int(np.argmin(lane_ends))
                lane_ends[lane] = max(lane_ends[lane], end)

        lane_values.append(lane)

    work["_lane"] = lane_values
    return work


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


def load_samplot_images(paths):
    """Use manifests, never filename guessing, to associate read evidence."""
    images = {}
    for path in paths:
        table = read_tsv(path)
        for row in table.to_dict("records"):
            if row.get("STATUS", "OK") != "OK":
                continue
            image = row.get("IMAGE", row.get("output"))
            if image and Path(str(image)).is_file():
                images.setdefault(str(row["SV_ID"]), []).append(str(image))
    return images


def plot_assembly_alignment(path, prefix):
    """Plot reported PAF alignment blocks from an existing local assembly.

    Separate reference chromosomes and contigs get separate axes. These lines
    connect block endpoints; they are not a base-by-base sequence dot matrix.
    """
    rows = []
    with open(path) as handle:
        for line in handle:
            f = line.rstrip().split("\t")
            if len(f) < 12:
                raise ValueError(f"Invalid PAF record in {path}")
            qlen, qs, qe, tlen, ts, te, mapq = map(int, [f[1], f[2], f[3], f[6], f[7], f[8], f[11]])
            if f[4] not in {"+", "-"} or not 0 <= qs < qe <= qlen or not 0 <= ts < te <= tlen:
                raise ValueError(f"Invalid PAF geometry in {path}")
            rows.append({"CONTIG":f[0],"CONTIG_LENGTH":qlen,"QUERY_START0":qs,"QUERY_END0":qe,
                         "STRAND":f[4],"REF_CHROM":f[5],"REF_START0":ts,"REF_END0":te,"MAPQ":mapq})
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError(f"No alignments in supplied PAF: {path}")
    groups = list(frame.groupby(["CONTIG", "REF_CHROM"], sort=False))
    fig, axes = plt.subplots(len(groups),1,figsize=(10,max(4,len(groups)*3)),squeeze=False)
    for ax, ((contig,chrom), group) in zip(axes.flat,groups):
        for r in group.itertuples():
            q = [r.QUERY_START0,r.QUERY_END0] if r.STRAND == "+" else [r.QUERY_END0,r.QUERY_START0]
            ax.plot([r.REF_START0,r.REF_END0],q,color="#0072B2" if r.STRAND == "+" else "#CC79A7",lw=2,alpha=.75)
        ax.set_title(f"{contig} aligned to {chrom}",loc="left",fontsize=11)
        ax.set_xlabel("Reference position (bp; 0-based)")
        ax.set_ylabel("Assembly contig position (bp)")
        ax.ticklabel_format(style="plain",useOffset=False)
    fig.suptitle("Local assembly alignment blocks",fontweight="bold")
    fig.text(.5,.01,"Blue: forward alignment; purple: reverse alignment. Lines join PAF block endpoints; no base-level identity is shown.\nA discordant block pattern is a review clue, not a validated reconstruction. Consult alignment quality and the assembly itself.",ha="center",fontsize=9)
    fig.subplots_adjust(left=.15,right=.97,top=.9,bottom=.17,hspace=.65)
    save_figure(fig,prefix)
    plt.close(fig)
    frame.to_csv(str(prefix)+"_source.tsv",sep="\t",index=False)
    return str(prefix)+".png"


def write_candidate_card(locus_image, samplots, phase_image, assembly_image, row, prefix):
    """One printable page; the full-resolution source figures remain available."""
    fig = plt.figure(figsize=(16.5,23.4))  # A2 portrait: one page with readable source panels.
    grid = fig.add_gridspec(3,2,height_ratios=[1.6,1,.8])
    axes = [fig.add_subplot(grid[0,:]),fig.add_subplot(grid[1,:]),fig.add_subplot(grid[2,0]),fig.add_subplot(grid[2,1])]
    panels = [("Locus, depth and methylation (when measured)", [locus_image]),
              ("Reads: existing Samplot", samplots[:2]),
              ("Gene phase: verified pair evidence", [phase_image] if phase_image else []),
              ("Existing local assembly alignment", [assembly_image] if assembly_image else [])]
    sources = []
    for ax, (label, paths) in zip(axes,panels):
        ax.axis("off")
        ax.set_title(label,loc="left",fontsize=11)
        if not paths:
            ax.text(.5,.5,"Evidence unavailable",transform=ax.transAxes,ha="center",fontsize=13,color="#777777")
        elif len(paths) == 1:
            ax.imshow(plt.imread(paths[0]))
        else:
            for i, path in enumerate(paths):
                inset = ax.inset_axes([0, .51 if i == 0 else 0, 1, .48])
                inset.imshow(plt.imread(path))
                inset.axis("off")
        for path in paths:
            sources.append({"SV_ID":row["SV_ID"],"GENE":row.get("GENE",row.get("GENES",".")),"PANEL":label,"SOURCE_IMAGE":path})
    fig.suptitle(f"Candidate review · {row.get('GENE',row.get('GENES','.'))} · {row['SV_ID']}",fontweight="bold",fontsize=16)
    fig.text(.5,.015,"A review summary of available evidence; missing panels are not negative findings. Open the separate original figures for full detail.\nCandidate selection, computational support and phase are not proof of disease causality.",ha="center",fontsize=10)
    fig.subplots_adjust(left=.025,right=.98,top=.93,bottom=.07,hspace=.15,wspace=.08)
    save_figure(fig,prefix,dpi=200)
    plt.close(fig)
    pd.DataFrame(sources).to_csv(str(prefix)+"_sources.tsv",sep="\t",index=False)


def main():
    args = parse_args()
    set_thesis_style()
    df = read_tsv(args.input)
    genes = load_gene_bed(args.gene_bed)
    samplot_images = load_samplot_images(args.samplot_manifest)
    assemblies = {}
    if args.assembly_manifest:
        assembly_table = read_tsv(args.assembly_manifest, required=["SV_ID", "PAF_PATH"])
        if assembly_table["SV_ID"].duplicated().any():
            raise ValueError("Assembly manifest must have one PAF per master SV ID")
        assemblies = dict(zip(assembly_table["SV_ID"].astype(str), assembly_table["PAF_PATH"].astype(str)))

    id_col = first_existing(df, ["SV_ID", "ID"])
    gene_col = first_existing(df, ["GENES", "ANNotsv_Gene", "Gene", "GENE"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    end_col = first_existing(df, ["END"])
    chr2_col = first_existing(df, ["CHR2"])
    pos2_col = first_existing(df, ["POS2"])
    type_col = first_existing(df, ["SVTYPE", "SV_type"])
    score_col = first_existing(
        df,
        [
            "FINAL_GENE_RELEVANCE_DISPLAY_SCORE",
            "GENE_RELEVANCE_DISPLAY_SCORE",
            "EVENT_GENE_RELEVANCE_SCORE",
            "INTEGRATED_DISCOVERY_SCORE",
            "integrated_discovery_score",
            "PHENOTYPE_SCORE",
            "ALLELE_RESEARCH_SCORE",
        ],
    )
    rank_col = first_existing(
        df,
        [
            "DIAGNOSTIC_REVIEW_RANK",
            "FINAL_EVENT_RANK_GLOBAL",
            "EVENT_RANK_GLOBAL",
            "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS",
            "EVENT_RANK_WITHIN_PANEL_STATUS",
        ],
    )
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
    work["_final_rank"] = (
        numeric(work[rank_col])
        if rank_col
        else pd.Series(np.nan, index=work.index)
    )
    work["_rank_missing"] = work["_final_rank"].isna().astype(int)
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
            [
                "_rank_missing",
                "_final_rank",
                "_priority",
                "_population_rank",
                "_caller_count",
                id_col,
                "_gene",
            ],
            ascending=[True, True, False, False, False, True, True],
            na_position="last",
        )
        .drop_duplicates([id_col, gene_col], keep="first")
        .copy()
    )

    if args.selection_table:
        shortlist = read_tsv(args.selection_table, required=["SV_ID", "GENE"])
        chosen_keys = set(zip(shortlist["SV_ID"].astype(str), shortlist["GENE"].astype(str)))
        work = work.loc[[(str(r[id_col]), str(r["_gene"])) in chosen_keys for _, r in work.iterrows()]].copy()
    bucket_col = first_existing(work, ["EVENT_REVIEW_BUCKET"]) if not args.selection_table else None
    bucket_rank_col = first_existing(work, ["EVENT_RANK_WITHIN_BUCKET"])
    bucket_order = [
        "BREAKPOINT_GENE_CANDIDATE",
        "LARGE_CNV_GENE_CANDIDATE",
        "VERY_LARGE_CNV_GENE_CONTEXT",
        "INSERTION_GENE_CANDIDATE",
        "SMALL_MEDIUM_CNV_GENE_CANDIDATE",
        "LARGE_COMPLEX_INTERVAL_CONTEXT",
    ]

    if bucket_col:
        work["_bucket_rank"] = (
            numeric(work[bucket_rank_col]).fillna(np.inf)
            if bucket_rank_col
            else np.inf
        )
        per_bucket = max(1, args.top_n // len(bucket_order))
        selected_indices = []

        for bucket in bucket_order:
            sub = work[work[bucket_col].eq(bucket)].sort_values(
                [
                    "_rank_missing",
                    "_final_rank",
                    "_bucket_rank",
                    "_priority",
                    "_population_rank",
                    "_caller_count",
                ],
                ascending=[True, True, True, False, False, False],
                na_position="last",
            )

            if "PANEL_STATUS" in sub.columns and per_bucket >= 2:
                panel_sub = sub[
                    sub["PANEL_STATUS"].fillna("").astype(str).eq("PANEL_GENE")
                ]
                nonpanel_sub = sub[
                    ~sub.index.isin(panel_sub.index)
                ]

                chosen = []
                if not panel_sub.empty:
                    chosen.extend(panel_sub.head(1).index.tolist())
                if not nonpanel_sub.empty:
                    chosen.extend(
                        nonpanel_sub.head(
                            per_bucket - len(chosen)
                        ).index.tolist()
                    )

                if len(chosen) < per_bucket:
                    chosen.extend(
                        sub.loc[~sub.index.isin(chosen)]
                        .head(per_bucket - len(chosen))
                        .index.tolist()
                    )

                selected_indices.extend(chosen)
            else:
                selected_indices.extend(
                    sub.head(per_bucket).index.tolist()
                )

        selected_indices = list(dict.fromkeys(selected_indices))

        if len(selected_indices) < args.top_n:
            remainder = work.loc[~work.index.isin(selected_indices)].sort_values(
                [
                    "_rank_missing",
                    "_final_rank",
                    "_priority",
                    "_population_rank",
                    "_caller_count",
                ],
                ascending=[True, True, False, False, False],
                na_position="last",
            )
            selected_indices.extend(
                remainder.head(args.top_n - len(selected_indices)).index.tolist()
            )

        work = work.loc[selected_indices].head(args.top_n).copy()
    else:
        work = work.head(args.top_n).copy()

    outdir = Path(args.out_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    work["SV_ID"] = work[id_col]
    work["GENE"] = work["_gene"]
    work = enrich_events(work, args.evidence_table)
    phase_images = {}
    if args.phase_pairs:
        from plot_phasing_qc import draw_gene_phase
        pairs = read_tsv(args.phase_pairs)
        pairs = pairs[pairs["GENE"].astype(str).isin(set(work["GENE"].astype(str)))]
        phase_dir = outdir / "gene_phase"
        draw_gene_phase(pairs, work, phase_dir, max(1,args.top_n))
        phase_manifest = read_tsv(phase_dir / "phase_gene_manifest.tsv")
        phase_images = dict(zip(phase_manifest["GENE"].astype(str),phase_manifest["FIGURE_PREFIX"].astype(str)+".png"))
    exons = load_selected_exons(args.gene_gtf, set(work["_gene"]))
    exons.to_csv(outdir / "candidate_exon_annotations.tsv", sep="\t", index=False)
    depth = read_tsv(args.depth_bins) if args.depth_bins else pd.DataFrame()
    measured_depth = pd.DataFrame()
    if args.depth_bed:
        measured_depth = pd.read_csv(args.depth_bed,sep="\t",header=None,comment="#")
        measured_depth = pd.DataFrame({"CHROM":measured_depth.iloc[:,0].map(normalize_chrom),
                                      "START0":pd.to_numeric(measured_depth.iloc[:,1],errors="raise"),
                                      "END0":pd.to_numeric(measured_depth.iloc[:,2],errors="raise"),
                                      "MEAN_DEPTH":pd.to_numeric(measured_depth.iloc[:,-1],errors="raise")})
        if measured_depth["MEAN_DEPTH"].lt(0).any() or measured_depth["END0"].le(measured_depth["START0"]).any():
            raise ValueError("Invalid measured mosdepth bins")
    manifest = []

    for rank, (_, row) in enumerate(work.iterrows(), 1):
        sv_id = str(row[id_col])
        gene = str(row["_gene"])
        svtype = str(row["_svtype"])

        chrom = normalize_chrom(row[chrom_col])
        start = int(row["_start"])
        end = int(row["_end"])

        target_gene_all = genes[
            genes["gene"].astype(str).eq(gene)
        ].copy()

        if svtype == "BND" and not target_gene_all.empty:
            target_chrom = str(target_gene_all.iloc[0]["chrom"])
            chrom2 = (
                normalize_chrom(row[chr2_col])
                if chr2_col and pd.notna(row.get(chr2_col))
                else None
            )
            pos2 = (
                numeric(pd.Series([row[pos2_col]])).iloc[0]
                if pos2_col
                else np.nan
            )

            if target_chrom == chrom:
                breakpoint_pos = start
            elif chrom2 and target_chrom == chrom2 and pd.notna(pos2):
                breakpoint_pos = int(pos2)
            else:
                gene_mid = float(
                    (
                        target_gene_all["start"].min()
                        + target_gene_all["end"].max()
                    )
                    / 2
                )
                candidates = [(chrom, start)]
                if chrom2 and pd.notna(pos2):
                    candidates.append((chrom2, int(pos2)))
                same_chrom = [
                    (c, p)
                    for c, p in candidates
                    if c == target_chrom
                ]
                breakpoint_pos = (
                    min(same_chrom, key=lambda x: abs(x[1] - gene_mid))[1]
                    if same_chrom
                    else start
                )

            chrom = target_chrom
            start = int(breakpoint_pos)
            end = start

        elif end < start:
            start, end = end, start

        is_breakpoint = svtype in {"INS", "BND"} or end == start

        event_span_col = first_existing(work, ["SV_SPAN_BP", "SV_EVENT_SPAN_BP"])
        event_span = (
            numeric(pd.Series([row[event_span_col]])).iloc[0]
            if event_span_col
            else np.nan
        )
        sv_span = (
            float(event_span)
            if pd.notna(event_span)
            else abs(end - start)
        )
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
                [nearby[keep_target], nearby[~keep_target].nsmallest(15, "_distance")]
            ).drop_duplicates().sort_values("start")

        methyl = methylation_records(args.methylation_bed, chrom, locus_start, locus_end, args.methylation_units)
        shown_exons, transcript, transcript_basis = representative_exons(exons, gene, chrom, evidence_value(row, ["VEP_CANONICAL_TRANSCRIPTS"]))
        depth_here = pd.DataFrame()
        raw_depth_track = False
        if not depth.empty and {"SV_ID", "CHROM", "PLOT_START", "PLOT_END", "NORMALIZED_DEPTH"}.issubset(depth.columns):
            depth_here = depth[depth["SV_ID"].astype(str).eq(sv_id) & depth["CHROM"].map(normalize_chrom).eq(chrom)].copy()
            depth_here["_pos"] = (numeric(depth_here["PLOT_START"]) + numeric(depth_here["PLOT_END"])) / 2
            depth_here["_depth"] = numeric(depth_here["NORMALIZED_DEPTH"])
            depth_here = depth_here[depth_here["_pos"].between(locus_start, locus_end)].dropna(subset=["_pos", "_depth"]).sort_values("_pos")
        if depth_here.empty and not measured_depth.empty and svtype in {"DEL","DUP"}:
            depth_here = measured_depth[measured_depth["CHROM"].eq(chrom) & measured_depth["END0"].gt(locus_start) & measured_depth["START0"].lt(locus_end)].copy()
            depth_here["_pos"] = (depth_here["START0"]+depth_here["END0"])/2 + 1
            depth_here["_depth"] = depth_here["MEAN_DEPTH"]
            raw_depth_track = True
        heights = [1.0, 1.55] + ([0.85] if not depth_here.empty else []) + ([0.85] if methyl else []) + [2.2]
        fig, axes = plt.subplots(len(heights), 1, figsize=(16.5, 10.8 + (len(heights)-3)*1.8), gridspec_kw={"height_ratios": heights})
        ax_sv, ax_gene, ax_ev = axes[0], axes[1], axes[-1]
        track_index = 2
        ax_depth = axes[track_index] if not depth_here.empty else None
        if ax_depth is not None:
            track_index += 1
        ax_methyl = axes[track_index] if methyl else None

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
            nearby = assign_gene_lanes(
                nearby,
                locus_start,
                locus_end,
                max_lanes=6,
            )
            for _, g in nearby.iterrows():
                y = float(g["_lane"])
                target = str(g["gene"]) == gene
                face = "#0072B2" if target else "#C7CDD3"
                gstart = float(g["start"])
                gend = float(g["end"])
                if target and not shown_exons.empty:
                    ax_gene.plot([gstart, gend], [y, y], color=face, linewidth=1.2)
                    for _, exon in shown_exons.iterrows():
                        ax_gene.add_patch(Rectangle((exon["start"], y-0.22), max(exon["end"]-exon["start"], 1), 0.44, facecolor=face, edgecolor="none"))
                else:
                    ax_gene.add_patch(Rectangle((gstart, y-0.22), max(gend-gstart, 1), 0.44, facecolor=face, edgecolor="none"))
                ax_gene.text(
                    (gstart + gend) / 2,
                    y + 0.27,
                    compact_gene_label(g["gene"]),
                    fontsize=8.0 if target else 7.2,
                    fontweight="bold" if target else "normal",
                    ha="center",
                    va="bottom",
                    clip_on=True,
                )
            max_lane = int(nearby["_lane"].max())
            ax_gene.set_ylim(-0.55, max_lane + 0.65)
        else:
            ax_gene.text(0.5, 0.5, "No gene BED entries in displayed interval", transform=ax_gene.transAxes, ha="center")
            ax_gene.set_ylim(0, 1)
        ax_gene.set_xlim(locus_start, locus_end)
        ax_gene.set_yticks([])
        ax_gene.set_ylabel("Genes")
        ax_gene.grid(axis="x", color="#E6E6E6", linewidth=0.7)
        exon_note = f"Target exons: {transcript}; representative transcript for display" if not shown_exons.empty else "Gene intervals only; exon annotation unavailable"
        ax_gene.set_title(exon_note, loc="left", fontsize=9)

        if ax_depth is not None:
            ax_depth.plot(depth_here["_pos"], depth_here["_depth"], color="#0072B2", linewidth=1)
            if not raw_depth_track:
                ax_depth.axhline(1, color="#777777", linestyle="--", linewidth=0.8)
            ax_depth.set_xlim(locus_start, locus_end)
            ax_depth.set_ylabel("Mean depth (reads)" if raw_depth_track else "Depth / flanks", fontsize=9)
            ax_depth.grid(axis="y", color="#EEEEEE", linewidth=0.6)
            ax_depth.set_title("Existing depth bins from the same BAM; no copy-number state inferred", loc="left", fontsize=9)

        if methyl:
            pos = [x[0] for x in methyl]
            pct = [x[1] for x in methyl]
            cov = [x[2] for x in methyl]
            ax_methyl.scatter(
                pos,
                pct,
                s=np.clip(np.sqrt(cov) * 4, 6, 36),
                alpha=0.65,
            )
            ax_methyl.set_xlim(locus_start, locus_end)
            ax_methyl.set_ylim(0, 100)
            ax_methyl.set_ylabel("5mC %", fontsize=8.5)
            ax_methyl.tick_params(labelsize=8)
            ax_methyl.grid(axis="y", color="#EEEEEE", linewidth=0.6)

        observations = event_evidence(row,platform=args.platform)
        evidence = [
            ("Overlapping gene", gene),
            ("Final event rank", (
                f"{int(row['_final_rank'])}"
                if pd.notna(row["_final_rank"])
                else "."
            )),
            ("Gene relevance tier", evidence_value(row, ["FINAL_GENE_RELEVANCE_TIER", "GENE_RELEVANCE_TIER"])),
            ("Event bucket", compact_text(evidence_value(row, ["EVENT_REVIEW_BUCKET"]), 52, 2)),
            ("Gene relationship", compact_text(evidence_value(row, ["SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP"]), 52, 2)),
            ("SV gene count", evidence_value(row, ["GENES_AFFECTED", "SV_GENE_COUNT"])),
            ("Callers", evidence_value(row, ["CALLERS"])),
            ("Caller count", evidence_value(row, ["CALLER_COUNT", "SUPP"])),
            ("Read support", evidence_value(row, ["READ_SUPPORT", "CALLER_READ_SUPPORT", "EVENT_MAX_CALLER_READ_SUPPORT"])),
            ("needLR AF", evidence_value(row, ["NEEDLR_AF"])),
            ("Population state", observations["Population"][1].replace("\n", " ")),
            ("gnomAD event AF", evidence_value(row, ["GNOMAD_SV_AF"])),
            ("Panel status", evidence_value(row, ["PANEL_STATUS"])),
            ("Phenotype scope", compact_text(evidence_value(row, ["FINAL_PHENOTYPE_RANKING_SCOPE", "PHENOTYPE_SCORE_SCOPE"]), 48, 2)),
            ("Mechanism", compact_text(evidence_value(row, ["FINAL_INHERITANCE_MECHANISM_CLASS", "INHERITANCE_MECHANISM_CLASS"]), 48, 2)),
            ("Depth evidence", observations["CNVpytor depth" if args.platform == "srs" else "Local depth"][1].replace("\n", " ")),
            ("Transcript shown", compact_text(transcript, 48, 2)),
            ("Exon display basis", "VEP canonical" if transcript_basis.startswith("VEP") else "Largest exon length" if not shown_exons.empty else "Unavailable"),
            ("GRIDSS" if args.platform == "srs" else "LongPhase", evidence_value(row, ["GRIDSS_SUPPORT"] if args.platform == "srs" else ["LONGPHASE_PHASED", "LONGPHASE", "LONGPHASE_MATCH"])),
            ("Repeat locus / MELT" if args.platform == "srs" else "Straglr / TLDR", evidence_value(row, ["EXPANSIONHUNTER_LOCUS_OVERLAP"] if args.platform == "srs" else ["STRAGLR", "STRAGLR_MATCH"]) + " / " + evidence_value(row, ["MELT_SUPPORT"] if args.platform == "srs" else ["TLDR", "TLDR_MATCH"])),
            ("SV phase" if args.platform == "srs" else "SNV + SV pairing", "Unavailable: small-variant phase only" if args.platform == "srs" else compact_text(evidence_value(row, ["RECESSIVE_PAIR_STATUS"]), 48, 2)),
            ("Methylation context", "Unavailable from SRS" if args.platform == "srs" else observations["Methylation"][1]),
        ]
        ax_ev.axis("off")
        left = evidence[:11]
        right = evidence[11:]
        for col_x, block in [(0.01, left), (0.52, right)]:
            y = 0.96
            for label, value in block:
                ax_ev.text(col_x, y, f"{label}:", fontweight="bold", fontsize=8.7, transform=ax_ev.transAxes, va="top")
                ax_ev.text(col_x + 0.19, y, compact_text(value, 48, 2), fontsize=8.2, transform=ax_ev.transAxes, va="top", wrap=True)
                y -= 0.088

        use_bp = locus_end-locus_start < 10000
        formatter = FuncFormatter(lambda x, _: f"{x:,.0f}" if use_bp else f"{x / 1e6:.3f}")
        ax_sv.xaxis.set_major_formatter(formatter)
        ax_gene.xaxis.set_major_formatter(formatter)
        if ax_depth is not None:
            ax_depth.xaxis.set_major_formatter(formatter)
            ax_gene.tick_params(axis="x", labelbottom=False)
        ax_sv.tick_params(axis="x", labelbottom=False)
        if methyl:
            ax_methyl.xaxis.set_major_formatter(formatter)
            ax_gene.tick_params(axis="x", labelbottom=False)
            if ax_depth is not None:
                ax_depth.tick_params(axis="x", labelbottom=False)
            ax_methyl.set_xlabel(f"{chrom} position ({'bp' if use_bp else 'Mb'})", fontsize=8.5)
        else:
            (ax_depth if ax_depth is not None else ax_gene).set_xlabel(f"{chrom} position ({'bp' if use_bp else 'Mb'})")

        fig.suptitle(f"{args.title_prefix} {rank}: {gene}", fontsize=16, fontweight="bold", y=0.995)
        fig.text(
            0.5,
            0.008,
            "The locus plot links the master SV to its annotated gene and available contextual evidence. Methylation and nearby phased small variants are context, not independent SV validation.",
            ha="center",
            fontsize=8.5,
        )
        fig.subplots_adjust(
            left=0.075,
            right=0.985,
            top=0.925,
            bottom=0.085,
            hspace=0.55,
        )

        stem = f"{rank:02d}_{safe_name(gene)}_{safe_name(sv_id)}"
        outputs = save_figure(fig, outdir / stem, dpi=300)
        depth_here.to_csv(outdir / (stem+"_depth_source.tsv"),sep="\t",index=False)
        plt.close(fig)

        assembly_image = None
        if sv_id in assemblies:
            assembly_image = plot_assembly_alignment(assemblies[sv_id], outdir / (stem+"_assembly"))
        card_prefix = outdir / "candidate_cards" / stem
        card_prefix.parent.mkdir(exist_ok=True)
        locus_png = next(
            path for path in outputs if path.suffix.lower() == ".png"
        )
        write_candidate_card(
            str(locus_png),
            samplot_images.get(sv_id, []),
            phase_images.get(gene),
            assembly_image,
            row,
            card_prefix,
        )

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
                "selection_rank": (
                    int(row["_final_rank"])
                    if pd.notna(row["_final_rank"])
                    else "."
                ),
                "selection_rank_field": rank_col or ".",
                "selection_score_field": score_col or ".",
                "needLR_AF": row["_af"] if pd.notna(row["_af"]) else ".",
                "large_gene_centered_view": "YES" if large_gene_centered else "NO",
                "sv_span_bp": sv_span,
                "figure_prefix": str(outdir / stem),
                "candidate_card_prefix": str(card_prefix),
                "methylation_records_plotted": len(methyl),
                "depth_bins_plotted": len(depth_here),
                "transcript_shown": transcript,
                "exon_count": len(shown_exons),
                "transcript_selection_basis": transcript_basis,
                "exon_annotation_scope": "GTF_EXONS" if not shown_exons.empty else "GENE_INTERVAL_ONLY",
                "evidence_join_status": row.get("EVIDENCE_JOIN_STATUS", "."),
                "nearby_genes_plotted": len(nearby),
            }
        )

    pd.DataFrame(manifest).to_csv(outdir / "candidate_locus_manifest.tsv", sep="\t", index=False)
    print(f"[OK] candidate_loci={len(manifest)} output={outdir}")


if __name__ == "__main__":
    main()

