#!/usr/bin/env python3
"""Show SV positions, spacing and optional repeat/depth context.

One row is one chromosome, with a common position scale. This is a view of
reported events, not a test of pathogenicity or a reconstruction of chromosomes.
Input SV positions are VCF positions (1-based); BED tracks are 0-based, half-open.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path as DrawPath
import numpy as np
import pandas as pd

from candidate_review import MISSING, select_events
from plot_utils import SVTYPE_COLORS, chromosome_sort_key, numeric, read_tsv, save_figure, set_thesis_style


def chrom_name(value):
    if value is None or str(value).upper() in MISSING:
        return "."
    value = str(value).removeprefix("chr").upper()
    return "chr" + ("M" if value == "MT" else value)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_lengths(path):
    """Use the index of the actual reference, rather than hard-coded lengths."""
    lengths = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            chrom, length = chrom_name(fields[0]), int(fields[1])
            if re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)", chrom):
                if chrom in lengths or length < 1:
                    raise ValueError(f"Duplicate chromosome alias or invalid length: {chrom}")
                lengths[chrom] = length
    if not lengths:
        raise ValueError("Reference index contains no canonical chromosomes")
    return dict(sorted(lengths.items(), key=lambda item: chromosome_sort_key(item[0])))


def read_bed(path, lengths):
    """Validate BED coordinates; never silently clip a resource from another build."""
    rows = []
    if path:
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.split()
                chrom, start, end = chrom_name(fields[0]), int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError(f"Invalid BED interval in {path}: {line.strip()}")
                if chrom not in lengths:
                    continue
                if end > lengths[chrom]:
                    raise ValueError(f"BED extends beyond reference chromosome: {path}: {line.strip()}")
                rows.append((chrom, start, end, fields[3] if len(fields) > 3 else "."))
    return pd.DataFrame(rows, columns=["CHROM", "START0", "END0", "LABEL"])


def interval_index(bed):
    """Merge overlapping intervals so breakpoint queries are quick and unambiguous."""
    index = {}
    for chrom, group in bed.groupby("CHROM"):
        merged = []
        for start, end in group.sort_values("START0")[["START0", "END0"]].itertuples(index=False, name=None):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        index[chrom] = np.asarray(merged, dtype=int)
    return index


def overlaps(index, chrom, start, end):
    intervals = index.get(chrom)
    if intervals is None or not len(intervals):
        return False
    i = np.searchsorted(intervals[:, 0], end, side="left") - 1
    return i >= 0 and intervals[i, 1] > start


def prepare_events(raw, lengths):
    """Count master IDs once; reject conflicting propagated coordinates."""
    required = ["SV_ID", "CHROM", "START", "SVTYPE"]
    if not set(required).issubset(raw.columns):
        raise ValueError(f"Genome map requires {required}")
    work = raw.copy()
    if work["SV_ID"].fillna(".").astype(str).str.upper().isin(MISSING).any():
        raise ValueError("Missing master SV IDs: positions cannot be deduplicated reliably")
    for col in ["END", "CHR2", "POS2"]:
        if col not in work:
            work[col] = "."
    work["CHROM"] = work["CHROM"].map(chrom_name)
    work["CHR2"] = work["CHR2"].map(chrom_name)
    work["SVTYPE"] = work["SVTYPE"].astype(str).str.upper().replace({"TRA": "BND", "CTX": "BND"})
    for col in ["START", "END", "POS2"]:
        work[col] = numeric(work[col])
    geometry = ["CHROM", "START", "END", "SVTYPE", "CHR2", "POS2"]
    conflicts = work.groupby("SV_ID")[geometry].nunique(dropna=False).gt(1).any(axis=1)
    bad = set(conflicts[conflicts].index)
    work = work.drop_duplicates("SV_ID").copy()
    work["MAP_STATUS"] = "PLOTTED"
    work.loc[work["SV_ID"].isin(bad), "MAP_STATUS"] = "CONFLICTING_COORDINATES"
    for i, row in work.iterrows():
        if row["SV_ID"] in bad:
            continue
        chrom, start = row["CHROM"], row["START"]
        if chrom not in lengths:
            work.at[i, "MAP_STATUS"] = "NONCANONICAL_OR_NOT_IN_REFERENCE"
        elif not np.isfinite(start) or start != int(start) or not 1 <= start <= lengths[chrom]:
            work.at[i, "MAP_STATUS"] = "INVALID_START"
        elif row["SVTYPE"] in {"DEL", "DUP", "INV", "CNV"} and (
            not np.isfinite(row["END"]) or row["END"] != int(row["END"]) or not start <= row["END"] <= lengths[chrom]
        ):
            work.at[i, "MAP_STATUS"] = "INVALID_INTERVAL"
    work["PARTNER_STATUS"] = "NOT_APPLICABLE"
    for i, row in work[work["SVTYPE"].eq("BND")].iterrows():
        pos2 = row["POS2"]
        valid = row["CHR2"] in lengths and np.isfinite(pos2) and pos2 == int(pos2) and 1 <= pos2 <= lengths[row["CHR2"]]
        work.at[i, "PARTNER_STATUS"] = "RESOLVED" if valid else "UNRESOLVED_OR_OFF_MAP"
    return work


def endpoints(row):
    points = [(row["CHROM"], int(row["START"]) - 1)]
    if row["SVTYPE"] == "BND" and row["PARTNER_STATUS"] == "RESOLVED":
        points.append((row["CHR2"], int(row["POS2"]) - 1))
    elif row["SVTYPE"] in {"DEL", "DUP", "INV", "CNV"}:
        points.append((row["CHROM"], int(row["END"]) - 1))
    return points


def rainfall_table(events):
    """Use one START anchor per SV; reset distances at each chromosome.

    Do not include the other end of the same deletion: that would confuse its
    length with the distance to the next event. Ties are real zero distances.
    """
    frame = events[events["MAP_STATUS"].eq("PLOTTED")][["SV_ID", "CHROM", "START", "SVTYPE"]].copy()
    frame = frame.sort_values(["CHROM", "START", "SV_ID"])
    frame["PREVIOUS_SV_ID"] = frame.groupby("CHROM")["SV_ID"].shift()
    frame["DISTANCE_BP"] = frame.groupby("CHROM")["START"].diff()
    return frame


def load_depth(path, lengths, exclude_index, display_bin=1000000):
    """Pool actual mosdepth means by length; normalise to a sample baseline.

    This ribbon is relative read depth, not GC-corrected or absolute copy number.
    X/Y use their own baseline; autosomes share an autosomal baseline.
    """
    if not path:
        return pd.DataFrame(columns=["CHROM", "START0", "END0", "MEAN_DEPTH", "RELATIVE_DEPTH"])
    pooled = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            chrom, start, end = chrom_name(fields[0]), int(fields[1]), int(fields[2])
            value = float(fields[-1])  # mosdepth adds a BED label before the depth when provided.
            if chrom not in lengths:
                continue
            if not np.isfinite(value) or value < 0 or start < 0 or end <= start or end > lengths[chrom]:
                raise ValueError(f"Invalid depth bin in {path}: {line.strip()}")
            # Mask reference/exclusion regions so zero depth in gaps cannot lower the baseline.
            if overlaps(exclude_index, chrom, start, end):
                continue
            for b in range(start // display_bin, (end - 1) // display_bin + 1):
                left, right = max(start, b * display_bin), min(end, (b + 1) * display_bin)
                total, width = pooled.get((chrom, b), (0., 0))
                pooled[(chrom, b)] = (total + value * (right - left), width + right - left)
    rows = [{"CHROM": chrom, "START0": b * display_bin, "END0": min((b + 1) * display_bin, lengths[chrom]),
             "MEAN_DEPTH": total / width, "COVERED_BP": width} for (chrom, b), (total, width) in pooled.items()]
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("Supplied depth file has no valid, unmasked canonical bins")
    autosomal = frame[frame["CHROM"].str.fullmatch(r"chr\d+")]
    auto_baseline = np.average(autosomal["MEAN_DEPTH"], weights=autosomal["COVERED_BP"]) if not autosomal.empty else np.nan
    frame["BASELINE_DEPTH"] = auto_baseline
    frame["BASELINE_SCOPE"] = "AUTOSOMAL_MEAN"
    for chrom in ["chrX", "chrY", "chrM"]:
        mask = frame["CHROM"].eq(chrom)
        if mask.any():
            frame.loc[mask, "BASELINE_DEPTH"] = np.average(frame.loc[mask, "MEAN_DEPTH"], weights=frame.loc[mask, "COVERED_BP"])
            frame.loc[mask, "BASELINE_SCOPE"] = "OWN_CHROMOSOME_MEAN"
    frame["RELATIVE_DEPTH"] = frame["MEAN_DEPTH"] / frame["BASELINE_DEPTH"].where(frame["BASELINE_DEPTH"] > 0)
    return frame


def draw_map(events, lengths, tracks, depth, labels, prefix):
    fig, ax = plt.subplots(figsize=(17, max(8, len(lengths) * .68 + 3)))
    y_of = {chrom: len(lengths) - i for i, chrom in enumerate(lengths)}
    maximum = max(lengths.values()) / 1e6
    for chrom, length in lengths.items():
        y = y_of[chrom]
        ax.plot([0, length / 1e6], [y, y], color="#E4E4E4", lw=15, solid_capstyle="butt", zorder=0)
    offsets = {kind: value for kind, value in zip(SVTYPE_COLORS, np.linspace(-.18, .18, len(SVTYPE_COLORS)))}
    for kind, group in events[events["MAP_STATUS"].eq("PLOTTED")].groupby("SVTYPE"):
        color = SVTYPE_COLORS.get(kind, "#999999")
        segments = []
        for _, row in group.iterrows():
            x, y = row["START"] / 1e6, y_of[row["CHROM"]] + offsets.get(kind, 0)
            end = row["END"] / 1e6 if kind in {"DEL", "DUP", "INV", "CNV"} else x
            segments.append([(x, y), (end, y)])
            if kind == "BND" and row["PARTNER_STATUS"] == "RESOLVED":
                x2, y2 = row["POS2"] / 1e6, y_of[row["CHR2"]] + offsets["BND"]
                bend = min(maximum * 1.02, max(x, x2) + maximum*.04)
                curve = DrawPath([(x, y), (bend, y + .32), (bend, y2 + .32), (x2, y2)],
                                 [DrawPath.MOVETO, DrawPath.CURVE4, DrawPath.CURVE4, DrawPath.CURVE4])
                ax.add_patch(PathPatch(curve, fill=False, color=color, alpha=.16, lw=.55, rasterized=True))
        ax.add_collection(LineCollection(segments, colors=color, linewidths=1.6, alpha=.75, rasterized=True))
        ax.scatter(group["START"] / 1e6, group["CHROM"].map(y_of) + offsets.get(kind, 0), s=5, color=color, alpha=.65, rasterized=True)
    track_colors = ["#494949", "#BDBDBD", "#C29162", "#7E57A2"]
    for (name, bed), color, offset in zip(tracks.items(), track_colors, [-.29, -.35, -.41, -.47]):
        segments = [[(r.START0 / 1e6, y_of[r.CHROM] + offset), (r.END0 / 1e6, y_of[r.CHROM] + offset)] for r in bed.itertuples()]
        ax.add_collection(LineCollection(segments, colors=color, linewidths=2, rasterized=True))
    if not depth.empty:
        cmap, norm = plt.get_cmap("RdBu_r"), Normalize(0, 2, clip=True)
        for chrom, group in depth.groupby("CHROM"):
            segs = [[(r.START0 / 1e6, y_of[chrom] - .57), (r.END0 / 1e6, y_of[chrom] - .57)] for r in group.itertuples()]
            lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=4, rasterized=True)
            lc.set_array(group["RELATIVE_DEPTH"].to_numpy())
            ax.add_collection(lc)
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=.018, pad=.03)
        cb.set_label("Relative depth (1 = baseline; colour capped at 2)")
    for n, (_, row) in enumerate(labels.iterrows(), 1):
        hit = events[events["SV_ID"].eq(row["SV_ID"]) & events["MAP_STATUS"].eq("PLOTTED")]
        if hit.empty:
            continue
        r = hit.iloc[0]
        # Labels live in the right margin, avoiding thousands of overlapping names.
        label_y = len(lengths) + .25 - n * min(.65,len(lengths)/max(len(labels),1))
        label = f"{n}. {str(row['GENE'])[:20]} · {str(r['SV_ID'])[:36]}"
        ax.annotate(label, xy=(r["START"] / 1e6, y_of[r["CHROM"]]),
                    xytext=(maximum * 1.04, label_y), fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": "#777777", "lw": .5})
    ax.set_yticks(list(y_of.values()), list(y_of))
    ax.set_xlim(0, maximum * (1.42 if not labels.empty else 1.03))
    ax.set_ylim(.15, len(lengths) + .7)
    ax.set_xlabel("Position on each chromosome (Mb; same scale for all chromosomes)")
    ax.set_title("Genome-wide SV positions", loc="left", fontweight="bold")
    handles = [Line2D([0], [0], color=color, lw=3, label=kind) for kind, color in SVTYPE_COLORS.items()]
    handles += [Line2D([0], [0], color=color, lw=3, label=name + (" (unavailable)" if bed.empty else ""))
                for (name, bed), color in zip(tracks.items(), track_colors)]
    height = fig.get_figheight()
    footer_height, legend_height = 36/(72*height), 32/(72*height)
    fig.legend(handles=handles, ncol=6, frameon=False, loc="lower left", bbox_to_anchor=(.06, footer_height+.035), fontsize=9)
    status = events["MAP_STATUS"].value_counts().to_dict()
    note = f"Each master SV is counted once; {status.get('PLOTTED', 0):,} anchors plotted; {len(events)-status.get('PLOTTED', 0):,} off-map/invalid/conflicting (see TSV)."
    note += "\nArcs connect reported BND/TRA partners; they do not reconstruct a rearrangement. Small events appear as points."
    note += "\nDepth ribbon: measured depth pooled into 1 Mb bins, not absolute copy number; use local depth/Samplot for individual DEL/DUP." if not depth.empty else "\nDepth ribbon unavailable: no binned measured coverage supplied."
    fig.text(.5, .015, note, ha="center", fontsize=9)
    fig.subplots_adjust(left=.065, right=.91, top=.95, bottom=footer_height+legend_height+.09)
    save_figure(fig, prefix)
    plt.close(fig)


def draw_rainfall(frame, lengths, prefix):
    fig, axes = plt.subplots(int(np.ceil(len(lengths) / 4)), 4, figsize=(16, max(5, np.ceil(len(lengths)/4)*2)), squeeze=False)
    for ax, chrom in zip(axes.flat, lengths):
        group = frame[frame["CHROM"].eq(chrom) & frame["DISTANCE_BP"].notna()]
        for kind, sub in group.groupby("SVTYPE"):
            ax.scatter(sub["START"] / 1e6, sub["DISTANCE_BP"], s=7, color=SVTYPE_COLORS.get(kind, "#999999"), rasterized=True)
        ax.set_yscale("symlog", linthresh=1)
        ax.set_xlim(0, lengths[chrom] / 1e6)
        ax.set_title(chrom, fontsize=10, loc="left")
        ax.set_xlabel("Position (Mb)", fontsize=8)
        ax.tick_params(labelsize=8)
        if group.empty:
            ax.text(.5, .5, "Fewer than two anchors", transform=ax.transAxes, ha="center", fontsize=8)
    for ax in list(axes.flat)[len(lengths):]:
        ax.axis("off")
    fig.supylabel("Distance from previous SV anchor (bp; log scale above 1 bp)")
    fig.suptitle("SV spacing along each chromosome", fontweight="bold")
    fig.text(.5, .01, "One START anchor per master SV; chromosomes restart independently. Zero-distance ties are retained.\nClusters may reflect repeats, mapping difficulties or biology; this plot cannot distinguish their cause.", ha="center", fontsize=9)
    fig.tight_layout(rect=[.02, .055, 1, .96])
    save_figure(fig, prefix)
    plt.close(fig)


def repeat_context(events, resources, lengths, mask, prefix, repeats=200, seed=27, window=100):
    """Compare breakpoint overlap with chromosome/geometry-matched random events.

    Both observed and random events must avoid the same mask. Random intervals
    retain their span; BND partners retain their respective chromosomes. This
    exploratory reference-genome null does not model caller/read callability.
    """
    rng = np.random.default_rng(seed)
    summary, sources = [], []
    for name, bed in resources.items():
        idx = interval_index(bed)
        for kind, group in events[events["MAP_STATUS"].eq("PLOTTED")].groupby("SVTYPE"):
            evaluable, observed = [], []
            for _, row in group.iterrows():
                points = endpoints(row)
                if kind == "BND" and row["PARTNER_STATUS"] != "RESOLVED":
                    continue
                if kind == "INS" or kind not in {"BND", "DEL", "DUP", "INV", "CNV"}:
                    span = 0
                else:
                    span = int(row["END"] - row["START"]) if kind != "BND" else None
                eligible = not any(overlaps(mask, chrom, max(0, pos-window), min(lengths[chrom], pos+window+1)) for chrom, pos in points)
                hit = any(overlaps(idx, chrom, max(0, pos-window), min(lengths[chrom], pos+window+1)) for chrom, pos in points)
                sources.append({"SV_ID": row["SV_ID"], "SVTYPE": kind, "REPEAT_CLASS": name,
                                "EVALUABLE_UNMASKED": eligible, "BREAKPOINT_OVERLAP": hit, "WINDOW_BP": window})
                if eligible:
                    evaluable.append((row, points, span))
                    observed.append(hit)
            if not evaluable:
                continue
            random_percent = []
            for _ in range(repeats):
                hits = []
                for row, points, span in evaluable:
                    accepted = False
                    for _attempt in range(1000):
                        if kind == "BND" and points[0][0] == points[1][0]:
                            distance = points[1][1] - points[0][1]
                            start = int(rng.integers(0, lengths[points[0][0]] - abs(distance)))
                            ordered = [start, start+abs(distance)] if distance >= 0 else [start+abs(distance), start]
                            random_points = [(points[0][0], position) for position in ordered]
                        elif kind == "BND":
                            random_points = [(chrom, int(rng.integers(0, lengths[chrom]))) for chrom, _pos in points]
                        else:
                            start = int(rng.integers(0, lengths[row["CHROM"]] - span))
                            random_points = [(row["CHROM"], start)] + ([(row["CHROM"], start+span)] if len(points)>1 else [])
                        if not any(overlaps(mask, chrom, max(0, pos-window), min(lengths[chrom], pos+window+1)) for chrom, pos in random_points):
                            accepted = True
                            break
                    if not accepted:
                        raise ValueError("Cannot draw an unmasked matched random event; check reference/mask")
                    hits.append(any(overlaps(idx, chrom, max(0, pos-window), min(lengths[chrom], pos+window+1)) for chrom, pos in random_points))
                random_percent.append(100*np.mean(hits))
            summary.append({"REPEAT_CLASS": name, "SVTYPE": kind, "N_EVALUABLE": len(evaluable),
                            "OBSERVED_PERCENT": 100*np.mean(observed), "RANDOM_MEAN_PERCENT": np.mean(random_percent),
                            "RANDOM_2_5_PERCENT": np.quantile(random_percent, .025), "RANDOM_97_5_PERCENT": np.quantile(random_percent, .975),
                            "RANDOM_DRAWS": repeats, "SEED": seed, "WINDOW_BP": window})
    pd.DataFrame(sources).to_csv(str(prefix)+"_breakpoints.tsv", sep="\t", index=False)
    frame = pd.DataFrame(summary)
    frame.to_csv(str(prefix)+"_summary.tsv", sep="\t", index=False)
    if frame.empty:
        return
    fig, ax = plt.subplots(figsize=(13, max(5, len(frame)*.34+2)))
    y = np.arange(len(frame))
    ax.barh(y, frame["OBSERVED_PERCENT"], color=[SVTYPE_COLORS.get(x,"#999999") for x in frame["SVTYPE"]], height=.7, label="Observed")
    center = frame["RANDOM_MEAN_PERCENT"].to_numpy()
    lo, hi = frame["RANDOM_2_5_PERCENT"].to_numpy(), frame["RANDOM_97_5_PERCENT"].to_numpy()
    ax.hlines(y, lo, hi, color="black", lw=2, label="Random 95% reference interval")
    ax.scatter(center, y, color="black", marker="|", s=100, label="Random mean")
    ax.set_yticks(y, [f"{r.REPEAT_CLASS} · {r.SVTYPE} (n={r.N_EVALUABLE})" for r in frame.itertuples()])
    ax.invert_yaxis()
    ax.set_xlim(0, 102)
    ax.set_xlabel(f"SVs with at least one breakpoint within {window} bp of a repeat (%)")
    ax.set_title("Repeat context compared with matched random positions", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0, -.12))
    fig.text(.5, .015, f"{repeats} draws; seed {seed}. Chromosome, SV type and interval span retained; observed/random endpoints use the same exclusion mask.\nClasses can overlap, so percentages are not additive. This reference-genome expectation does not model mappability or caller sensitivity; no mechanism is proven.", ha="center", fontsize=9)
    fig.subplots_adjust(left=.29, right=.97, top=.92, bottom=.24)
    save_figure(fig, prefix)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True)
    p.add_argument("--reference-fai", required=True)
    p.add_argument("--candidate-table", help="Completed review shortlist; labels use its existing order")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--label-top-n", type=int, default=10)
    for name in ["centromere-bed", "gap-bed", "exclude-bed", "panel-bed", "depth-bed"]:
        p.add_argument("--"+name)
    p.add_argument("--repeat-bed", action="append", default=[], metavar="CLASS=PATH", help="Optional BED; repeat classes stay separate")
    p.add_argument("--random-draws", type=int, default=200)
    p.add_argument("--seed", type=int, default=27)
    p.add_argument("--breakpoint-window", type=int, default=100)
    args = p.parse_args()
    if args.label_top_n < 0 or args.random_draws < 20 or args.breakpoint_window < 0:
        p.error("label-top-n/window must be nonnegative and random-draws at least 20")
    set_thesis_style()
    lengths = read_lengths(args.reference_fai)
    events = prepare_events(read_tsv(args.input), lengths)
    tracks = {"Centromeres": read_bed(args.centromere_bed, lengths), "Assembly gaps": read_bed(args.gap_bed, lengths),
              "DELLY excluded": read_bed(args.exclude_bed, lengths), "Panel genes": read_bed(args.panel_bed, lengths)}
    mask = interval_index(pd.concat(list(tracks.values())[:3], ignore_index=True))
    labels = pd.DataFrame()
    if args.candidate_table and args.label_top_n:
        selected = read_tsv(args.candidate_table)
        labels = select_events(selected, args.label_top_n).drop_duplicates("SV_ID").head(args.label_top_n)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, bed in tracks.items():
        idx = interval_index(bed)
        events[name.upper().replace(" ", "_")+"_BREAKPOINT"] = [
            any(overlaps(idx, chrom, pos, pos+1) for chrom, pos in endpoints(row)) if row["MAP_STATUS"] == "PLOTTED" and not bed.empty else pd.NA
            for _, row in events.iterrows()]
        bed.to_csv(str(prefix)+"_"+name.lower().replace(" ", "_")+".tsv", sep="\t", index=False)
    depth = load_depth(args.depth_bed, lengths, mask)
    fields = [col for col in ["SV_ID","CHROM","START","END","SVTYPE","SVLEN","CHR2","POS2",
                              "CALLER_COUNT","NEEDLR_AF","GNOMAD_SV_AF","MAP_STATUS","PARTNER_STATUS"] if col in events]
    fields += [col for col in events if col.endswith("_BREAKPOINT")]
    events[fields].to_csv(str(prefix)+"_events.tsv", sep="\t", index=False)
    depth.to_csv(str(prefix)+"_depth.tsv", sep="\t", index=False)
    spacing = rainfall_table(events)
    spacing.to_csv(str(prefix)+"_rainfall.tsv", sep="\t", index=False)
    draw_map(events, lengths, tracks, depth, labels, prefix)
    draw_rainfall(spacing, lengths, str(prefix)+"_rainfall")
    resources = {}
    for item in args.repeat_bed:
        name, path = item.split("=", 1)
        if name in resources:
            raise ValueError(f"Repeat class specified twice: {name}")
        resources[name] = read_bed(path, lengths)
    if resources:
        repeat_context(events, resources, lengths, mask, str(prefix)+"_repeat_context", args.random_draws, args.seed, args.breakpoint_window)
    metadata = vars(args) | {"reference_lengths": lengths, "sv_coordinates": "VCF_1_BASED", "bed_coordinates": "0_BASED_HALF_OPEN",
                              "depth_display_bin_bp": 1000000, "map_status_counts": events["MAP_STATUS"].value_counts().to_dict()}
    Path(str(prefix)+"_settings.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"[OK] {len(events)} master SVs; genome map, spacing and available context: {prefix}")


if __name__ == "__main__":
    main()
