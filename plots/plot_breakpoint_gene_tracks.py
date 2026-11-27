#!/usr/bin/env python3
"""Compact measured-coverage and SV-breakpoint gene-track figures.

Original graphics inspired by scientific genomic track conventions (not
a reproduction of external images). No BAM, Samplot or simulated coverage.
Precise caller endpoints are shown ONLY when independently matched to a
record in an accessible caller VCF; missing traces stay unassessed.
"""
import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def chrom_name(value):
    text = str(value).strip()
    return text if text.startswith("chr") else "chr" + text


def number(value):
    try:
        n = float(value)
        return int(n) if np.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def choose_candidates(path, top_n):
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    required = {"SV_ID", "SVTYPE"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing columns in final SV candidate table: {sorted(required-set(frame.columns))}")
    frame = frame.copy()
    rank = next((name for name in ["FINAL_EVENT_RANK_GLOBAL", "EVENT_RANK_GLOBAL"]
                 if name in frame), None)
    if rank:
        frame["_rank"] = pd.to_numeric(frame[rank], errors="coerce")
        frame = frame.sort_values("_rank", na_position="last")
    frame = frame.drop_duplicates("SV_ID")
    records = []
    for _, row in frame.iterrows():
        c = next((row[x] for x in ["CHROM", "CHR", "CHROMOSOME"]
                  if x in row and pd.notna(row[x])), None)
        start = next((number(row[x]) for x in ["START", "POS", "SV_START"]
                      if x in row and number(row[x]) is not None), None)
        end = next((number(row[x]) for x in ["END", "SV_END"]
                    if x in row and number(row[x]) is not None), None)
        svtype = str(row["SVTYPE"]).upper()
        if c is None or start is None or end is None or end < start:
            continue
        if svtype not in {"DEL", "DUP", "INV", "INS"}:
            continue
        genes = str(row.get("GENE", row.get("GENES", "."))).split(";")
        records.append(dict(SV_ID=str(row["SV_ID"]), chrom=chrom_name(c),
                            start=start, end=end, svtype=svtype,
                            gene=genes[0], row=row))
        if len(records) >= top_n:
            break
    return records


def coverage_tracks(path, candidates):
    tracks = defaultdict(list)
    if not path or not Path(path).exists():
        return tracks
    windows = {}
    for i, event in enumerate(candidates):
        size = max(event["end"]-event["start"], 1)
        pad = min(max(20000, size // 3), 500000)
        windows[i] = (event["chrom"], max(0,event["start"]-pad), event["end"]+pad)
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            try:
                chrom, a, b, d = chrom_name(parts[0]), int(parts[1]), int(parts[2]), float(parts[3])
            except (ValueError, IndexError):
                continue
            if not np.isfinite(d) or d < 0:
                continue
            for i, (wc, left, right) in windows.items():
                if chrom == wc and b > left and a < right:
                    tracks[i].append((a, b, d))
    return tracks


def gene_models(path, events):
    models = defaultdict(list)
    if not path or not Path(path).is_file():
        return models
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                chrom, a, b = chrom_name(f[0]), int(f[1]), int(f[2])
            except ValueError:
                continue
            for idx, e in enumerate(events):
                span = max(20000, int((e["end"]-e["start"])*0.25))
                if chrom == e["chrom"] and a <= e["end"]+span and b >= e["start"]-span:
                    models[idx].append((a,b,f[3]))
    return models


def exon_models(path, events):
    """True exons from GTF; BED4 tracks are only whole-gene spans."""
    exons = defaultdict(list)
    if not path or not Path(path).is_file():
        return exons
    wanted = {e["gene"] for e in events}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "exon":
                continue
            gene = re.search(r'gene_name\s+"([^"]+)"', f[8])
            if not gene or gene.group(1) not in wanted:
                continue
            t = re.search(r'transcript_id\s+"([^"]+)"', f[8])
            try:
                chrom, a, b = chrom_name(f[0]), int(f[3])-1, int(f[4])
            except ValueError:
                continue
            for i,e in enumerate(events):
                if chrom == e["chrom"] and gene.group(1) == e["gene"] and b >= e["start"]-500000 and a <= e["end"]+500000:
                    exons[i].append((a,b,t.group(1) if t else "."))
    return exons


def raw_callers(sample_dir, platform):
    s = sample_dir.name
    if platform == "lrs":
        patterns = {
            "Sniffles2": f"sv/sniffles2/{s}_sniffles.vcf*",
            "cuteSV": f"sv/cutesv/{s}_cutesv.vcf*",
            "DELLY": f"sv/delly/{s}_delly.vcf*",
            "pbsv": f"sv/pbsv/{s}_pbsv.vcf*",
            "SVIM": f"sv/svim/{s}_svim.vcf*",
        }
    else:
        patterns = {
            "Manta": f"sv/manta/{s}_manta.vcf*",
            "DELLY": f"sv/delly/{s}_delly.vcf*",
            "GRIDSS": f"sv/gridss/{s}*vcf*",
        }
    result = {}
    for label, pattern in patterns.items():
        matches = [p for p in sample_dir.glob(pattern)
                   if p.is_file() and p.suffix not in {".tbi", ".csi"}]
        if matches:
            result[label] = sorted(matches, key=lambda x:(not str(x).endswith(".vcf.gz"), len(x.name)))[0]
    return result


def caller_breakpoints(samples, events):
    """Read each caller file once; match same chr/type and both breakpoints.

    Multiple nearby calls are not assumed independent or equivalent.
    """
    matches = defaultdict(dict)
    for name, path in samples.items():
        best = {}
        with open_text(path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                chrom = chrom_name(f[0])
                start = number(f[1])
                info = {}
                for item in f[7].split(";"):
                    if "=" in item:
                        key, val = item.split("=", 1)
                        info[key] = val
                svtype = info.get("SVTYPE", "")
                if not svtype and f[4].startswith("<"):
                    svtype = f[4].strip("<>").split(":")[0]
                end = number(info.get("END"))
                if not svtype or start is None or end is None or end < start:
                    continue
                for i, event in enumerate(events):
                    if chrom != event["chrom"] or svtype.upper() != event["svtype"]:
                        continue
                    tolerance = min(2000, max(250,int((event["end"]-event["start"])*.01)))
                    delta = abs(event["start"]-start)+abs(event["end"]-end)
                    if abs(event["start"]-start) <= tolerance and abs(event["end"]-end) <= tolerance:
                        if i not in best or delta < best[i][0]:
                            best[i] = (delta,start,end)
        for i, (_,start,end) in best.items():
            matches[i][name] = (start,end)
    return matches


def draw(event, depth, genes, exons, callers, outprefix):
    start, end = event["start"], event["end"]
    span = max(1,end-start)
    padding = min(max(20000, span//3),500000)
    xleft, xright = max(0,start-padding),end+padding
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 5.4), sharex=True,
                             gridspec_kw={"height_ratios":[2.0,1.3]})
    cov, ax = axes
    if depth:
        depth = sorted(depth)
        xs = [(a+b)/2 for a,b,v in depth]
        ys = [v for a,b,v in depth]
        cov.fill_between(xs, ys, step="mid", color="#1B7517", alpha=.92)
        cov.set_ylabel("Measured read depth")
        overlap = [v for a,b,v in depth if b > start and a < end]
        flanks = [v for a,b,v in depth if b <= start or a >= end]
        # A single coarse 10-kb bin cannot resolve a small exon-scale SV.
        if len(overlap) >= 2 and len(flanks) >= 2:
            ratio = np.median(overlap)/max(np.median(flanks),1e-10)
            cov.text(.01,.95, f"Inside/flank depth ratio: {ratio:.2f} (descriptive)",
                     transform=cov.transAxes, ha="left", va="top", fontsize=9)
        elif overlap:
            cov.text(.01,.95,"Coverage bin size limits CNV interpretation",
                     transform=cov.transAxes,ha="left",va="top",fontsize=9)
    else:
        cov.text(.5,.5,"Measured binned depth unavailable",transform=cov.transAxes,
                 ha="center",va="center")
        cov.set_ylabel("Depth unavailable")
    cov.axvspan(start,end,color="#777777",alpha=.09)
    cov.plot([start,end], [1.035,1.035], transform=cov.get_xaxis_transform(),
             linestyle="--",color="gray",clip_on=False)
    cov.scatter([start,end], [1.035,1.035],transform=cov.get_xaxis_transform(),
                marker="D",s=32,color="gray",clip_on=False)
    for j,(name,(a,b)) in enumerate(sorted(callers.items())):
        y = 1.14+j*.08
        cov.plot([a,b],[y,y],transform=cov.get_xaxis_transform(),alpha=.85,
                 linewidth=1.2,linestyle="--",clip_on=False)
        cov.scatter([a,b],[y,y], transform=cov.get_xaxis_transform(),s=24,marker="|",
                    clip_on=False)
        cov.text(b,y,name,transform=cov.get_xaxis_transform(),fontsize=8,
                 va="center",ha="left",clip_on=False)
    cov.set_title(f'{event["svtype"]} {event["chrom"]}:{start:,}-{end:,} | {event["gene"]}',
                  fontsize=12,weight="bold",pad=32+len(callers)*11)
    shown = 0
    for a,b,g in sorted(genes,key=lambda z:(z[0],z[2])):
        if shown >= 7:
            break
        if b < xleft or a > xright:
            continue
        y = 0.9-shown*.13
        ax.plot([a,b],[y,y],color="#302724",linewidth=2)
        gx = [e for e in exons if a <= e[0] and e[1] <= b]
        # Multiple isoforms can overlap: show the union, not one invented transcript.
        for u,v in sorted(set((x,y) for x,y,_ in gx)):
            ax.add_patch(plt.Rectangle((u,y-.027),max(1,v-u),.054,color="#302724"))
        ax.text((a+b)/2,y-.06,g,ha="center",va="top",fontsize=9,weight="bold")
        shown+=1
    if not shown:
        ax.text(.5,.5,"No local gene BED annotation available",
                transform=ax.transAxes,ha="center")
    ax.axvspan(start,end,color="#777777",alpha=.11)
    ax.set_ylim(-.12,1.2)
    ax.set_yticks([])
    ax.set_xlabel(f'{event["chrom"]} position (GRCh38, bp)')
    ax.set_xlim(xleft,xright)
    fig.text(.5,.012,
             "Diamond endpoints: merged SV. Dashed caller endpoints: VCF interval matches, not read-level validation. "
             "Exon rectangles require a supplied GTF; otherwise only whole-gene spans are shown.",
             ha="center",fontsize=8)
    fig.tight_layout(rect=(0,.06,1,.95))
    outprefix.parent.mkdir(parents=True,exist_ok=True)
    for ext in (".pdf",".png"):
        fig.savefig(str(outprefix)+ext,dpi=240,bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--platform", choices=["lrs","srs"], required=True)
    parser.add_argument("--gene-bed")
    parser.add_argument("--gene-gtf")
    parser.add_argument("--depth")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--top-n",type=int,default=10)
    args = parser.parse_args()
    if args.top_n < 1:
        parser.error("--top-n must be at least 1")
    sample_dir=Path(args.sample_dir)
    candidates=choose_candidates(args.candidates,args.top_n)
    folder=Path(args.out_dir)
    folder.mkdir(parents=True,exist_ok=True)
    depth=coverage_tracks(args.depth,candidates)
    genes=gene_models(args.gene_bed,candidates)
    exons=exon_models(args.gene_gtf,candidates)
    caller_sources=raw_callers(sample_dir,args.platform)
    endpoints=caller_breakpoints(caller_sources,candidates)
    rows=[]
    for i,event in enumerate(candidates):
        name=f'{i+1:02d}_{event["svtype"]}_{event["chrom"]}_{event["start"]}_{event["end"]}'
        draw(event, depth.get(i,[]),genes.get(i,[]),exons.get(i,[]),
             endpoints.get(i,{}),folder/name)
        rows.append({
            "SV_ID":event["SV_ID"],"GENE":event["gene"],
            "SVTYPE":event["svtype"],"CHROM":event["chrom"],
            "START":event["start"],"END":event["end"],
            "BINNED_DEPTH_STATUS":"AVAILABLE" if depth.get(i) else "NOT_AVAILABLE",
            "CALLER_ENDPOINT_MATCHES":";".join(sorted(endpoints.get(i,{}))) or ".",
            "CALLER_VCF_FILES_SEARCHED":";".join(sorted(caller_sources)) or ".",
            "EXONS_SHOWN":"YES" if exons.get(i) else "NO_GTF_OR_NO_MATCH",
            "FIGURE":name+".pdf",
        })
    pd.DataFrame(rows,columns=[
        "SV_ID","GENE","SVTYPE","CHROM","START","END","BINNED_DEPTH_STATUS",
        "CALLER_ENDPOINT_MATCHES","CALLER_VCF_FILES_SEARCHED","EXONS_SHOWN","FIGURE"
    ]).to_csv(folder/"breakpoint_tracks_manifest.tsv",sep="\t",index=False)
    print(f"[OK] breakpoint gene tracks={len(rows)} at {folder}")


if __name__=="__main__":
    main()
