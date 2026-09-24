#!/usr/bin/env python3
"""Attach methylation and nearby small-variant phasing context to master SV rows.

This script does not treat methylation or WhatsHap small variants as SV
confirmations. It adds local breakpoint-context measurements to the already
integrated structural-variant table.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from plot_utils import first_existing, read_tsv


def parse_args():
    p = argparse.ArgumentParser(description="Add methylation and WhatsHap context to integrated SVs.")
    p.add_argument("--integrated", required=True)
    p.add_argument("--methylation-bed", default=None)
    p.add_argument("--whatshap-vcf", default=None)
    p.add_argument("--methylation-window", type=int, default=2000)
    p.add_argument("--whatshap-window", type=int, default=5000)
    p.add_argument("--min-methylation-coverage", type=int, default=5)
    p.add_argument("--output", required=True)
    return p.parse_args()


def normalize_chrom(value):
    text = str(value)
    return text if text.startswith("chr") else "chr" + text


def to_int(value):
    try:
        return int(float(value))
    except Exception:
        return None


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw):
    out = {}
    for item in str(raw).split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def load_whatshap(path):
    """Index phased heterozygous small variants by chromosome and position."""
    by_chrom = {}
    temp = {}
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
            if "|" not in gt:
                continue
            alleles = gt.replace("|", "/").split("/")
            if len(alleles) != 2 or "." in alleles or alleles[0] == alleles[1]:
                continue
            chrom = normalize_chrom(fields[0])
            pos = to_int(fields[1])
            if pos is None:
                continue
            temp.setdefault(chrom, []).append((pos, data.get("PS", "."), gt))

    for chrom, records in temp.items():
        records.sort(key=lambda x: x[0])
        by_chrom[chrom] = {
            "records": records,
            "positions": [x[0] for x in records],
        }
    return by_chrom


def query_whatshap(index, chrom, centers, window):
    if chrom not in index:
        return 0, set()
    positions = index[chrom]["positions"]
    records = index[chrom]["records"]
    hits = {}
    for center in centers:
        if center is None:
            continue
        lo = bisect.bisect_left(positions, max(0, center - window))
        hi = bisect.bisect_right(positions, center + window)
        for pos, ps, gt in records[lo:hi]:
            hits[(pos, ps, gt)] = True
    phase_sets = {ps for _, ps, _ in hits if ps not in {"", ".", "NA"}}
    return len(hits), phase_sets


def methylation_scale(path):
    values = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                values.append(float(fields[10]))
            except Exception:
                continue
            if len(values) >= 1000:
                break
    if values and max(values) <= 1.0:
        return 100.0
    return 1.0


def normalize_mod_code(value):
    text = str(value).strip().lower()
    if text in {"m", "5mc", "c+m"}:
        return "5mC"
    if text in {"h", "5hmc", "c+h"}:
        return "5hmC"
    return str(value).strip() or "other"


def open_tabix(path):
    try:
        import pysam
    except ImportError:
        return None
    index_candidates = [Path(str(path) + ".tbi"), Path(str(path) + ".csi")]
    if not any(p.exists() for p in index_candidates):
        return None
    try:
        return pysam.TabixFile(str(path))
    except Exception:
        return None


def parse_bedmethyl_record(line, scale, min_cov):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return None
    try:
        start = int(fields[1])
        end = int(fields[2])
        coverage = float(fields[9])
        percent = float(fields[10]) * scale
    except Exception:
        return None
    if coverage < min_cov:
        return None
    return {
        "start": start,
        "end": end,
        "mod": normalize_mod_code(fields[3]),
        "coverage": coverage,
        "percent": percent,
    }


def query_methylation_tabix(tabix, chrom, centers, window, scale, min_cov):
    records = {}
    for center in centers:
        if center is None:
            continue
        start = max(0, center - window)
        end = center + window + 1
        try:
            fetched = tabix.fetch(chrom, start, end)
        except Exception:
            continue
        for line in fetched:
            rec = parse_bedmethyl_record(line, scale, min_cov)
            if rec is None:
                continue
            key = (rec["start"], rec["end"], rec["mod"])
            records[key] = rec
    return list(records.values())


def summarize_methylation(records):
    result = {
        "METHYLATION_CONTEXT": "NO_CPG",
        "METHYLATION_CPG_RECORD_COUNT": 0,
        "METHYLATION_MEAN_COVERAGE": ".",
        "METHYLATION_5MC_CPG_COUNT": 0,
        "METHYLATION_5MC_MEAN_PERCENT": ".",
        "METHYLATION_5HMC_CPG_COUNT": 0,
        "METHYLATION_5HMC_MEAN_PERCENT": ".",
    }
    if not records:
        return result

    result["METHYLATION_CONTEXT"] = "EVALUATED"
    result["METHYLATION_CPG_RECORD_COUNT"] = len(records)
    result["METHYLATION_MEAN_COVERAGE"] = f"{np.mean([r['coverage'] for r in records]):.3f}"

    for mod, prefix in [("5mC", "5MC"), ("5hmC", "5HMC")]:
        vals = [r["percent"] for r in records if r["mod"] == mod]
        result[f"METHYLATION_{prefix}_CPG_COUNT"] = len(vals)
        result[f"METHYLATION_{prefix}_MEAN_PERCENT"] = f"{np.mean(vals):.3f}" if vals else "."
    return result


def main():
    args = parse_args()
    df = read_tsv(args.integrated)
    if df.empty:
        raise ValueError("Integrated table is empty.")

    id_col = first_existing(df, ["SV_ID", "ID"])
    chrom_col = first_existing(df, ["CHROM", "chrom"])
    start_col = first_existing(df, ["START", "POS"])
    end_col = first_existing(df, ["END"])
    chr2_col = first_existing(df, ["CHR2"])
    pos2_col = first_existing(df, ["POS2"])
    if None in (id_col, chrom_col, start_col):
        raise ValueError("Integrated table needs SV_ID, CHROM and START/POS.")

    whatshap = load_whatshap(args.whatshap_vcf) if args.whatshap_vcf else {}
    methylation_path = Path(args.methylation_bed) if args.methylation_bed else None
    methylation_available = bool(methylation_path and methylation_path.exists())
    scale = methylation_scale(methylation_path) if methylation_available else 1.0
    tabix = open_tabix(methylation_path) if methylation_available else None

    out = df.copy()
    defaults = {
        "WHATSHAP_CONTEXT": "NOT_AVAILABLE" if not args.whatshap_vcf else "NO_NEARBY_PHASED_HET",
        "WHATSHAP_WINDOW_BP": args.whatshap_window,
        "WHATSHAP_PHASED_HET_COUNT": 0,
        "WHATSHAP_PHASE_SET_COUNT": 0,
        "WHATSHAP_PHASE_SETS": ".",
        "METHYLATION_CONTEXT": "NOT_AVAILABLE" if not methylation_available else ("INDEX_MISSING" if tabix is None else "NO_CPG"),
        "METHYLATION_WINDOW_BP": args.methylation_window,
        "METHYLATION_CPG_RECORD_COUNT": 0,
        "METHYLATION_MEAN_COVERAGE": ".",
        "METHYLATION_5MC_CPG_COUNT": 0,
        "METHYLATION_5MC_MEAN_PERCENT": ".",
        "METHYLATION_5HMC_CPG_COUNT": 0,
        "METHYLATION_5HMC_MEAN_PERCENT": ".",
    }
    for col, value in defaults.items():
        out[col] = value

    cache = {}
    for idx, row in out.iterrows():
        sv_id = str(row[id_col])
        if sv_id in cache:
            context = cache[sv_id]
        else:
            chrom = normalize_chrom(row[chrom_col])
            start = to_int(row[start_col])
            end = to_int(row[end_col]) if end_col else start
            centers = [start]
            if end is not None and end != start:
                centers.append(end)

            context = {}
            if args.whatshap_vcf:
                count, phase_sets = query_whatshap(whatshap, chrom, centers, args.whatshap_window)
                context.update(
                    {
                        "WHATSHAP_CONTEXT": "EVALUATED" if count else "NO_NEARBY_PHASED_HET",
                        "WHATSHAP_PHASED_HET_COUNT": count,
                        "WHATSHAP_PHASE_SET_COUNT": len(phase_sets),
                        "WHATSHAP_PHASE_SETS": ";".join(sorted(phase_sets)) if phase_sets else ".",
                    }
                )

            if tabix is not None:
                records = query_methylation_tabix(
                    tabix,
                    chrom,
                    centers,
                    args.methylation_window,
                    scale,
                    args.min_methylation_coverage,
                )
                if chr2_col and pos2_col:
                    chr2 = row[chr2_col]
                    pos2 = to_int(row[pos2_col])
                    if pos2 is not None and str(chr2) not in {"", ".", "nan", "None"}:
                        records += query_methylation_tabix(
                            tabix,
                            normalize_chrom(chr2),
                            [pos2],
                            args.methylation_window,
                            scale,
                            args.min_methylation_coverage,
                        )
                unique = {}
                for rec in records:
                    unique[(rec["start"], rec["end"], rec["mod"], rec["percent"])] = rec
                context.update(summarize_methylation(list(unique.values())))

            cache[sv_id] = context

        for key, value in context.items():
            out.at[idx, key] = value

    if tabix is not None:
        tabix.close()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    unique = out.drop_duplicates(id_col)
    methyl_eval = int((unique["METHYLATION_CONTEXT"] == "EVALUATED").sum())
    wh_eval = int((unique["WHATSHAP_CONTEXT"] == "EVALUATED").sum())
    print(
        f"[OK] master_SVs={len(unique)} methylation_context={methyl_eval} "
        f"whatshap_nearby_phased={wh_eval} output={output}"
    )
    if methylation_available and tabix is None:
        print(
            "[WARN] methylation file was found but no .tbi/.csi index was usable; "
            "global methylation plotting will still work, but per-SV methylation context was not queried."
        )


if __name__ == "__main__":
    main()
