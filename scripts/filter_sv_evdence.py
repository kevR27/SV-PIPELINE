#!/usr/bin/env python3
"""
filter_sv_evidence.py

Applies configurable, transparent pre-merge quality control to a single
caller's SV evidence table (the TSV produced by parse_sv_caller_vcf.py)
before that caller's VCF is handed to Jasmine for merging.

Design principle: hard-filter only on evidence that is cheap, caller-agnostic,
and unlikely to discard a real clinical SV (FILTER, read support, SV size).
Everything softer / more caller-specific (precision flag, low genotype
quality, overlap with a low-mappability/segdup blacklist) is written as an
annotation column rather than used to drop the call, so a human reviewer can
still see and evaluate it in the final integrated table.

Outputs:
  --output-tsv   the input table with two new columns: EVIDENCE_STATUS
                 (PASS/FAIL) and EVIDENCE_FLAGS (semicolon-separated soft flags)
  --output-ids   a plain-text list of SV_IDs with EVIDENCE_STATUS=PASS,
                 one per line, suitable for `bcftools view -i 'ID=@file'`
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

MISSING = "."


def to_float(v):
    if v in (None, "", MISSING):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def load_blacklist(bed_path):
    """Minimal BED interval loader: {chrom: [(start,end), ...]} sorted by start."""
    intervals = {}
    if not bed_path:
        return intervals
    with open(bed_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def overlaps_blacklist(intervals, chrom, pos, end):
    ivs = intervals.get(chrom)
    if not ivs:
        return False
    pos, end = int(pos), int(end)
    lo, hi = min(pos, end), max(pos, end)
    # linear scan is fine at panel/sample scale; swap for bisect if this
    # becomes a bottleneck on genome-wide callsets with a large blacklist
    for s, e in ivs:
        if s > hi:
            break
        if e >= lo:
            return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description="Pre-merge SV evidence filter.")
    ap.add_argument("--caller-tsv", required=True, help="Output of parse_sv_caller_vcf.py")
    ap.add_argument("--caller", required=True, help="Caller name, for logging only")
    ap.add_argument("--min-support", type=float, default=2,
                     help="Minimum CALLER_SUPPORT to pass (default: 2)")
    ap.add_argument("--min-svlen", type=float, default=50,
                     help="Minimum absolute SV length in bp to pass (default: 50)")
    ap.add_argument("--max-svlen", type=float, default=5_000_000,
                     help="Maximum absolute SV length in bp to pass (default: 5,000,000); "
                          "set higher/None-like huge value for panels expecting large CNVs")
    ap.add_argument("--require-pass", action="store_true", default=True,
                     help="Require FILTER in {PASS,.} to pass (default: on)")
    ap.add_argument("--allow-any-filter", dest="require_pass", action="store_false",
                     help="Disable the FILTER=PASS/. requirement")
    ap.add_argument("--min-gq", type=float, default=None,
                     help="Optional: flag (not drop) records below this genotype quality")
    ap.add_argument("--blacklist-bed", default=None,
                     help="Optional BED of low-mappability/segdup regions; overlap is "
                          "flagged, not excluded")
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--output-ids", required=True)
    a = ap.parse_args()

    blacklist = load_blacklist(a.blacklist_bed)

    with open(a.caller_tsv, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    n_total = len(rows)
    n_pass = 0
    for row in rows:
        flags = []
        fail_reasons = []

        filt = row.get("FILTER", MISSING)
        if a.require_pass and filt not in ("PASS", MISSING, "."):
            fail_reasons.append("NON_PASS_FILTER")

        support = to_float(row.get("CALLER_SUPPORT"))
        if support is None or support < a.min_support:
            fail_reasons.append("LOW_SUPPORT")

        svlen = to_float(row.get("SVLEN"))
        svlen_abs = abs(svlen) if svlen is not None else None
        if svlen_abs is None:
            # BND / TRA and similar records often have no SVLEN; don't fail
            # them on size, just note it
            flags.append("NO_SVLEN")
        elif svlen_abs < a.min_svlen or svlen_abs > a.max_svlen:
            fail_reasons.append("SIZE_OUT_OF_RANGE")

        imprecise = row.get("CALLER_IMPRECISE", MISSING)
        precise = row.get("CALLER_PRECISE", MISSING)
        if imprecise not in (MISSING, "", "False") or precise in ("false", "False", "0"):
            flags.append("IMPRECISE")

        gq = to_float(row.get("CALLER_GQ"))
        if a.min_gq is not None and gq is not None and gq < a.min_gq:
            flags.append("LOW_GQ")

        if blacklist:
            chrom = row.get("CHROM", MISSING)
            start = row.get("START", MISSING)
            end = row.get("END", MISSING)
            end = end if end not in (MISSING, "") else start
            if chrom not in (MISSING, "") and start not in (MISSING, ""):
                try:
                    if overlaps_blacklist(blacklist, chrom, start, end):
                        flags.append("BLACKLIST_REGION")
                except ValueError:
                    pass

        status = "FAIL" if fail_reasons else "PASS"
        if status == "PASS":
            n_pass += 1
        row["EVIDENCE_STATUS"] = status
        row["EVIDENCE_FAIL_REASONS"] = ";".join(fail_reasons) if fail_reasons else MISSING
        row["EVIDENCE_FLAGS"] = ";".join(flags) if flags else MISSING

    out_cols = fieldnames + ["EVIDENCE_STATUS", "EVIDENCE_FAIL_REASONS", "EVIDENCE_FLAGS"]
    out_tsv = Path(a.output_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    out_ids = Path(a.output_ids)
    out_ids.parent.mkdir(parents=True, exist_ok=True)
    with out_ids.open("w", encoding="utf-8") as fh:
        for row in rows:
            if row["EVIDENCE_STATUS"] == "PASS":
                fh.write(row["SV_ID"] + "\n")

    print(f"[OK] caller={a.caller} total={n_total} pass={n_pass} "
          f"fail={n_total - n_pass} output_tsv={out_tsv} output_ids={out_ids}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
