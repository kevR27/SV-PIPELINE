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
  --output-tsv   the input table with EVIDENCE_STATUS, EVIDENCE_FAIL_REASONS,
                 and EVIDENCE_FLAGS appended
  --output-ids   a plain-text list of SV_IDs with EVIDENCE_STATUS=PASS,
                 one per line, suitable for `bcftools view -i 'ID=@file'`
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path
from sv_evidence_common import breakend, number, sv_length


def set_csv_field_size_limit() -> None:
    """Raise Python's CSV field-size limit as high as the platform permits."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


set_csv_field_size_limit()

MISSING = "."

CANONICAL_CHROMS = {
    *(f"chr{i}" for i in range(1, 23)),
    "chrX",
    "chrY",
    "chrM",
}


def normalize_chrom(value):
    text = str(value or "").strip()
    if text in {"", MISSING}:
        return None
    if text.startswith("chr"):
        return text
    if text == "MT":
        return "chrM"
    return "chr" + text


def to_float(v):
    if v in (None, "", MISSING):
        return None
    try:
        return number(v)
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
    lo, hi = min(pos, end) - 1, max(pos, end)
    # linear scan is fine at panel/sample scale; swap for bisect if this
    # becomes a bottleneck on genome-wide callsets with a large blacklist
    for s, e in ivs:
        if s >= hi:
            break
        if e > lo:
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
    ap.add_argument("--max-svlen", type=float, default=None,
                     help="Optional maximum absolute SV length in bp. "
                          "Omit this option to apply no upper SV-length filter.")
    ap.add_argument("--require-pass", action="store_true", default=True,
                     help="Require FILTER in {PASS,.} to pass (default: on)")
    ap.add_argument("--allow-any-filter", dest="require_pass", action="store_false",
                     help="Disable the FILTER=PASS/. requirement")
    ap.add_argument(
        "--rescue-cov-var",
        action="store_true",
        help=(
            "For Sniffles2 only, allow FILTER=COV_VAR calls of any SVTYPE "
            "to bypass the FILTER=PASS gate. They must still satisfy the "
            "normal support, size and chromosome filters and remain flagged "
            "RESCUED_COV_VAR for manual review."
        ),
    )
    ap.add_argument("--min-gq", type=float, default=None,
                     help="Optional: flag (not drop) records below this genotype quality")
    ap.add_argument("--blacklist-bed", default=None,
                     help="Optional BED of low-mappability/segdup regions; overlap is "
                          "flagged, not excluded")
    ap.add_argument(
        "--canonical-only",
        action="store_true",
        help="Keep only chr1-22, chrX, chrY and chrM. For BND/TRA both breakends must be canonical.",
    )
    ap.add_argument("--output-tsv", required=True)
    ap.add_argument("--output-ids", required=True)
    a = ap.parse_args()

    blacklist = load_blacklist(a.blacklist_bed)
    out_tsv = Path(a.output_tsv)
    out_ids = Path(a.output_ids)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_ids.parent.mkdir(parents=True, exist_ok=True)

    # Write to temporary files and atomically replace the requested outputs
    # only after the entire input TSV has been processed successfully. This
    # avoids leaving apparently complete partial outputs after an interruption.
    tmp_tsv = Path(str(out_tsv) + ".tmp")
    tmp_ids = Path(str(out_ids) + ".tmp")

    n_total = 0
    n_pass = 0
    n_rescued_cov_var = 0

    try:
        with (
            open(a.caller_tsv, "r", encoding="utf-8", errors="replace", newline="") as in_fh,
            tmp_tsv.open("w", encoding="utf-8", newline="") as out_fh,
            tmp_ids.open("w", encoding="utf-8") as ids_fh,
        ):
            reader = csv.DictReader(in_fh, delimiter="\t")
            fieldnames = list(reader.fieldnames or [])
            out_cols = list(dict.fromkeys(fieldnames + ["SVLEN_REPORTED", "SVLEN_SOURCE"])) + [
                "EVIDENCE_STATUS",
                "EVIDENCE_FAIL_REASONS",
                "EVIDENCE_FLAGS",
            ]
            writer = csv.DictWriter(
                out_fh,
                fieldnames=out_cols,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in reader:
                n_total += 1
                flags = []
                fail_reasons = []

                filt = row.get("FILTER", MISSING)
                support = to_float(row.get("CALLER_SUPPORT"))
                svtype = (row.get("SVTYPE", MISSING) or MISSING).upper()
                row.setdefault("SVLEN_REPORTED", row.get("SVLEN", MISSING))
                length, source = sv_length(row)
                row["SVLEN"] = length
                if row.get("SVLEN_SOURCE", MISSING) in (MISSING, ""):
                    row["SVLEN_SOURCE"] = source
                svlen = to_float(length)
                svlen_abs = abs(svlen) if svlen is not None else None

                rescued_cov_var = (
                    a.rescue_cov_var
                    and a.caller.lower() == "sniffles2"
                    and filt == "COV_VAR"
                )

                if a.require_pass and filt not in ("PASS", MISSING, "."):
                    if rescued_cov_var:
                        flags.append("RESCUED_COV_VAR")
                    else:
                        fail_reasons.append("NON_PASS_FILTER")

                if support is None or support < a.min_support:
                    fail_reasons.append("LOW_SUPPORT")

                if a.canonical_only:
                    chrom1 = normalize_chrom(row.get("CHROM"))
                    if chrom1 not in CANONICAL_CHROMS:
                        fail_reasons.append("NON_CANONICAL_CHROM")
                    if svtype in {"BND", "TRA"}:
                        chrom2 = normalize_chrom(row.get("CHR2"))
                        if chrom2 is None:
                            fail_reasons.append("UNRESOLVED_BREAKEND_CHROM")
                        elif chrom2 not in CANONICAL_CHROMS:
                            fail_reasons.append("NON_CANONICAL_BREAKEND_CHROM")

                if svtype in {"BND", "TRA"}:
                    flags.append("NO_SVLEN")
                elif svlen_abs is None:
                    # BND / TRA and similar records often have no SVLEN; do
                    # not fail them on size, just record the missing length.
                    flags.append("NO_SVLEN")
                elif svlen_abs < a.min_svlen:
                    fail_reasons.append("SIZE_BELOW_MIN")
                elif a.max_svlen is not None and svlen_abs > a.max_svlen:
                    fail_reasons.append("SIZE_ABOVE_MAX")

                imprecise = row.get("CALLER_IMPRECISE", MISSING)
                precise = row.get("CALLER_PRECISE", MISSING)
                if (
                    imprecise not in (MISSING, "", "False")
                    or precise in ("false", "False", "0")
                ):
                    flags.append("IMPRECISE")

                gq = to_float(row.get("CALLER_GQ"))
                if a.min_gq is not None and gq is not None and gq < a.min_gq:
                    flags.append("LOW_GQ")

                if blacklist:
                    chrom = row.get("CHROM", MISSING)
                    start = row.get("START", MISSING)
                    end_pos = row.get("END", MISSING)
                    end_pos = (
                        end_pos if end_pos not in (MISSING, "") else start
                    )
                    if (
                        chrom not in (MISSING, "")
                        and start not in (MISSING, "")
                    ):
                        try:
                            if svtype in {"BND", "TRA"}:
                                c1, p1, c2, p2, _ = breakend(row)
                                hit = any(c not in ("", MISSING) and p is not None and
                                          overlaps_blacklist(blacklist, c, p, p)
                                          for c, p in ((c1, p1), (c2, p2)))
                            else:
                                hit = overlaps_blacklist(blacklist, chrom, start, end_pos)
                            if hit:
                                flags.append("BLACKLIST_REGION")
                        except ValueError:
                            pass

                status = "FAIL" if fail_reasons else "PASS"
                if status == "PASS":
                    n_pass += 1
                    if rescued_cov_var:
                        n_rescued_cov_var += 1

                row["EVIDENCE_STATUS"] = status
                row["EVIDENCE_FAIL_REASONS"] = (
                    ";".join(fail_reasons) if fail_reasons else MISSING
                )
                row["EVIDENCE_FLAGS"] = (
                    ";".join(flags) if flags else MISSING
                )

                writer.writerow(row)

                if status == "PASS":
                    sv_id = row.get("SV_ID", MISSING)
                    if sv_id not in (MISSING, ""):
                        ids_fh.write(sv_id + "\n")

        tmp_tsv.replace(out_tsv)
        tmp_ids.replace(out_ids)
    except Exception:
        tmp_tsv.unlink(missing_ok=True)
        tmp_ids.unlink(missing_ok=True)
        raise

    print(f"[OK] caller={a.caller} total={n_total} pass={n_pass} "
          f"fail={n_total - n_pass} rescued_cov_var={n_rescued_cov_var} "
          f"output_tsv={out_tsv} output_ids={out_ids}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
