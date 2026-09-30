#!/usr/bin/env python3
"""Find conservatively matching SVs across the study cohort.

This is an internal recurrence analysis, not a population allele-frequency
calculation. Matching is intentionally strict so recurrent technical artifacts
can be recognized without pretending that loosely similar events are the same
allele.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import pandas as pd


@dataclass
class SV:
    sample: str
    sv_id: str
    svtype: str
    chrom: str
    start: int
    end: int
    chrom2: str
    pos2: int | None
    size: int | None
    orientation: str


def integer(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def ratio_close(a: int | None, b: int | None, minimum: float = 0.5) -> bool:
    if not a or not b:
        return True
    return min(abs(a), abs(b)) / max(abs(a), abs(b)) >= minimum


def same_sv(a: SV, b: SV, breakpoint_bp: int, insertion_bp: int) -> bool:
    if a.svtype != b.svtype:
        return False

    if a.svtype in {"BND", "TRA"}:
        if a.pos2 is None or b.pos2 is None:
            return False

        direct = (
            a.chrom == b.chrom
            and a.chrom2 == b.chrom2
            and abs(a.start - b.start) <= breakpoint_bp
            and abs(a.pos2 - b.pos2) <= breakpoint_bp
        )
        swapped = (
            a.chrom == b.chrom2
            and a.chrom2 == b.chrom
            and abs(a.start - b.pos2) <= breakpoint_bp
            and abs(a.pos2 - b.start) <= breakpoint_bp
        )
        if not (direct or swapped):
            return False

        if a.orientation not in {"", "."} and b.orientation not in {"", "."}:
            return a.orientation == b.orientation
        return True

    if a.chrom != b.chrom:
        return False

    if a.svtype == "INS":
        return (
            abs(a.start - b.start) <= insertion_bp
            and ratio_close(a.size, b.size)
        )

    if a.svtype in {"DEL", "DUP", "INV"}:
        return (
            abs(a.start - b.start) <= breakpoint_bp
            and abs(a.end - b.end) <= breakpoint_bp
        )

    return (
        abs(a.start - b.start) <= breakpoint_bp
        and ratio_close(a.size, b.size)
    )


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def nearby_keys(sv: SV, breakpoint_bp: int, insertion_bp: int):
    """Generate only nearby coordinate bins, avoiding all-by-all comparison."""
    if sv.svtype == "INS":
        step = max(insertion_bp, 1)
        base = sv.start // step
        for shift in (-1, 0, 1):
            yield (sv.svtype, sv.chrom, base + shift)
        return

    step = max(breakpoint_bp, 1)

    if sv.svtype in {"BND", "TRA"} and sv.pos2 is not None:
        start_bin = sv.start // step
        pos2_bin = sv.pos2 // step
        for left_shift, right_shift in product((-1, 0, 1), repeat=2):
            yield (
                sv.svtype,
                sv.chrom,
                sv.chrom2,
                start_bin + left_shift,
                pos2_bin + right_shift,
            )
            yield (
                sv.svtype,
                sv.chrom2,
                sv.chrom,
                pos2_bin + right_shift,
                start_bin + left_shift,
            )
        return

    start_bin = sv.start // step
    end_bin = sv.end // step
    for left_shift, right_shift in product((-1, 0, 1), repeat=2):
        yield (
            sv.svtype,
            sv.chrom,
            start_bin + left_shift,
            end_bin + right_shift,
        )


def primary_key(sv: SV, breakpoint_bp: int, insertion_bp: int):
    if sv.svtype == "INS":
        step = max(insertion_bp, 1)
        return (sv.svtype, sv.chrom, sv.start // step)

    step = max(breakpoint_bp, 1)
    if sv.svtype in {"BND", "TRA"} and sv.pos2 is not None:
        return (
            sv.svtype,
            sv.chrom,
            sv.chrom2,
            sv.start // step,
            sv.pos2 // step,
        )

    return (
        sv.svtype,
        sv.chrom,
        sv.start // step,
        sv.end // step,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-input",
        action="append",
        required=True,
        help="SAMPLE=/path/to/sv_gene_candidates.tsv",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--members-output", required=True)
    parser.add_argument("--breakpoint-bp", type=int, default=1000)
    parser.add_argument("--insertion-bp", type=int, default=500)
    args = parser.parse_args()

    variants: list[SV] = []

    for item in args.sample_input:
        sample, path = item.split("=", 1)
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        df = df.drop_duplicates("SV_ID")

        for _, row in df.iterrows():
            start = integer(row.get("START"))
            end = integer(row.get("END"))
            if start is None:
                continue
            end = end if end is not None else start

            variants.append(
                SV(
                    sample=sample,
                    sv_id=str(row.get("SV_ID", ".")),
                    svtype=str(row.get("SVTYPE", ".")).upper(),
                    chrom=str(row.get("CHROM", ".")),
                    start=start,
                    end=end,
                    chrom2=str(row.get("CHR2", ".")),
                    pos2=integer(row.get("POS2")),
                    size=integer(row.get("SVLEN")),
                    orientation=str(row.get("BND_ORIENTATION", ".")),
                )
            )

    uf = UnionFind(len(variants))
    index: dict[tuple, list[int]] = {}

    for i, sv in enumerate(variants):
        candidate_indices = set()

        for key in nearby_keys(sv, args.breakpoint_bp, args.insertion_bp):
            candidate_indices.update(index.get(key, []))

        for j in candidate_indices:
            other = variants[j]
            if other.sample == sv.sample:
                continue
            if same_sv(other, sv, args.breakpoint_bp, args.insertion_bp):
                uf.union(i, j)

        key = primary_key(sv, args.breakpoint_bp, args.insertion_bp)
        index.setdefault(key, []).append(i)

    groups: dict[int, list[SV]] = {}
    for i, sv in enumerate(variants):
        groups.setdefault(uf.find(i), []).append(sv)

    summary_rows = []
    member_rows = []

    for cluster_number, members in enumerate(groups.values(), 1):
        samples = sorted({sv.sample for sv in members})
        representative = members[0]
        recurrence = "RECURRENT" if len(samples) > 1 else "UNIQUE_IN_COHORT"
        cohort_id = f"COHORT_SV_{cluster_number:06d}"

        summary_rows.append(
            {
                "COHORT_SV_ID": cohort_id,
                "SVTYPE": representative.svtype,
                "CHROM": representative.chrom,
                "START": representative.start,
                "END": representative.end,
                "CHR2": representative.chrom2,
                "POS2": representative.pos2 if representative.pos2 is not None else ".",
                "SAMPLE_COUNT": len(samples),
                "SAMPLES": ";".join(samples),
                "COHORT_RECURRENCE": recurrence,
            }
        )

        for sv in members:
            member_rows.append(
                {
                    "COHORT_SV_ID": cohort_id,
                    "SAMPLE": sv.sample,
                    "SV_ID": sv.sv_id,
                    "SVTYPE": sv.svtype,
                    "CHROM": sv.chrom,
                    "START": sv.start,
                    "END": sv.end,
                    "COHORT_RECURRENCE": recurrence,
                }
            )

    output = Path(args.output)
    members_output = Path(args.members_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    members_output.parent.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame(
        summary_rows,
        columns=[
            "COHORT_SV_ID",
            "SVTYPE",
            "CHROM",
            "START",
            "END",
            "CHR2",
            "POS2",
            "SAMPLE_COUNT",
            "SAMPLES",
            "COHORT_RECURRENCE",
        ],
    )
    members = pd.DataFrame(
        member_rows,
        columns=[
            "COHORT_SV_ID",
            "SAMPLE",
            "SV_ID",
            "SVTYPE",
            "CHROM",
            "START",
            "END",
            "COHORT_RECURRENCE",
        ],
    )

    if not summary.empty:
        summary = summary.sort_values(
            ["SAMPLE_COUNT", "SVTYPE", "CHROM"],
            ascending=[False, True, True],
        )

    summary.to_csv(output, sep="\t", index=False)
    members.to_csv(members_output, sep="\t", index=False)

    print(f"[OK] cohort_sv_groups={len(summary)} output={output}")


if __name__ == "__main__":
    main()
