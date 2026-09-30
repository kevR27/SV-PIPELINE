#!/usr/bin/env python3
"""Find conservatively matching SVs across the study cohort.

This is an internal recurrence analysis, not a population allele-frequency
calculation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
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
    return min(a, b) / max(a, b) >= minimum


def reciprocal_overlap(a: SV, b: SV) -> float:
    left = max(a.start, b.start)
    right = min(a.end, b.end)
    overlap = max(0, right - left)
    len_a = max(1, a.end - a.start)
    len_b = max(1, b.end - b.start)
    return min(overlap / len_a, overlap / len_b)


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
        return abs(a.start - b.start) <= insertion_bp and ratio_close(a.size, b.size)

    if a.svtype == "INV":
        return (
            abs(a.start - b.start) <= breakpoint_bp
            and abs(a.end - b.end) <= breakpoint_bp
        )

    if a.svtype in {"DEL", "DUP"}:
        close_breakpoints = (
            abs(a.start - b.start) <= breakpoint_bp
            and abs(a.end - b.end) <= breakpoint_bp
        )
        return close_breakpoints or reciprocal_overlap(a, b) >= 0.8

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-input", action="append", required=True, help="SAMPLE=/path/to/sv_gene_candidates.tsv")
    parser.add_argument("--output", required=True)
    parser.add_argument("--members-output", required=True)
    parser.add_argument("--breakpoint-bp", type=int, default=1000)
    parser.add_argument("--insertion-bp", type=int, default=500)
    args = parser.parse_args()

    variants: list[SV] = []
    for item in args.sample_input:
        sample, path = item.split("=", 1)
        df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).drop_duplicates("SV_ID")
        for _, row in df.iterrows():
            start = integer(row.get("START"))
            end = integer(row.get("END"))
            if start is None:
                continue
            end = end if end is not None else start
            variants.append(SV(
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
            ))

    uf = UnionFind(len(variants))
    for i in range(len(variants)):
        for j in range(i + 1, len(variants)):
            if variants[i].sample == variants[j].sample:
                continue
            if same_sv(variants[i], variants[j], args.breakpoint_bp, args.insertion_bp):
                uf.union(i, j)

    groups: dict[int, list[SV]] = {}
    for i, sv in enumerate(variants):
        groups.setdefault(uf.find(i), []).append(sv)

    summary_rows = []
    member_rows = []
    for cluster_number, members in enumerate(groups.values(), 1):
        samples = sorted({sv.sample for sv in members})
        representative = members[0]
        recurrence = "RECURRENT" if len(samples) > 1 else "UNIQUE_IN_COHORT"

        summary_rows.append({
            "COHORT_SV_ID": f"COHORT_SV_{cluster_number:06d}",
            "SVTYPE": representative.svtype,
            "CHROM": representative.chrom,
            "START": representative.start,
            "END": representative.end,
            "CHR2": representative.chrom2,
            "POS2": representative.pos2 if representative.pos2 is not None else ".",
            "SAMPLE_COUNT": len(samples),
            "SAMPLES": ";".join(samples),
            "COHORT_RECURRENCE": recurrence,
        })

        for sv in members:
            member_rows.append({
                "COHORT_SV_ID": f"COHORT_SV_{cluster_number:06d}",
                "SAMPLE": sv.sample,
                "SV_ID": sv.sv_id,
                "SVTYPE": sv.svtype,
                "CHROM": sv.chrom,
                "START": sv.start,
                "END": sv.end,
                "COHORT_RECURRENCE": recurrence,
            })

    output = Path(args.output)
    members_output = Path(args.members_output)
    output.parent.mkdir(parents=True, exist_ok=True)
    members_output.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(summary_rows).sort_values(
        ["SAMPLE_COUNT", "SVTYPE", "CHROM"], ascending=[False, True, True]
    ).to_csv(output, sep="\t", index=False)
    pd.DataFrame(member_rows).to_csv(members_output, sep="\t", index=False)

    print(f"[OK] cohort_sv_groups={len(summary_rows)} output={output}")


if __name__ == "__main__":
    main()
