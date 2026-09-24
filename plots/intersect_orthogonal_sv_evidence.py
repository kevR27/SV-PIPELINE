#!/usr/bin/env python3
"""Attach orthogonal and phasing evidence to the integrated Jasmine SV-gene table.

The Jasmine master SV universe is never replaced. Straglr and TLDR are
coordinate-aware orthogonal discovery layers; LongPhase is used as SV phasing
evidence. Independent Straglr/TLDR findings without a master-SV match can be
written to a separate table so they are not discarded.
"""

from __future__ import annotations

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from plot_utils import first_existing, normalize_svtype, read_tsv


def parse_args():
    p = argparse.ArgumentParser(description="Integrate orthogonal evidence with master SVs.")
    p.add_argument("--integrated", required=True)
    p.add_argument("--straglr", default=None)
    p.add_argument("--tldr", default=None)
    p.add_argument("--longphase-vcf", default=None)
    p.add_argument("--straglr-breakpoint-tol", type=int, default=500)
    p.add_argument("--tldr-breakpoint-tol", type=int, default=500)
    p.add_argument("--longphase-breakpoint-tol", type=int, default=500)
    p.add_argument("--include-nonpass-tldr", action="store_true")
    p.add_argument("--independent-output", default=None)
    p.add_argument("--output", required=True)
    return p.parse_args()


def to_int(value):
    try:
        return int(float(value))
    except Exception:
        return None


def normalize_chrom(value) -> str:
    text = str(value)
    return text if text.startswith("chr") else "chr" + text


def interval_overlap(a_start, a_end, b_start, b_end) -> int:
    if None in (a_start, a_end, b_start, b_end):
        return 0
    alo, ahi = sorted((a_start, a_end))
    blo, bhi = sorted((b_start, b_end))
    return max(0, min(ahi, bhi) - max(alo, blo) + 1)


def split_ids(value) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if text in {"", ".", "NA", "nan", "None"}:
        return []
    return [x for x in re.split(r"[;,|\s]+", text) if x and x != "."]


def load_straglr(path: str) -> list[dict]:
    df = read_tsv(path)
    chrom = first_existing(df, ["chrom", "CHROM", "Chr"])
    start = first_existing(df, ["start", "START", "POS"])
    end = first_existing(df, ["end", "END"])
    locus = first_existing(df, ["locus"])
    gene = first_existing(df, ["overlapping_genes"])
    copy_number = first_existing(df, ["copy_number"])
    support = first_existing(df, ["supporting_reads"])
    genotype = first_existing(df, ["genotype"])

    if chrom is None or start is None or end is None:
        raise ValueError("Straglr table needs chromosome/start/end columns.")

    rows = []
    for i, row in df.iterrows():
        rows.append(
            {
                "source_index": i,
                "chrom": normalize_chrom(row[chrom]),
                "start": to_int(row[start]),
                "end": to_int(row[end]),
                "label": str(row[locus]) if locus else f"STRAGLR_{i+1}",
                "genes": str(row[gene]) if gene else ".",
                "copy_number": str(row[copy_number]) if copy_number else ".",
                "supporting_reads": str(row[support]) if support else ".",
                "genotype": str(row[genotype]) if genotype else ".",
            }
        )
    return rows


def load_tldr(path: str, include_nonpass: bool = False) -> list[dict]:
    df = read_tsv(path)
    chrom = first_existing(df, ["Chrom", "chrom", "CHROM", "chr", "Chr"])
    start = first_existing(df, ["Start", "start", "START", "pos", "POS"])
    end = first_existing(df, ["End", "end", "END"])
    family = first_existing(df, ["Family", "family", "FAMILY"])
    subfamily = first_existing(df, ["Subfamily", "subfamily", "SUBFAMILY"])
    used_reads = first_existing(df, ["UsedReads", "used_reads", "support", "SUPPORT"])
    span_reads = first_existing(df, ["SpanReads", "span_reads", "spanning_reads"])
    uuid = first_existing(df, ["UUID", "uuid"])
    filt = first_existing(df, ["Filter", "FILTER", "filter"])
    length = first_existing(df, ["LengthIns", "length_ins", "length"])

    if chrom is None or start is None:
        raise ValueError("TLDR table needs chromosome and insertion-position columns.")

    rows = []
    for i, row in df.iterrows():
        filter_value = str(row[filt]) if filt else "."
        if not include_nonpass and filt and filter_value.upper() != "PASS":
            continue

        label_parts = [
            str(row[x])
            for x in [family, subfamily]
            if x and str(row[x]) not in {"", ".", "NA", "nan"}
        ]
        rows.append(
            {
                "source_index": i,
                "chrom": normalize_chrom(row[chrom]),
                "start": to_int(row[start]),
                "end": to_int(row[end]) if end else to_int(row[start]),
                "label": ":".join(label_parts) if label_parts else f"TLDR_{i+1}",
                "uuid": str(row[uuid]) if uuid else ".",
                "used_reads": str(row[used_reads]) if used_reads else ".",
                "span_reads": str(row[span_reads]) if span_reads else ".",
                "filter": filter_value,
                "length": str(row[length]) if length else ".",
            }
        )
    return rows


def open_vcf(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_info(raw: str) -> dict[str, str]:
    out = {}
    for token in str(raw).split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            out[key] = value
        else:
            out[token] = "True"
    return out


def load_longphase(path: str) -> list[dict]:
    records = []
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            info = parse_info(fields[7])
            fmt = fields[8].split(":")
            sample = fields[9].split(":")
            sample_data = dict(zip(fmt, sample))

            svtype = info.get("SVTYPE", "")
            pos = to_int(fields[1])
            end = to_int(info.get("END", fields[1]))
            gt = sample_data.get("GT", ".")
            ps = sample_data.get("PS", ".")
            records.append(
                {
                    "id": fields[2],
                    "chrom": normalize_chrom(fields[0]),
                    "start": pos,
                    "end": end if end is not None else pos,
                    "svtype": svtype.upper(),
                    "gt": gt,
                    "ps": ps,
                    "phased": "YES" if "|" in gt else "NO",
                }
            )
    return records


def chrom_index(records: list[dict]) -> dict[str, list[dict]]:
    out = defaultdict(list)
    for record in records:
        out[record["chrom"]].append(record)
    return out


def longphase_indexes(records: list[dict], bin_size: int):
    by_id = {}
    by_bin = defaultdict(list)
    for record in records:
        if record["id"] not in {"", ".", "NA"}:
            by_id[record["id"]] = record
        if record["start"] is not None:
            key = (record["chrom"], record["svtype"], record["start"] // bin_size)
            by_bin[key].append(record)
    return by_id, by_bin


def find_longphase_match(row, svtype, chrom, start, end, id_columns, by_id, by_bin, tol):
    ids = []
    for col in id_columns:
        if col in row.index:
            ids.extend(split_ids(row[col]))

    for candidate_id in ids:
        if candidate_id in by_id:
            rec = by_id[candidate_id]
            if rec["chrom"] == chrom and (not svtype or rec["svtype"] == svtype):
                return rec, 0

    if svtype not in {"DEL", "INS", "DUP", "INV"} or start is None:
        return None, None

    center = start // tol
    candidates = []
    for b in [center - 1, center, center + 1]:
        candidates.extend(by_bin.get((chrom, svtype, b), []))

    best = None
    best_distance = None
    for rec in candidates:
        if rec["start"] is None:
            continue
        start_dist = abs(start - rec["start"])
        if svtype == "INS":
            distance = start_dist
            valid = start_dist <= tol
        else:
            rec_end = rec["end"] if rec["end"] is not None else rec["start"]
            end_dist = abs(end - rec_end) if end is not None else start_dist
            distance = max(start_dist, end_dist)
            valid = start_dist <= tol and end_dist <= tol

        if valid and (best_distance is None or distance < best_distance):
            best = rec
            best_distance = distance

    return best, best_distance


def main():
    args = parse_args()
    df = read_tsv(args.integrated)

    chrom_col = first_existing(df, ["CHROM", "Chr", "chrom"])
    start_col = first_existing(df, ["START", "POS", "SV_start"])
    end_col = first_existing(df, ["END", "SV_end"])
    type_col = first_existing(df, ["SVTYPE", "SV_type", "Type"])
    id_col = first_existing(df, ["SV_ID", "ID"])
    if None in (chrom_col, start_col, type_col, id_col):
        raise ValueError("Integrated table needs SV_ID, CHROM, START/POS and SVTYPE columns.")

    straglr = load_straglr(args.straglr) if args.straglr else []
    tldr = load_tldr(args.tldr, args.include_nonpass_tldr) if args.tldr else []
    longphase = load_longphase(args.longphase_vcf) if args.longphase_vcf else []

    straglr_by_chrom = chrom_index(straglr)
    tldr_by_chrom = chrom_index(tldr)
    lp_by_id, lp_by_bin = longphase_indexes(longphase, max(args.longphase_breakpoint_tol, 1))

    out = df.copy()
    for column, default in [
        ("STRAGLR_MATCH", "NO"),
        ("STRAGLR_LOCI", "."),
        ("STRAGLR_GENES", "."),
        ("STRAGLR_COPY_NUMBER", "."),
        ("STRAGLR_SUPPORTING_READS", "."),
        ("TLDR_MATCH", "NO"),
        ("TLDR_INSERTIONS", "."),
        ("TLDR_UUID", "."),
        ("TLDR_USED_READS", "."),
        ("TLDR_SPAN_READS", "."),
        ("LONGPHASE_MATCH", "NO"),
        ("LONGPHASE_ID", "."),
        ("LONGPHASE_GT", "."),
        ("LONGPHASE_PS", "."),
        ("LONGPHASE_PHASED", "NO"),
        ("LONGPHASE_MATCH_DISTANCE", "."),
    ]:
        out[column] = default

    svtypes = normalize_svtype(out[type_col])
    id_columns = [
        c for c in ["SV_ID", "INFO_IDLIST", "INFO_IDLIST_EXT", "INFO_INTRASAMPLE_IDLIST"]
        if c in out.columns
    ]

    matched_straglr = set()
    matched_tldr = set()
    cache = {}

    for idx, row in out.iterrows():
        master_id = str(row[id_col])
        if master_id in cache:
            evidence = cache[master_id]
        else:
            chrom = normalize_chrom(row[chrom_col])
            start = to_int(row[start_col])
            end = to_int(row[end_col]) if end_col else start
            svtype = svtypes.loc[idx]
            if end is None:
                end = start

            smatches, sgenes, scn, ssupport = [], [], [], []
            if start is not None:
                for s in straglr_by_chrom.get(chrom, []):
                    if s["start"] is None or s["end"] is None:
                        continue
                    overlap = interval_overlap(start, end, s["start"], s["end"])
                    near_bp = (
                        min(
                            abs(start - s["start"]),
                            abs(start - s["end"]),
                            abs(end - s["start"]),
                            abs(end - s["end"]),
                        )
                        <= args.straglr_breakpoint_tol
                    )
                    if overlap > 0 or (svtype in {"INS", "BND"} and near_bp):
                        matched_straglr.add(s["source_index"])
                        smatches.append(s["label"])
                        if s["genes"] not in {"", ".", "nan"}:
                            sgenes.append(s["genes"])
                        if s["copy_number"] not in {"", ".", "nan"}:
                            scn.append(s["copy_number"])
                        if s["supporting_reads"] not in {"", ".", "nan"}:
                            ssupport.append(s["supporting_reads"])

            tmatches, tuuids, tused, tspan = [], [], [], []
            if start is not None and svtype in {"INS", "BND"}:
                for t in tldr_by_chrom.get(chrom, []):
                    if t["start"] is None:
                        continue
                    distance = min(abs(start - t["start"]), abs(end - t["start"]))
                    if distance <= args.tldr_breakpoint_tol:
                        matched_tldr.add(t["source_index"])
                        tmatches.append(t["label"])
                        if t["uuid"] not in {"", ".", "nan"}:
                            tuuids.append(t["uuid"])
                        if t["used_reads"] not in {"", ".", "nan"}:
                            tused.append(t["used_reads"])
                        if t["span_reads"] not in {"", ".", "nan"}:
                            tspan.append(t["span_reads"])

            lp_rec, lp_dist = (None, None)
            if args.longphase_vcf and start is not None:
                lp_rec, lp_dist = find_longphase_match(
                    row, svtype, chrom, start, end, id_columns,
                    lp_by_id, lp_by_bin, max(args.longphase_breakpoint_tol, 1),
                )

            evidence = {
                "STRAGLR_MATCH": "YES" if smatches else "NO",
                "STRAGLR_LOCI": ";".join(sorted(set(smatches))) if smatches else ".",
                "STRAGLR_GENES": ";".join(sorted(set(sgenes))) if sgenes else ".",
                "STRAGLR_COPY_NUMBER": ";".join(sorted(set(scn))) if scn else ".",
                "STRAGLR_SUPPORTING_READS": ";".join(sorted(set(ssupport))) if ssupport else ".",
                "TLDR_MATCH": "YES" if tmatches else "NO",
                "TLDR_INSERTIONS": ";".join(sorted(set(tmatches))) if tmatches else ".",
                "TLDR_UUID": ";".join(sorted(set(tuuids))) if tuuids else ".",
                "TLDR_USED_READS": ";".join(sorted(set(tused))) if tused else ".",
                "TLDR_SPAN_READS": ";".join(sorted(set(tspan))) if tspan else ".",
                "LONGPHASE_MATCH": "YES" if lp_rec else "NO",
                "LONGPHASE_ID": lp_rec["id"] if lp_rec else ".",
                "LONGPHASE_GT": lp_rec["gt"] if lp_rec else ".",
                "LONGPHASE_PS": lp_rec["ps"] if lp_rec else ".",
                "LONGPHASE_PHASED": lp_rec["phased"] if lp_rec else "NO",
                "LONGPHASE_MATCH_DISTANCE": str(lp_dist) if lp_dist is not None else ".",
            }
            cache[master_id] = evidence

        for key, value in evidence.items():
            out.at[idx, key] = value

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, sep="\t", index=False)

    if args.independent_output:
        independent_rows = []
        for s in straglr:
            if s["source_index"] in matched_straglr:
                continue
            independent_rows.append(
                {
                    "RECORD_CLASS": "STRAGLR_ONLY",
                    "SOURCE_ID": s["label"],
                    "CHROM": s["chrom"],
                    "START": s["start"],
                    "END": s["end"],
                    "SOURCE_TYPE": "TANDEM_REPEAT",
                    "GENES": s["genes"],
                    "SUPPORT": s["supporting_reads"],
                    "DETAILS": f"copy_number={s['copy_number']};genotype={s['genotype']}",
                    "FILTER": ".",
                }
            )

        for t in tldr:
            if t["source_index"] in matched_tldr:
                continue
            independent_rows.append(
                {
                    "RECORD_CLASS": "TLDR_ONLY",
                    "SOURCE_ID": t["uuid"],
                    "CHROM": t["chrom"],
                    "START": t["start"],
                    "END": t["end"],
                    "SOURCE_TYPE": t["label"],
                    "GENES": ".",
                    "SUPPORT": t["used_reads"],
                    "DETAILS": f"span_reads={t['span_reads']};length={t['length']}",
                    "FILTER": t["filter"],
                }
            )

        independent = pd.DataFrame(
            independent_rows,
            columns=[
                "RECORD_CLASS", "SOURCE_ID", "CHROM", "START", "END",
                "SOURCE_TYPE", "GENES", "SUPPORT", "DETAILS", "FILTER",
            ],
        )
        independent_path = Path(args.independent_output)
        independent_path.parent.mkdir(parents=True, exist_ok=True)
        independent.to_csv(independent_path, sep="\t", index=False)

    unique = out.drop_duplicates(id_col)
    print(
        f"[OK] master_SVs={len(unique)} "
        f"straglr_matches={(unique['STRAGLR_MATCH'] == 'YES').sum()} "
        f"tldr_matches={(unique['TLDR_MATCH'] == 'YES').sum()} "
        f"longphase_matches={(unique['LONGPHASE_MATCH'] == 'YES').sum()} "
        f"output={output}"
    )


if __name__ == "__main__":
    main()
