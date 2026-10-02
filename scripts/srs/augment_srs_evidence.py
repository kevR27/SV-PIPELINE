#!/usr/bin/env python3
"""Attach SRS-specific RD/BAF and GRIDSS breakpoint evidence to the master table.

The master Manta/DELLY Jasmine callset remains row-defining. CNVpytor and
GRIDSS are independent computational evidence branches. Unmatched CNVpytor and
GRIDSS events are emitted separately so potentially useful events are not
silently discarded.

No confidence score or pathogenicity class is generated here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


MISSING = "."


def read_tsv(path):
    if not path:
        return []
    with open(path, encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def integer(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def numeric(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def known(value):
    return value not in (None, "", MISSING, "NA", "N/A", "nan", "None")


def norm_type(value):
    text = str(value or "").upper()
    if text == "TRA":
        return "BND"
    for svtype in ("DEL","DUP","INV","INS","BND","CNV"):
        if svtype in text:
            return svtype
    return text or MISSING


def interval_overlap_fraction(a_start, a_end, b_start, b_end):
    values = [a_start, a_end, b_start, b_end]
    if any(v is None for v in values):
        return 0.0
    alo, ahi = sorted((a_start, a_end))
    blo, bhi = sorted((b_start, b_end))
    overlap = max(0, min(ahi, bhi) - max(alo, blo) + 1)
    if overlap <= 0:
        return 0.0
    return min(overlap / max(1, ahi - alo + 1), overlap / max(1, bhi - blo + 1))


def cnv_match(master, cnv):
    mtype = norm_type(master.get("SVTYPE"))
    ctype = norm_type(cnv.get("SVTYPE"))
    if mtype not in {"DEL","DUP","CNV"}:
        return None
    if mtype != "CNV" and ctype != mtype:
        return None
    if master.get("CHROM") != cnv.get("CHROM"):
        return None
    ms, me = integer(master.get("START")), integer(master.get("END"))
    cs, ce = integer(cnv.get("START")), integer(cnv.get("END"))
    if None in (ms, me, cs, ce):
        return None
    reciprocal = interval_overlap_fraction(ms, me, cs, ce)
    mlen, clen = abs(me-ms)+1, abs(ce-cs)+1
    tolerance = max(500, int(0.20 * max(mlen, clen)))
    boundary_distance = abs(ms-cs) + abs(me-ce)
    if reciprocal >= 0.5 or (abs(ms-cs) <= tolerance and abs(me-ce) <= tolerance):
        return (-reciprocal, boundary_distance)
    return None


def gridss_match(master, event, tolerance):
    mtype = norm_type(master.get("SVTYPE"))
    chrom = master.get("CHROM")
    start = integer(master.get("START"))
    end = integer(master.get("END"))
    c1, p1 = event.get("CHROM1"), integer(event.get("POS1"))
    c2, p2 = event.get("CHROM2"), integer(event.get("POS2"))
    if start is None or p1 is None:
        return None

    if mtype == "BND":
        chr2 = master.get("CHR2")
        pos2 = integer(master.get("POS2") or master.get("END"))
        if pos2 is None or p2 is None:
            return None
        direct = chrom == c1 and chr2 == c2 and abs(start-p1) <= tolerance and abs(pos2-p2) <= tolerance
        swapped = chrom == c2 and chr2 == c1 and abs(start-p2) <= tolerance and abs(pos2-p1) <= tolerance
        if direct:
            return abs(start-p1) + abs(pos2-p2)
        if swapped:
            return abs(start-p2) + abs(pos2-p1)
        return None

    if mtype == "INS":
        distances = []
        if chrom == c1:
            distances.append(abs(start-p1))
        if p2 is not None and chrom == c2:
            distances.append(abs(start-p2))
        best = min(distances) if distances else None
        return best if best is not None and best <= tolerance else None

    if mtype in {"DEL","DUP","INV","CNV"} and end is not None and p2 is not None and chrom == c1 == c2:
        left, right = sorted((start, end))
        ep1, ep2 = sorted((p1, p2))
        if abs(left-ep1) <= tolerance and abs(right-ep2) <= tolerance:
            return abs(left-ep1) + abs(right-ep2)
    return None


def annotsv_gene_index(rows):
    index = defaultdict(set)
    for row in rows:
        sv_id = row.get("SV_ID") or row.get("ID") or row.get("AnnotSV_ID")
        if not known(sv_id):
            continue
        for key in ("Gene_name","Gene","GENE","Genes","SYMBOL"):
            value = row.get(key)
            if not known(value):
                continue
            for gene in str(value).replace("|",";").replace(",",";").split(";"):
                gene = gene.strip()
                if gene and gene != MISSING:
                    index[sv_id].add(gene)
    return index


def choose_unique(scored):
    if not scored:
        return None, "NO_MATCH"
    scored = sorted(scored, key=lambda x: x[0])
    if len(scored) > 1 and scored[0][0] == scored[1][0]:
        return None, "AMBIGUOUS_MATCH"
    return scored[0][1], "MATCHED"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--integrated", required=True)
    parser.add_argument("--cnvpytor", required=True)
    parser.add_argument("--cnvpytor-annotsv", required=True)
    parser.add_argument("--gridss-events")
    parser.add_argument("--breakpoint-tolerance", type=int, default=500)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cnvpytor-only-output", required=True)
    parser.add_argument("--gridss-only-output", required=True)
    args = parser.parse_args()

    integrated = read_tsv(args.integrated)
    cnvs = read_tsv(args.cnvpytor)
    gridss = read_tsv(args.gridss_events)
    cnv_genes = annotsv_gene_index(read_tsv(args.cnvpytor_annotsv))

    master_by_id = {}
    for row in integrated:
        master_by_id.setdefault(row.get("SV_ID", MISSING), row)

    cnv_matched_ids = set()
    gridss_matched_ids = set()
    extra_columns = [
        "CNVPYTOR_MATCH_STATUS","CNVPYTOR_ID","CNVPYTOR_LEVEL","CNVPYTOR_EVAL1",
        "CNVPYTOR_Q0","CNVPYTOR_PN","CNVPYTOR_DG","CNVPYTOR_BAF_SHIFT",
        "CNVPYTOR_BAF_PVALUE","CNVPYTOR_HET_COUNT","CNVPYTOR_EVIDENCE_FLAGS",
        "GRIDSS_MATCH_STATUS","GRIDSS_EVENT_IDS","GRIDSS_QUAL_MAX","GRIDSS_FILTERS",
        "GRIDSS_AS","GRIDSS_RAS","GRIDSS_SR","GRIDSS_RP","GRIDSS_VF",
        "SRS_EVIDENCE_MODALITIES","SRS_EVIDENCE_MODALITY_COUNT",
    ]

    output_rows = []
    for row in integrated:
        cnv_scored = []
        for cnv in cnvs:
            score = cnv_match(row, cnv)
            if score is not None:
                cnv_scored.append((score, cnv))
        cnv, cnv_status = choose_unique(cnv_scored)
        if cnv:
            cnv_matched_ids.add(cnv["SV_ID"])

        gridss_scored = []
        for event in gridss:
            distance = gridss_match(row, event, args.breakpoint_tolerance)
            if distance is not None:
                gridss_scored.append(((distance,), event))
        event, gridss_status = choose_unique(gridss_scored)
        if event:
            gridss_matched_ids.add(event["GRIDSS_EVENT_ID"])

        row = dict(row)
        row.update({
            "CNVPYTOR_MATCH_STATUS": cnv_status,
            "CNVPYTOR_ID": cnv.get("SV_ID", MISSING) if cnv else MISSING,
            "CNVPYTOR_LEVEL": cnv.get("CNVPYTOR_LEVEL", MISSING) if cnv else MISSING,
            "CNVPYTOR_EVAL1": cnv.get("CNVPYTOR_EVAL1", MISSING) if cnv else MISSING,
            "CNVPYTOR_Q0": cnv.get("CNVPYTOR_Q0", MISSING) if cnv else MISSING,
            "CNVPYTOR_PN": cnv.get("CNVPYTOR_PN", MISSING) if cnv else MISSING,
            "CNVPYTOR_DG": cnv.get("CNVPYTOR_DG", MISSING) if cnv else MISSING,
            "CNVPYTOR_BAF_SHIFT": cnv.get("CNVPYTOR_BAF_SHIFT", MISSING) if cnv else MISSING,
            "CNVPYTOR_BAF_PVALUE": cnv.get("CNVPYTOR_BAF_PVALUE", MISSING) if cnv else MISSING,
            "CNVPYTOR_HET_COUNT": cnv.get("CNVPYTOR_HET_COUNT", MISSING) if cnv else MISSING,
            "CNVPYTOR_EVIDENCE_FLAGS": cnv.get("EVIDENCE_FLAGS", MISSING) if cnv else MISSING,
            "GRIDSS_MATCH_STATUS": gridss_status if args.gridss_events else "NOT_RUN",
            "GRIDSS_EVENT_IDS": event.get("GRIDSS_EVENT_ID", MISSING) if event else MISSING,
            "GRIDSS_QUAL_MAX": event.get("QUAL_MAX", MISSING) if event else MISSING,
            "GRIDSS_FILTERS": event.get("FILTERS", MISSING) if event else MISSING,
            "GRIDSS_AS": event.get("AS", MISSING) if event else MISSING,
            "GRIDSS_RAS": event.get("RAS", MISSING) if event else MISSING,
            "GRIDSS_SR": event.get("SR", MISSING) if event else MISSING,
            "GRIDSS_RP": event.get("RP", MISSING) if event else MISSING,
            "GRIDSS_VF": event.get("VF", MISSING) if event else MISSING,
        })

        modalities = []
        if known(row.get("CALLERS")):
            modalities.append("BREAKPOINT_CALLER_PE_SR")
        if cnv:
            modalities.append("READ_DEPTH_CNV")
            if known(cnv.get("CNVPYTOR_BAF_SHIFT")) and (numeric(cnv.get("CNVPYTOR_HET_COUNT")) or 0) > 0:
                modalities.append("BAF_CONTEXT")
        if event:
            modalities.append("GRIDSS_LOCAL_ASSEMBLY_BREAKPOINT")
        if str(row.get("MELT_MATCH", "")).upper() == "YES":
            modalities.append("MOBILE_ELEMENT_CALL")
        if str(row.get("EXPANSIONHUNTER_LOCUS_OVERLAP", "")).upper() == "YES":
            modalities.append("REPEAT_LOCUS_CONTEXT")
        row["SRS_EVIDENCE_MODALITIES"] = ";".join(modalities) if modalities else MISSING
        row["SRS_EVIDENCE_MODALITY_COUNT"] = str(len(modalities))
        output_rows.append(row)

    base_columns = list(integrated[0]) if integrated else []
    columns = base_columns + [c for c in extra_columns if c not in base_columns]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    cnv_only = []
    for cnv in cnvs:
        if cnv["SV_ID"] in cnv_matched_ids:
            continue
        record = dict(cnv)
        record["MASTER_MATCH_STATUS"] = "NO_GEOMETRIC_MATCH"
        record["ANNOTSV_GENES"] = ";".join(sorted(cnv_genes.get(cnv["SV_ID"], set()))) or MISSING
        cnv_only.append(record)

    cnv_columns = list(cnvs[0]) if cnvs else [
        "SAMPLE","CALLER","SV_ID","CHROM","START","END","SVTYPE","SVLEN"
    ]
    cnv_columns += [x for x in ("MASTER_MATCH_STATUS","ANNOTSV_GENES") if x not in cnv_columns]
    cnv_output = Path(args.cnvpytor_only_output)
    cnv_output.parent.mkdir(parents=True, exist_ok=True)
    with cnv_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=cnv_columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(cnv_only)

    gridss_only = []
    for event in gridss:
        if event["GRIDSS_EVENT_ID"] in gridss_matched_ids:
            continue
        record = dict(event)
        record["MASTER_MATCH_STATUS"] = "NO_GEOMETRIC_MATCH"
        gridss_only.append(record)
    grid_columns = list(gridss[0]) if gridss else [
        "GRIDSS_EVENT_ID","GRIDSS_RECORD_IDS","CHROM1","POS1","CHROM2","POS2"
    ]
    if "MASTER_MATCH_STATUS" not in grid_columns:
        grid_columns.append("MASTER_MATCH_STATUS")
    grid_output = Path(args.gridss_only_output)
    grid_output.parent.mkdir(parents=True, exist_ok=True)
    with grid_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=grid_columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(gridss_only)

    print(
        f"[OK] master_rows={len(output_rows)} cnvpytor_only={len(cnv_only)} "
        f"gridss_only={len(gridss_only)} output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
