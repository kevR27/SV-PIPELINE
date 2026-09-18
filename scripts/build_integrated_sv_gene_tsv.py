#!/usr/bin/env python3
"""Build the final per-patient SV/gene evidence table.

Source-of-truth hierarchy
-------------------------
1. Merged master VCF (Jasmine or SURVIVOR): one row-defining SV universe; no SV is removed because
   another annotation tool cannot evaluate it.
2. AnnotSV: functional/gene annotation of that same master VCF.
3. Merge caller-support summary: exact join by the merged SV_ID.
4. Filtered per-caller evidence: linked through merged IDLIST when available,
   with conservative coordinate/type fallback only for evidence display.
5. needLR: optional supplementary ONT population-frequency evidence. needLR is run on a
   Sniffles2-v2.6.2-compatible query, so its rows are matched back to the master
   SV by SV type and coordinates; they are never joined merely because the same
   gene is present.
   Long-read workflows may provide it; short-read workflows intentionally do
   not.  When no needLR table is supplied, the common output schema is retained
   and NEEDLR_STATUS is set to NOT_APPLICABLE_SRS.
6. Monarch/ranking: gene-level phenotype evidence.

One output row is written per (master SV, overlapping gene).  BNDs and SVs
>=10 Mb remain in the table and are explicitly labelled as not evaluable by
needLR rather than being silently lost.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2**31 - 1)

MISSING = "."
NEEDLR_MAX_SIZE = 10_000_000


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict[str, str]:
    info: dict[str, str] = {}
    if not raw or raw == MISSING:
        return info
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = "True"
    return info


def first(row: dict, names: list[str]) -> str:
    for name in names:
        value = row.get(name, MISSING)
        if value not in ("", MISSING, None):
            return str(value)
    lower = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lower.get(name.lower(), MISSING)
        if value not in ("", MISSING, None):
            return str(value)
    return MISSING


def as_int(value):
    if value in (None, "", MISSING):
        return None
    try:
        return int(float(str(value).replace(",", "")))
    except (TypeError, ValueError):
        return None


def as_float(value):
    if value in (None, "", MISSING):
        return None
    try:
        return float(str(value))
    except (TypeError, ValueError):
        return None


def normalize_svtype(value: str) -> str:
    x = str(value or MISSING).upper()
    if "BND" in x or x == "TRA":
        return "BND"
    for svtype in ("DEL", "DUP", "INS", "INV", "CNV"):
        if svtype in x:
            return svtype
    return x


def split_values(value) -> list[str]:
    if value in (None, "", MISSING):
        return []
    return [x.strip() for x in re.split(r"[,;|]", str(value)) if x.strip()]


def split_genes(value) -> list[str]:
    return sorted({x.upper() for x in split_values(value) if x not in {".", "NA"}})


def read_tsv(path: str | None) -> list[dict]:
    if not path:
        return []
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def read_panel(path: str | None) -> set[str]:
    if not path:
        return set()
    with open_text(path) as fh:
        return {
            line.strip().upper()
            for line in fh
            if line.strip() and not line.startswith("#")
        }


def read_vcf(path: str):
    rows: list[dict] = []
    info_keys: OrderedDict[str, None] = OrderedDict()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<"):
                match = re.search(r"ID=([^,>]+)", line)
                if match:
                    info_keys.setdefault(match.group(1), None)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")
            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)
            for key in info:
                info_keys.setdefault(key, None)

            svtype = normalize_svtype(first(info, ["SVTYPE"]))
            if svtype == MISSING:
                if "[" in alt or "]" in alt:
                    svtype = "BND"
                elif alt.startswith("<") and alt.endswith(">"):
                    svtype = normalize_svtype(alt[1:-1])

            end = first(info, ["END"])
            if svtype == "BND":
                # Keep second-breakpoint fields separately; END is not treated
                # as an interval endpoint for BND matching.
                end = first(info, ["POS2", "END"])

            row = OrderedDict(
                [
                    ("SV_ID", sv_id),
                    ("CHROM", chrom),
                    ("START", pos),
                    ("END", end),
                    ("CHR2", first(info, ["CHR2"])),
                    ("POS2", first(info, ["POS2"])),
                    ("SVTYPE", svtype),
                    ("SVLEN", first(info, ["SVLEN"])),
                    ("QUAL", qual),
                    ("FILTER", filt),
                    ("REF", ref),
                    ("ALT", alt),
                    ("INFO_RAW", info_raw),
                ]
            )
            for key, value in info.items():
                clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", key)
                row[f"INFO_{clean}"] = value
            rows.append(row)
    return rows, info_keys


def coord_tuple_from_annotsv(row: dict):
    chrom = first(row, ["SV_chrom", "CHROM", "Chr", "chrom"])
    start = as_int(first(row, ["SV_start", "START", "Start", "POS"]))
    end = as_int(first(row, ["SV_end", "END", "End"]))
    svtype = normalize_svtype(first(row, ["SV_type", "SVTYPE", "Type"]))
    if chrom == MISSING or start is None:
        return None
    return chrom, start, end, svtype


def build_annotsv_indexes(rows: list[dict]):
    by_id: dict[str, list[dict]] = defaultdict(list)
    by_coord: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        sv_id = first(row, ["SV_ID", "ID", "AnnotSV_ID", "AnnotSV ID"])
        if sv_id != MISSING:
            # Some historical needLR-derived AnnotSV files used a prefix.
            if sv_id.startswith("needLR."):
                sv_id = sv_id[len("needLR.") :]
            by_id[sv_id].append(row)
        coord = coord_tuple_from_annotsv(row)
        if coord:
            by_coord[coord].append(row)
    return by_id, by_coord


def annotsv_matches_for_sv(sv: dict, by_id, by_coord):
    exact = by_id.get(sv["SV_ID"], [])
    if exact:
        return exact, "SV_ID"
    key = (
        sv["CHROM"],
        as_int(sv["START"]),
        as_int(sv["END"]),
        normalize_svtype(sv["SVTYPE"]),
    )
    coord = by_coord.get(key, [])
    if coord:
        return coord, "COORDINATE"
    return [], "NO_MATCH"


def genes_from_annot(row: dict) -> list[str]:
    genes: list[str] = []
    for column in [
        "Gene_name",
        "Gene",
        "GENE",
        "Genes",
        "gene",
        "SYMBOL",
        "GeneID",
        "AnnotSV_Gene",
    ]:
        genes.extend(split_genes(first(row, [column])))
    return sorted(set(genes))


def load_ranking(rows: list[dict]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for row in rows:
        gene = first(row, ["gene", "Gene", "GENE", "SYMBOL"])
        if gene != MISSING:
            index[gene.upper()] = row
    return index


def load_caller_summary(rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for row in rows:
        sv_id = first(row, ["SV_ID", "ID"])
        if sv_id != MISSING:
            out[sv_id] = row
    return out


def load_caller_evidence(paths: list[str]):
    by_id: dict[str, dict] = {}
    by_type_chrom: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for path in paths:
        for row in read_tsv(path):
            sv_id = first(row, ["SV_ID", "ID"])
            if sv_id != MISSING:
                by_id[sv_id] = row
            chrom = first(row, ["CHROM", "Chr", "chrom"])
            svtype = normalize_svtype(first(row, ["SVTYPE", "SV_type", "Type"]))
            if chrom != MISSING and svtype != MISSING:
                by_type_chrom[(chrom, svtype)].append(row)
    return by_id, by_type_chrom


def interval_overlap(a_start, a_end, b_start, b_end):
    if None in (a_start, a_end, b_start, b_end):
        return 0.0
    alo, ahi = sorted((a_start, a_end))
    blo, bhi = sorted((b_start, b_end))
    overlap = max(0, min(ahi, bhi) - max(alo, blo) + 1)
    if overlap <= 0:
        return 0.0
    alen = max(1, ahi - alo + 1)
    blen = max(1, bhi - blo + 1)
    return min(overlap / alen, overlap / blen)


def evidence_match_score(sv: dict, row: dict):
    start1 = as_int(sv.get("START"))
    end1 = as_int(sv.get("END"))
    start2 = as_int(first(row, ["START", "Start_Pos", "SV_start", "POS"]))
    end2 = as_int(first(row, ["END", "End_Pos", "SV_end"]))
    svtype = normalize_svtype(sv.get("SVTYPE", MISSING))
    if start1 is None or start2 is None:
        return None

    len1 = abs(as_int(sv.get("SVLEN")) or ((end1 - start1) if end1 is not None else 0))
    len2_raw = as_int(first(row, ["SVLEN", "SV_Length", "SV_length"]))
    len2 = abs(len2_raw or ((end2 - start2) if end2 is not None else 0))
    tolerance = max(500, int(0.20 * max(len1, len2, 1)))

    if svtype == "INS":
        dist = abs(start1 - start2)
        if dist > tolerance:
            return None
        if len1 and len2 and min(len1, len2) / max(len1, len2) < 0.5:
            return None
        return float(dist)

    if svtype == "BND":
        dist = abs(start1 - start2)
        return float(dist) if dist <= tolerance else None

    if end1 is None or end2 is None:
        dist = abs(start1 - start2)
        return float(dist) if dist <= tolerance else None

    reciprocal = interval_overlap(start1, end1, start2, end2)
    bp_dist = abs(start1 - start2) + abs(end1 - end2)
    if reciprocal >= 0.5 or (
        abs(start1 - start2) <= tolerance and abs(end1 - end2) <= tolerance
    ):
        return float(bp_dist)
    return None


def caller_evidence_for_sv(sv: dict, by_id, by_type_chrom):
    ids = split_values(sv.get("INFO_IDLIST", MISSING))
    matches = [by_id[x] for x in ids if x in by_id]
    method = "IDLIST" if matches else MISSING

    if not matches:
        candidates = by_type_chrom.get((sv["CHROM"], normalize_svtype(sv["SVTYPE"])), [])
        scored = []
        for row in candidates:
            score = evidence_match_score(sv, row)
            if score is not None:
                scored.append((score, row))
        # Keep the best coordinate match per caller.
        best_by_caller: dict[str, tuple[float, dict]] = {}
        for score, row in scored:
            caller = first(row, ["CALLER"])
            if caller == MISSING:
                continue
            if caller not in best_by_caller or score < best_by_caller[caller][0]:
                best_by_caller[caller] = (score, row)
        matches = [item[1] for item in best_by_caller.values()]
        if matches:
            method = "COORDINATE_FALLBACK"

    details = []
    flags = []
    for row in sorted(matches, key=lambda r: first(r, ["CALLER"])):
        caller = first(row, ["CALLER"])
        support = first(row, ["CALLER_SUPPORT"])
        status = first(row, ["EVIDENCE_STATUS"])
        details.append(f"{caller}:{support}")
        flag = first(row, ["EVIDENCE_FLAGS"])
        if flag != MISSING:
            flags.append(f"{caller}:{flag}")
    return (
        ";".join(details) if details else MISSING,
        ";".join(flags) if flags else MISSING,
        method,
    )


def build_needlr_index(rows: list[dict]):
    index: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        chrom = first(row, ["Chr", "CHROM", "chrom", "chromosome"])
        svtype = normalize_svtype(first(row, ["SV_Type", "SVTYPE", "SV_type", "Type"]))
        if chrom != MISSING and svtype != MISSING:
            index[(chrom, svtype)].append(row)
    return index


def match_needlr(sv: dict, needlr_index):
    svtype = normalize_svtype(sv["SVTYPE"])
    svlen = abs(as_int(sv["SVLEN"]) or 0)
    if svtype == "BND":
        return None, "NOT_EVALUABLE_BND", MISSING
    if svlen >= NEEDLR_MAX_SIZE:
        return None, "NOT_EVALUABLE_GE_10MB", MISSING

    candidates = needlr_index.get((sv["CHROM"], svtype), [])
    scored = []
    for row in candidates:
        score = evidence_match_score(sv, row)
        if score is not None:
            scored.append((score, row))
    if not scored:
        return None, "NO_MATCH", MISSING
    scored.sort(key=lambda x: x[0])
    best_score, best_row = scored[0]
    return best_row, "MATCHED", str(int(best_score))


def annotsv_row_for_gene(matches: list[dict], gene: str):
    if not matches:
        return {}
    if gene == MISSING:
        return matches[0]
    for row in matches:
        if gene in genes_from_annot(row):
            return row
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Create integrated SV/gene evidence TSV")
    parser.add_argument("--vcf", required=True, help="Merged master SV VCF")
    parser.add_argument("--annotsv", required=True, help="Genome-wide AnnotSV TSV")
    parser.add_argument("--needlr", help="needLR RESULTS TSV")
    parser.add_argument("--ranking", help="Genome-wide candidate ranking TSV")
    parser.add_argument("--panel", help="Candidate gene panel list")
    parser.add_argument("--caller-summary", help="Summary produced from the same merged VCF")
    parser.add_argument(
        "--caller-tsv",
        action="append",
        default=[],
        help="Filtered per-caller evidence TSV; repeat once per caller",
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sv_rows, _ = read_vcf(args.vcf)
    annotsv_rows = read_tsv(args.annotsv)
    needlr_enabled = bool(args.needlr)
    needlr_rows = read_tsv(args.needlr)
    ranking = load_ranking(read_tsv(args.ranking))
    panel = read_panel(args.panel)
    caller_summary = load_caller_summary(read_tsv(args.caller_summary))
    caller_by_id, caller_by_type_chrom = load_caller_evidence(args.caller_tsv)

    annotsv_by_id, annotsv_by_coord = build_annotsv_indexes(annotsv_rows)
    needlr_index = build_needlr_index(needlr_rows)

    output_rows: list[dict] = []
    annotsv_match_count = 0
    needlr_match_count = 0

    for sv in sv_rows:
        sv_id = sv["SV_ID"]
        ann_matches, ann_match_method = annotsv_matches_for_sv(
            sv, annotsv_by_id, annotsv_by_coord
        )
        if ann_matches:
            annotsv_match_count += 1

        genes = sorted(
            {gene for match in ann_matches for gene in genes_from_annot(match)}
        )
        if not genes:
            genes = [MISSING]

        summary = caller_summary.get(sv_id, {})
        callers = first(summary, ["CALLERS"])
        caller_count = first(summary, ["CALLER_COUNT"])
        supp = first(summary, ["SUPP"])
        supp_vec = first(summary, ["SUPP_VEC"])
        support_class = first(summary, ["CALLER_SUPPORT_CLASS"])

        support_detail, evidence_flags, evidence_match_method = caller_evidence_for_sv(
            sv, caller_by_id, caller_by_type_chrom
        )

        if needlr_enabled:
            needlr_row, needlr_status, needlr_distance = match_needlr(
                sv, needlr_index
            )
            if needlr_row:
                needlr_match_count += 1
            else:
                needlr_row = {}
        else:
            needlr_row = {}
            needlr_status = "NOT_APPLICABLE_SRS"
            needlr_distance = MISSING

        for gene in genes:
            ann_row = annotsv_row_for_gene(ann_matches, gene)
            ranking_row = ranking.get(gene, {}) if gene != MISSING else {}

            row = OrderedDict(
                [
                    ("SV_ID", sv["SV_ID"]),
                    ("CHROM", sv["CHROM"]),
                    ("START", sv["START"]),
                    ("END", sv["END"]),
                    ("CHR2", sv["CHR2"]),
                    ("POS2", sv["POS2"]),
                    ("SVTYPE", sv["SVTYPE"]),
                    ("SVLEN", sv["SVLEN"]),
                    ("QUAL", sv["QUAL"]),
                    ("FILTER", sv["FILTER"]),
                    ("CALLERS", callers),
                    ("CALLER_COUNT", caller_count),
                    ("CALLER_SUPPORT_CLASS", support_class),
                    ("SUPP", supp),
                    ("SUPP_VEC", supp_vec),
                    ("CALLER_READ_SUPPORT", support_detail),
                    ("CALLER_EVIDENCE_FLAGS", evidence_flags),
                    ("CALLER_EVIDENCE_MATCH", evidence_match_method),
                    ("GENES", gene),
                    ("ANNOTSV_MATCH", ann_match_method),
                    (
                        "NEEDLR_AF",
                        first(
                            needlr_row,
                            [
                                "Allele_Freq_ALL",
                                "AlleleFreqAll",
                                "AF",
                                "MAX_AF",
                                "AF_MAX",
                                "SV_AF",
                                "AF_1KGP",
                                "1KGP_AF",
                            ],
                        ),
                    ),
                    ("NEEDLR_STATUS", needlr_status),
                    ("NEEDLR_MATCH_DISTANCE", needlr_distance),
                    ("NEEDLR_GENES", first(needlr_row, ["Genes", "Gene", "GENE"])),
                    ("NEEDLR_HPO", first(needlr_row, ["HPO", "HPO_terms", "HPO_Terms"])),
                    ("NEEDLR_OMIM", first(needlr_row, ["OMIM", "OMIM_phenotypes"])),
                    ("NEEDLR_GENCC", first(needlr_row, ["GenCC", "GENCC", "GenCC_phenotypes"])),
                    ("OMIM", first(ann_row, ["OMIM", "AnnotSV_OMIM", "AnnotSV_OMIM_evidence"])),
                    ("GENCC", first(ann_row, ["GENCC", "GenCC", "AnnotSV_GENCC", "AnnotSV_GENCC_evidence"])),
                    (
                        "ANNotsv_Gene",
                        first(
                            ann_row,
                            ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"],
                        ),
                    ),
                    (
                        "ANNotsv_Classification",
                        first(
                            ann_row,
                            [
                                "AnnotSV ranking",
                                "AnnotSV_rank",
                                "AnnotSV_ranking_score",
                                "AnnotSV_Classification",
                                "ACMG_class",
                            ],
                        ),
                    ),
                    (
                        "PANEL_STATUS",
                        "PANEL_GENE"
                        if gene in panel
                        else ("UNRESOLVED" if gene == MISSING else "NONPANEL_GENE"),
                    ),
                    ("PHENOTYPE_SCORE", first(ranking_row, ["phenotype_score", "PHENOTYPE_SCORE"])),
                    ("CANDIDATE_CLASS", first(ranking_row, ["classification", "CANDIDATE_CLASS"])),
                ]
            )

            for key, value in sv.items():
                if key.startswith("INFO_"):
                    row[key] = value
            output_rows.append(row)

    n_sv = len(sv_rows)
    ann_rate = (annotsv_match_count / n_sv * 100.0) if n_sv else 0.0
    print(
        f"[INFO] AnnotSV master-SV match: {annotsv_match_count}/{n_sv} "
        f"({ann_rate:.1f}%)",
        file=sys.stderr,
    )
    if annotsv_rows and n_sv and ann_rate < 80.0:
        print(
            "[WARN] AnnotSV matching is below 80%. AnnotSV is run directly on "
            "the merged master VCF, so inspect AnnotSV_ID/SV_chrom/SV_start/SV_end "
            "columns before interpreting apparently unannotated calls.",
            file=sys.stderr,
        )

    if needlr_enabled:
        eligible_needlr = sum(
            1
            for sv in sv_rows
            if normalize_svtype(sv["SVTYPE"]) != "BND"
            and abs(as_int(sv["SVLEN"]) or 0) < NEEDLR_MAX_SIZE
        )
        print(
            f"[INFO] needLR matches: {needlr_match_count}/{eligible_needlr} "
            "needLR-eligible master SVs",
            file=sys.stderr,
        )
    else:
        print(
            "[INFO] needLR not supplied; marked NOT_APPLICABLE_SRS.",
            file=sys.stderr,
        )

    fixed_columns = [
        "SV_ID",
        "CHROM",
        "START",
        "END",
        "CHR2",
        "POS2",
        "SVTYPE",
        "SVLEN",
        "QUAL",
        "FILTER",
        "CALLERS",
        "CALLER_COUNT",
        "CALLER_SUPPORT_CLASS",
        "SUPP",
        "SUPP_VEC",
        "CALLER_READ_SUPPORT",
        "CALLER_EVIDENCE_FLAGS",
        "CALLER_EVIDENCE_MATCH",
        "GENES",
        "ANNOTSV_MATCH",
        "NEEDLR_AF",
        "NEEDLR_STATUS",
        "NEEDLR_MATCH_DISTANCE",
        "NEEDLR_GENES",
        "NEEDLR_HPO",
        "NEEDLR_OMIM",
        "NEEDLR_GENCC",
        "OMIM",
        "GENCC",
        "ANNotsv_Gene",
        "ANNotsv_Classification",
        "PANEL_STATUS",
        "PHENOTYPE_SCORE",
        "CANDIDATE_CLASS",
    ]
    info_columns = sorted(
        {key for row in output_rows for key in row if key.startswith("INFO_")}
    )
    columns = fixed_columns + info_columns

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in output_rows:
            writer.writerow({column: row.get(column, MISSING) for column in columns})

    print(f"[OK] integrated_rows={len(output_rows)} output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
