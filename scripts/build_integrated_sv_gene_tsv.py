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
SPATIAL_BIN_SIZE = 1_000_000


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


ACMG_CLASS_LABELS = {
    1: "BENIGN",
    2: "LIKELY_BENIGN",
    3: "VUS",
    4: "LIKELY_PATHOGENIC",
    5: "PATHOGENIC",
}


def acmg_class_label(value) -> str:
    """Normalize AnnotSV ACMG class while retaining the original value separately."""
    if value in (None, "", MISSING):
        return MISSING

    raw = str(value).strip()
    match = re.search(r"(?:class\s*)?([1-5])(?:\.0)?$", raw, re.IGNORECASE)
    if match:
        return ACMG_CLASS_LABELS[int(match.group(1))]

    upper = raw.upper().replace(" ", "_").replace("-", "_")
    aliases = {
        "BENIGN": "BENIGN",
        "LIKELY_BENIGN": "LIKELY_BENIGN",
        "VARIANT_OF_UNKNOWN_SIGNIFICANCE": "VUS",
        "VARIANT_OF_UNCERTAIN_SIGNIFICANCE": "VUS",
        "VUS": "VUS",
        "LIKELY_PATHOGENIC": "LIKELY_PATHOGENIC",
        "PATHOGENIC": "PATHOGENIC",
    }
    return aliases.get(upper, raw)


def dosage_score_label(value) -> str:
    """Translate ClinGen HI/TS score to a readable evidence term."""
    if value in (None, "", MISSING):
        return "NOT_AVAILABLE"

    raw = str(value).strip().split(";")[0].split("|")[0].strip()
    try:
        score = int(float(raw))
    except ValueError:
        return f"UNPARSED_{raw}"

    return {
        3: "SUFFICIENT_EVIDENCE",
        2: "EMERGING_EVIDENCE",
        1: "LITTLE_EVIDENCE",
        0: "NO_EVIDENCE",
        40: "DOSAGE_SENSITIVITY_UNLIKELY",
        30: "AUTOSOMAL_RECESSIVE_GENE",
        -1: "NOT_EVALUATED",
    }.get(score, f"SCORE_{score}")


def dosage_relevance(svtype: str, hi_value, ts_value) -> str:
    """Choose the ClinGen dosage mechanism that matches the SV direction."""
    svtype = normalize_svtype(svtype)
    if svtype == "DEL":
        return "HI_" + dosage_score_label(hi_value)
    if svtype == "DUP":
        return "TS_" + dosage_score_label(ts_value)
    if svtype == "CNV":
        return "NOT_APPLICABLE_CNV_DIRECTION_UNKNOWN"
    return "NOT_APPLICABLE_NON_CNV"


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
    rows: list[dict] = []
    for path in paths:
        for row in read_tsv(path):
            rows.append(row)
            sv_id = first(row, ["SV_ID", "ID"])
            if sv_id != MISSING:
                by_id[sv_id] = row
    return by_id, build_spatial_index(rows)


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


def candidate_coordinates(row: dict):
    start = as_int(first(row, ["START", "Start_Pos", "SV_start", "POS"]))
    end = as_int(first(row, ["END", "End_Pos", "SV_end"]))
    svlen_raw = as_int(first(row, ["SVLEN", "SV_Length", "SV_length"]))
    if start is None:
        return None, None, 0
    span = abs((end - start) if end is not None else 0)
    svlen = abs(svlen_raw) if svlen_raw is not None else span
    return start, end, svlen


def build_spatial_index(rows: list[dict]):
    bins: dict[tuple[str, str], dict[int, list[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    max_len: dict[tuple[str, str], int] = defaultdict(int)

    for row in rows:
        chrom = first(row, ["CHROM", "Chr", "chrom", "chromosome"])
        svtype = normalize_svtype(
            first(row, ["SVTYPE", "SV_Type", "SV_type", "Type"])
        )
        start, end, svlen = candidate_coordinates(row)
        if chrom == MISSING or svtype == MISSING or start is None:
            continue

        key = (chrom, svtype)
        max_len[key] = max(max_len[key], svlen)

        if svtype in {"INS", "BND"} or end is None:
            lo = hi = start
        else:
            lo, hi = sorted((start, end))

        first_bin = max(0, lo // SPATIAL_BIN_SIZE)
        last_bin = max(0, hi // SPATIAL_BIN_SIZE)
        for bin_id in range(first_bin, last_bin + 1):
            bins[key][bin_id].append(row)

    return {"bins": bins, "max_len": max_len}


def spatial_candidates_for_sv(sv: dict, spatial_index) -> list[dict]:
    svtype = normalize_svtype(sv["SVTYPE"])
    key = (sv["CHROM"], svtype)
    key_bins = spatial_index["bins"].get(key)
    if not key_bins:
        return []

    start = as_int(sv.get("START"))
    end = as_int(sv.get("END"))
    if start is None:
        return []

    svlen = abs(
        as_int(sv.get("SVLEN"))
        or ((end - start) if end is not None else 0)
    )
    max_candidate_len = spatial_index["max_len"].get(key, 0)

    # evidence_match_score uses 20% of the larger SV length as breakpoint
    # tolerance.  Use the largest candidate length in this chromosome/type
    # bucket so the positional prefilter cannot exclude a valid match.
    padding = max(500, int(0.20 * max(svlen, max_candidate_len, 1)))

    if svtype in {"INS", "BND"} or end is None:
        lo = max(0, start - padding)
        hi = start + padding
    else:
        left, right = sorted((start, end))
        lo = max(0, left - padding)
        hi = right + padding

    seen: set[int] = set()
    candidates: list[dict] = []
    first_bin = lo // SPATIAL_BIN_SIZE
    last_bin = hi // SPATIAL_BIN_SIZE

    for bin_id in range(first_bin, last_bin + 1):
        for row in key_bins.get(bin_id, []):
            row_id = id(row)
            if row_id in seen:
                continue
            seen.add(row_id)
            candidates.append(row)

    return candidates


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


def caller_evidence_for_sv(sv: dict, by_id, spatial_index):
    ids = split_values(sv.get("INFO_IDLIST", MISSING))
    matches = [by_id[x] for x in ids if x in by_id]
    method = "IDLIST" if matches else MISSING

    if not matches:
        candidates = spatial_candidates_for_sv(sv, spatial_index)
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
    return build_spatial_index(rows)


def match_needlr(sv: dict, needlr_index):
    svtype = normalize_svtype(sv["SVTYPE"])
    svlen = abs(as_int(sv["SVLEN"]) or 0)
    if svtype == "BND":
        return None, "NOT_EVALUABLE_BND", MISSING
    if svlen >= NEEDLR_MAX_SIZE:
        return None, "NOT_EVALUABLE_GE_10MB", MISSING

    candidates = spatial_candidates_for_sv(sv, needlr_index)
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
    caller_by_id, caller_spatial_index = load_caller_evidence(args.caller_tsv)

    annotsv_by_id, annotsv_by_coord = build_annotsv_indexes(annotsv_rows)
    needlr_index = build_needlr_index(needlr_rows)

    output_rows: list[dict] = []
    annotsv_match_count = 0
    needlr_match_count = 0
    n_sv = len(sv_rows)

    print(
        f"[INFO] Spatial matching enabled: bin_size={SPATIAL_BIN_SIZE:,} bp, "
        f"master_SVs={n_sv}, needLR_rows={len(needlr_rows)}",
        file=sys.stderr,
    )

    for sv_number, sv in enumerate(sv_rows, start=1):
        if sv_number == 1 or sv_number % 1000 == 0 or sv_number == n_sv:
            print(
                f"[INFO] Processing master SV {sv_number}/{n_sv}",
                file=sys.stderr,
            )
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
            sv, caller_by_id, caller_spatial_index
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
                    (
                        "OMIM",
                        first(
                            ann_row,
                            [
                                "OMIM_phenotype",
                                "OMIM",
                                "AnnotSV_OMIM",
                                "AnnotSV_OMIM_evidence",
                            ],
                        ),
                    ),
                    (
                        "OMIM_INHERITANCE",
                        first(ann_row, ["OMIM_inheritance", "OMIM inheritance"]),
                    ),
                    (
                        "OMIM_MORBID",
                        first(ann_row, ["OMIM_morbid", "OMIM_morbid_candidate"]),
                    ),
                    (
                        "GENCC",
                        first(
                            ann_row,
                            [
                                "GenCC_classification",
                                "GENCC_classification",
                                "GENCC",
                                "GenCC",
                                "AnnotSV_GENCC",
                                "AnnotSV_GENCC_evidence",
                            ],
                        ),
                    ),
                    (
                        "GENCC_DISEASE",
                        first(ann_row, ["GenCC_disease", "GENCC_disease"]),
                    ),
                    (
                        "GENCC_MOI",
                        first(ann_row, ["GenCC_moi", "GENCC_moi"]),
                    ),
                    (
                        "GENE_DISEASE_EVIDENCE_SCORE",
                        first(
                            ranking_row,
                            [
                                "gene_disease_evidence_score",
                                "GENE_DISEASE_EVIDENCE_SCORE",
                            ],
                        ),
                    ),
                    (
                        "GENE_DISEASE_EVIDENCE_LEVEL",
                        first(
                            ranking_row,
                            [
                                "gene_disease_evidence_level",
                                "GENE_DISEASE_EVIDENCE_LEVEL",
                            ],
                        ),
                    ),
                    (
                        "GENE_DISEASE_EVIDENCE_SOURCE",
                        first(
                            ranking_row,
                            [
                                "gene_disease_evidence_source",
                                "GENE_DISEASE_EVIDENCE_SOURCE",
                            ],
                        ),
                    ),
                    (
                        "GENE_DISEASE_EVIDENCE_CONFLICT",
                        first(
                            ranking_row,
                            [
                                "gene_disease_evidence_conflict",
                                "GENE_DISEASE_EVIDENCE_CONFLICT",
                            ],
                        ),
                    ),
                    ("CLINGEN_HI", first(ann_row, ["HI"])),
                    ("CLINGEN_TS", first(ann_row, ["TS"])),
                    (
                        "DOSAGE_RELEVANCE",
                        dosage_relevance(
                            sv["SVTYPE"],
                            first(ann_row, ["HI"]),
                            first(ann_row, ["TS"]),
                        ),
                    ),
                    (
                        "ANNOTSV_RANKING_SCORE",
                        first(
                            ann_row,
                            ["AnnotSV_ranking_score", "AnnotSV ranking score"],
                        ),
                    ),
                    (
                        "ANNOTSV_RANKING_CRITERIA",
                        first(
                            ann_row,
                            ["AnnotSV_ranking_criteria", "AnnotSV ranking criteria"],
                        ),
                    ),
                    (
                        "ANNOTSV_ACMG_CLASS_RAW",
                        first(ann_row, ["ACMG_class", "ACMG class"]),
                    ),
                    (
                        "ACMG_CNV_CLASS",
                        (
                            acmg_class_label(first(ann_row, ["ACMG_class", "ACMG class"]))
                            if normalize_svtype(sv["SVTYPE"]) in {"DEL", "DUP"}
                            else "NOT_APPLICABLE_NON_CNV"
                        ),
                    ),
                    (
                        "ACMG_CNV_SCORE",
                        (
                            first(
                                ann_row,
                                ["AnnotSV_ranking_score", "AnnotSV ranking score"],
                            )
                            if normalize_svtype(sv["SVTYPE"]) in {"DEL", "DUP"}
                            else MISSING
                        ),
                    ),
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
                                "ACMG_class",
                                "AnnotSV_Classification",
                                "AnnotSV ranking",
                                "AnnotSV_rank",
                                "AnnotSV_ranking_score",
                            ],
                        ),
                    ),
                    (
                        "ANNOTSV_GENERAL_CLASSIFICATION",
                        (
                            acmg_class_label(
                                first(ann_row, ["ACMG_class", "ACMG class"])
                            )
                            if first(ann_row, ["ACMG_class", "ACMG class"]) != MISSING
                            else MISSING
                        ),
                    ),
                    (
                        "ANNOTSV_CLASSIFICATION_SCOPE",
                        (
                            "GAIN_LOSS_CNV_FRAMEWORK"
                            if normalize_svtype(sv["SVTYPE"]) in {"DEL", "DUP"}
                            else "NO_FORMAL_GAIN_LOSS_CLASS_FOR_SVTYPE"
                        ),
                    ),
                    (
                        "PANEL_STATUS",
                        "PANEL_GENE"
                        if gene in panel
                        else ("UNRESOLVED" if gene == MISSING else "NONPANEL_GENE"),
                    ),
                    ("PHENOTYPE_SCORE", first(ranking_row, ["phenotype_score", "PHENOTYPE_SCORE"])),
                    (
                        "SV_EVIDENCE_SCORE",
                        first(ranking_row, ["SV_evidence_score", "SV_EVIDENCE_SCORE"]),
                    ),
                    (
                        "INTEGRATED_DISCOVERY_SCORE",
                        first(
                            ranking_row,
                            [
                                "integrated_discovery_score",
                                "INTEGRATED_DISCOVERY_SCORE",
                            ],
                        ),
                    ),
                    ("CANDIDATE_CLASS", first(ranking_row, ["classification", "CANDIDATE_CLASS"])),
                ]
            )

            for key, value in sv.items():
                if key.startswith("INFO_"):
                    row[key] = value
            output_rows.append(row)

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
        "OMIM_INHERITANCE",
        "OMIM_MORBID",
        "GENCC",
        "GENCC_DISEASE",
        "GENCC_MOI",
        "GENE_DISEASE_EVIDENCE_SCORE",
        "GENE_DISEASE_EVIDENCE_LEVEL",
        "GENE_DISEASE_EVIDENCE_SOURCE",
        "GENE_DISEASE_EVIDENCE_CONFLICT",
        "CLINGEN_HI",
        "CLINGEN_TS",
        "DOSAGE_RELEVANCE",
        "ANNOTSV_RANKING_SCORE",
        "ANNOTSV_RANKING_CRITERIA",
        "ANNOTSV_ACMG_CLASS_RAW",
        "ACMG_CNV_CLASS",
        "ACMG_CNV_SCORE",
        "ANNotsv_Gene",
        "ANNotsv_Classification",
        "ANNOTSV_GENERAL_CLASSIFICATION",
        "ANNOTSV_CLASSIFICATION_SCOPE",
        "PANEL_STATUS",
        "PHENOTYPE_SCORE",
        "SV_EVIDENCE_SCORE",
        "INTEGRATED_DISCOVERY_SCORE",
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
