#!/usr/bin/env python3
"""Build the integrated SV-gene analysis table.

Biological/data-lineage rules implemented here:
  * the Jasmine merged VCF is the master SV callset;
  * AnnotSV and caller-concordance annotations are joined to that master callset;
  * needLR contributes population-frequency evidence only;
  * needLR is NEVER joined by gene, because different SVs in the same gene can
    have different population frequencies;
  * BNDs and SVs >=10 Mb remain in the master analysis even though needLR does
    not evaluate those classes;
  * candidate-gene / phenotype ranking remains independent of needLR.

One output row is written per (master SV, overlapping gene) pair.
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
NEEDLR_MAX_SUPPORTED_SVLEN = 10_000_000


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
    lowered = {str(k).lower(): v for k, v in row.items()}
    for name in names:
        value = lowered.get(name.lower(), MISSING)
        if value not in ("", MISSING, None):
            return str(value)
    return MISSING


def to_float(value) -> float | None:
    if value in (None, "", MISSING):
        return None
    try:
        return float(str(value).split(",")[0])
    except (TypeError, ValueError):
        return None


def to_int(value) -> int | None:
    x = to_float(value)
    return int(round(x)) if x is not None and math.isfinite(x) else None


def norm_chrom(value: str) -> str:
    value = str(value).strip()
    return value[3:] if value.lower().startswith("chr") else value


def norm_svtype(value: str) -> str:
    value = str(value).upper().strip()
    aliases = {"TRA": "BND", "BREAKEND": "BND", "DELETION": "DEL", "DUPLICATION": "DUP", "INSERTION": "INS", "INVERSION": "INV"}
    return aliases.get(value, value)


def split_genes(value) -> list[str]:
    if value in ("", MISSING, None):
        return []
    return sorted({g.strip().upper() for g in re.split(r"[,;|/]", str(value)) if g.strip() and g.strip() != MISSING})


def read_tsv(path: str | None) -> list[dict]:
    if not path:
        return []
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def read_panel(path: str | None) -> set[str]:
    if not path:
        return set()
    with open_text(path) as fh:
        return {line.strip().upper() for line in fh if line.strip() and not line.startswith("#")}


def read_vcf(path: str) -> list[dict]:
    rows: list[dict] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF record: {line.rstrip()}")
            chrom, pos, sv_id, ref, alt, qual, filt, info_raw = fields[:8]
            info = parse_info(info_raw)
            svtype = first(info, ["SVTYPE"])
            if svtype == MISSING:
                if "[" in alt or "]" in alt:
                    svtype = "BND"
                elif alt.startswith("<") and alt.endswith(">"):
                    svtype = alt[1:-1]
            row = OrderedDict(
                SV_ID=sv_id,
                CHROM=chrom,
                START=pos,
                END=first(info, ["END"]),
                SVTYPE=norm_svtype(svtype),
                SVLEN=first(info, ["SVLEN"]),
                QUAL=qual,
                FILTER=filt,
                REF=ref,
                ALT=alt,
                INFO_RAW=info_raw,
            )
            for key, value in info.items():
                clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", key)
                row[f"INFO_{clean}"] = value
            rows.append(row)
    return rows


def index_by_id(rows: list[dict]) -> dict[str, list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        sv_id = first(row, ["SV_ID", "ID", "AnnotSV_ID", "SV_ID_AnnotSV"])
        if sv_id != MISSING:
            index[sv_id].append(row)
    return index


def genes_from_annot(row: dict) -> list[str]:
    genes: list[str] = []
    for column in ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"]:
        genes.extend(split_genes(first(row, [column])))
    return sorted(set(genes))


def load_ranking(rows: list[dict]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for row in rows:
        gene = first(row, ["gene", "Gene", "GENE", "SYMBOL"])
        if gene != MISSING:
            index[gene.upper()] = row
    return index


def annotsv_coord(row: dict) -> tuple[str, int | None, int | None, str]:
    return (
        norm_chrom(first(row, ["SV_chrom", "CHROM", "Chr", "chrom"])),
        to_int(first(row, ["SV_start", "START", "Start", "POS"])),
        to_int(first(row, ["SV_end", "END", "End"])),
        norm_svtype(first(row, ["SV_type", "SVTYPE", "Type"])),
    )


def master_coord(sv: dict) -> tuple[str, int | None, int | None, str, int | None]:
    start = to_int(sv.get("START"))
    end = to_int(sv.get("END"))
    svlen = to_int(sv.get("SVLEN"))
    if end is None and start is not None and svlen is not None and norm_svtype(sv.get("SVTYPE", MISSING)) not in {"INS", "BND"}:
        end = start + abs(svlen)
    return norm_chrom(sv.get("CHROM", MISSING)), start, end, norm_svtype(sv.get("SVTYPE", MISSING)), svlen


def annotsv_matches_for_sv(sv: dict, by_id: dict[str, list[dict]], all_rows: list[dict], tolerance: int = 10) -> list[dict]:
    direct = by_id.get(sv["SV_ID"], [])
    if direct:
        return direct
    chrom, start, end, svtype, _ = master_coord(sv)
    if start is None:
        return []
    matches = []
    for row in all_rows:
        rchrom, rstart, rend, rtype = annotsv_coord(row)
        if rchrom != chrom or (rtype not in {MISSING, svtype} and svtype != MISSING):
            continue
        if rstart is None or abs(rstart - start) > tolerance:
            continue
        if end is not None and rend is not None and abs(rend - end) > tolerance:
            continue
        matches.append(row)
    return matches


def parse_needlr_row(row: dict) -> dict | None:
    chrom = first(row, ["Chr", "CHROM", "Chrom", "chromosome"])
    start = to_int(first(row, ["Start_Pos", "START", "Start", "POS"])); end = to_int(first(row, ["End_Pos", "END", "End"])); svtype = norm_svtype(first(row, ["SV_Type", "SVTYPE", "Type"])); svlen = to_int(first(row, ["SV_Length", "SVLEN", "Length"])); query_id = first(row, ["Query ID", "Query_ID", "QueryID", "ID", "SV_ID"])
    if chrom == MISSING or start is None or svtype == MISSING:
        return None
    return {
        "row": row,
        "chrom": norm_chrom(chrom),
        "start": start,
        "end": end,
        "svtype": svtype,
        "svlen": svlen,
        "query_id": query_id,
        "af": first(row, ["Allele_Freq_ALL", "AF", "MAX_AF", "AF_MAX", "SV_AF", "AF_1KGP", "1KGP_AF"]),
    }


def load_needlr(rows: list[dict]) -> list[dict]:
    parsed = []
    for row in rows:
        item = parse_needlr_row(row)
        if item is not None:
            parsed.append(item)
    return parsed


def relative_size_difference(a: int | None, b: int | None) -> float | None:
    if a is None or b is None or a == 0 or b == 0:
        return None
    aa, bb = abs(a), abs(b)
    return abs(aa - bb) / max(aa, bb)


def match_needlr(sv: dict, needlr: list[dict], bp_tolerance: int, rel_size_tolerance: float) -> dict:
    chrom, start, end, svtype, svlen = master_coord(sv)
    abs_len = abs(svlen) if svlen is not None else None

    if svtype == "BND":
        return {"status": "NOT_EVALUABLE_BND"}
    if abs_len is not None and abs_len >= NEEDLR_MAX_SUPPORTED_SVLEN:
        return {"status": "NOT_EVALUABLE_GE_10MB"}
    if start is None:
        return {"status": "NO_MASTER_START"}

    candidates = []
    for item in needlr:
        if item["chrom"] != chrom or item["svtype"] != svtype:
            continue
        start_delta = abs(item["start"] - start)
        if start_delta > bp_tolerance:
            continue

        end_delta = None
        if svtype != "INS" and end is not None and item["end"] is not None:
            end_delta = abs(item["end"] - end)
            if end_delta > bp_tolerance:
                continue

        rel_size = relative_size_difference(svlen, item["svlen"])
        if rel_size is not None and rel_size > rel_size_tolerance:
            continue

        score = start_delta + (end_delta or 0)
        if rel_size is not None:
            score += rel_size * bp_tolerance
        candidates.append((score, start_delta, end_delta, rel_size, item))

    if not candidates:
        return {"status": "NO_MATCH"}

    candidates.sort(key=lambda x: x[0])
    best = candidates[0]
    if len(candidates) > 1 and math.isclose(candidates[1][0], best[0], rel_tol=0.0, abs_tol=1e-9):
        return {"status": "AMBIGUOUS_MATCH"}

    item = best[4]
    caution = "SEX_CHROMOSOME_AF_DENOMINATOR" if chrom in {"X", "Y"} and item["af"] not in (MISSING, "0", "0.0") else MISSING
    return {
        "status": "MATCHED",
        "row": item["row"],
        "query_id": item["query_id"],
        "af": item["af"],
        "bp_distance": int(best[1] + (best[2] or 0)),
        "relative_size_difference": best[3] if best[3] is not None else MISSING,
        "match_method": "CHR_TYPE_BREAKPOINT_SIZE",
        "caution": caution,
    }


def load_caller_summary(rows: list[dict]) -> dict[str, dict]:
    out = {}
    for row in rows:
        sv_id = first(row, ["SV_ID", "ID"])
        if sv_id != MISSING:
            out[sv_id] = row
    return out


def caller_evidence_index(paths: list[str]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = defaultdict(list)
    for path in paths:
        for row in read_tsv(path):
            sv_id = first(row, ["SV_ID", "ID"])
            if sv_id != MISSING:
                out[sv_id].append(row)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Create an integrated SV-gene analysis TSV")
    ap.add_argument("--vcf", required=True, help="Master Jasmine merged VCF")
    ap.add_argument("--annotsv", required=True, help="Genome-wide AnnotSV TSV from the same master VCF")
    ap.add_argument("--needlr", help="needLR RESULTS TSV; population AF only")
    ap.add_argument("--ranking", help="Genome-wide candidate ranking TSV")
    ap.add_argument("--panel", help="Candidate gene panel list")
    ap.add_argument("--caller-summary", help="Jasmine caller-support summary TSV")
    ap.add_argument("--caller-tsv", action="append", default=[], help="Filtered per-caller evidence TSV; repeat once per caller")
    ap.add_argument("--needlr-bp-tolerance", type=int, default=1000)
    ap.add_argument("--needlr-rel-size-tolerance", type=float, default=0.5)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    sv_rows = read_vcf(args.vcf)
    annotsv_rows = read_tsv(args.annotsv)
    annotsv_by_id = index_by_id(annotsv_rows)
    needlr = load_needlr(read_tsv(args.needlr))
    ranking = load_ranking(read_tsv(args.ranking))
    panel = read_panel(args.panel)
    caller_summary = load_caller_summary(read_tsv(args.caller_summary))
    caller_by_id = caller_evidence_index(args.caller_tsv)

    n_direct_annot = sum(1 for sv in sv_rows if annotsv_by_id.get(sv["SV_ID"]))
    print(f"[INFO] AnnotSV direct ID join: {n_direct_annot}/{len(sv_rows)}", file=sys.stderr)

    output_rows: list[OrderedDict] = []
    needlr_status_counts: dict[str, int] = defaultdict(int)

    for sv in sv_rows:
        sv_id = sv["SV_ID"]
        annotsv_matches = annotsv_matches_for_sv(sv, annotsv_by_id, annotsv_rows)
        genes = sorted({gene for match in annotsv_matches for gene in genes_from_annot(match)}) or [MISSING]
        annotsv_row = annotsv_matches[0] if annotsv_matches else {}

        cs = caller_summary.get(sv_id, {})
        callers = first(cs, ["CALLERS"])
        caller_count = first(cs, ["CALLER_COUNT", "SUPP"])
        supp_vec = first(cs, ["SUPP_VEC"])
        support_class = first(cs, ["CALLER_SUPPORT_CLASS"])

        exact_evidence = caller_by_id.get(sv_id, [])
        read_support = []
        for match in exact_evidence:
            caller = first(match, ["CALLER"])
            value = first(match, ["CALLER_SUPPORT"])
            if caller != MISSING and value != MISSING:
                read_support.append(f"{caller}:{value}")

        nl = match_needlr(
            sv, needlr,
            bp_tolerance=args.needlr_bp_tolerance,
            rel_size_tolerance=args.needlr_rel_size_tolerance,
        )
        needlr_status_counts[nl["status"]] += 1

        for gene in genes:
            ranking_row = ranking.get(gene, {}) if gene != MISSING else {}
            row = OrderedDict([
                ("SV_ID", sv["SV_ID"]),
                ("CHROM", sv["CHROM"]),
                ("START", sv["START"]),
                ("END", sv["END"]),
                ("SVTYPE", sv["SVTYPE"]),
                ("SVLEN", sv["SVLEN"]),
                ("QUAL", sv["QUAL"]),
                ("FILTER", sv["FILTER"]),
                ("CALLERS", callers),
                ("CALLER_COUNT", caller_count),
                ("CALLER_SUPPORT_CLASS", support_class),
                ("SUPP_VEC", supp_vec),
                ("EXACT_ID_READ_SUPPORT", ";".join(sorted(read_support)) if read_support else MISSING),
                ("GENES", gene),
                ("NEEDLR_STATUS", nl["status"]),
                ("NEEDLR_AF", nl.get("af", MISSING)),
                ("NEEDLR_QUERY_ID", nl.get("query_id", MISSING)),
                ("NEEDLR_MATCH_METHOD", nl.get("match_method", MISSING)),
                ("NEEDLR_BP_DISTANCE", nl.get("bp_distance", MISSING)),
                ("NEEDLR_REL_SIZE_DIFFERENCE", nl.get("relative_size_difference", MISSING)),
                ("NEEDLR_AF_CAUTION", nl.get("caution", MISSING)),
                ("OMIM", first(annotsv_row, ["OMIM", "AnnotSV_OMIM", "AnnotSV_OMIM_evidence"])),
                ("GENCC", first(annotsv_row, ["GENCC", "GenCC", "AnnotSV_GENCC", "AnnotSV_GENCC_evidence"])),
                ("ANNotsv_Gene", first(annotsv_row, ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"])),
                ("ANNotsv_Classification", first(annotsv_row, ["AnnotSV_ranking", "AnnotSV ranking", "AnnotSV_rank", "AnnotSV_Classification"])),
                ("PANEL_STATUS", "PANEL_GENE" if gene in panel else ("UNRESOLVED" if gene == MISSING else "NONPANEL_GENE")),
                ("PHENOTYPE_SCORE", first(ranking_row, ["phenotype_score", "PHENOTYPE_SCORE"])),
                ("CANDIDATE_CLASS", first(ranking_row, ["classification", "CANDIDATE_CLASS"])),
            ])
            for key, value in sv.items():
                if key.startswith("INFO_"):
                    row[key] = value
            output_rows.append(row)

    fixed_columns = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN", "QUAL", "FILTER",
        "CALLERS", "CALLER_COUNT", "CALLER_SUPPORT_CLASS", "SUPP_VEC", "EXACT_ID_READ_SUPPORT",
        "GENES", "NEEDLR_STATUS", "NEEDLR_AF", "NEEDLR_QUERY_ID", "NEEDLR_MATCH_METHOD",
        "NEEDLR_BP_DISTANCE", "NEEDLR_REL_SIZE_DIFFERENCE", "NEEDLR_AF_CAUTION",
        "OMIM", "GENCC", "ANNotsv_Gene", "ANNotsv_Classification", "PANEL_STATUS",
        "PHENOTYPE_SCORE", "CANDIDATE_CLASS",
    ]
    info_columns = sorted({k for row in output_rows for k in row if k.startswith("INFO_")})
    columns = fixed_columns + info_columns

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in output_rows:
            writer.writerow({column: row.get(column, MISSING) for column in columns})

    status_text = ", ".join(f"{k}={v}" for k, v in sorted(needlr_status_counts.items()))
    print(f"[INFO] needLR matching: {status_text}", file=sys.stderr)
    print(f"[OK] integrated_rows={len(output_rows)} info_columns={len(info_columns)} output={out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
