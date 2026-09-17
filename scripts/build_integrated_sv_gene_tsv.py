#!/usr/bin/env python3
"""Build the master per-SV/per-gene evidence table.

The Jasmine VCF is the source of truth for the integrated SV callset.  Variant
provenance is recovered from Jasmine SUPP_VEC + IDLIST, whose positions follow
the caller VCF order supplied to Jasmine.  This is important because a Jasmine
merged SV ID is not guaranteed to equal the original ID from every caller.

needLR is treated as a population-frequency side annotation.  It is joined to
Jasmine through the original Sniffles2 query ID, not through gene symbol.  A
gene-level needLR join is unsafe because two different SVs can overlap the same
gene but have different population frequencies.

AnnotSV is joined to the Jasmine master callset by input ID when available and
falls back to coordinates/SVTYPE.  The output contains one row per (SV, gene)
pair and retains caller-specific evidence and raw Jasmine INFO fields.

This script does not assign pathogenicity. Candidate ranking columns are copied
from the separate exploratory ranking table when available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
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


def clean_colname(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())


def split_genes(value: object) -> list[str]:
    if value in ("", MISSING, None):
        return []
    genes = {
        g.strip().upper()
        for g in re.split(r"[,;|]", str(value))
        if g.strip() and g.strip() != MISSING
    }
    return sorted(genes)


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
            x.strip().upper()
            for x in fh
            if x.strip() and not x.startswith("#")
        }


def to_int(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def infer_svtype(info: dict[str, str], alt: str) -> str:
    svtype = first(info, ["SVTYPE"])
    if svtype != MISSING:
        return svtype.upper()
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].upper()
    return MISSING


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
            row = OrderedDict(
                [
                    ("SV_ID", sv_id),
                    ("CHROM", chrom),
                    ("START", pos),
                    ("END", first(info, ["END"])),
                    ("SVTYPE", infer_svtype(info, alt)),
                    ("SVLEN", first(info, ["SVLEN"])),
                    ("QUAL", qual),
                    ("FILTER", filt),
                    ("REF", ref),
                    ("ALT", alt),
                    ("INFO_RAW", info_raw),
                ]
            )
            for key, value in info.items():
                row[f"INFO_{clean_colname(key)}"] = value
            rows.append(row)
    return rows


def normalize_caller(value: str) -> str:
    x = value.strip().lower()
    if x in {"sniffles", "sniffles2"}:
        return "Sniffles2"
    if x in {"cutesv", "cute_sv"}:
        return "cuteSV"
    if x == "delly":
        return "delly"
    return value.strip()


def caller_prefix(caller: str) -> str:
    return {
        "Sniffles2": "SNIFFLES",
        "cuteSV": "CUTESV",
        "delly": "DELLY",
    }.get(caller, clean_colname(caller).upper())


def load_caller_evidence(paths: list[str]) -> tuple[dict[str, dict[str, dict]], dict[str, set[str]]]:
    by_caller: dict[str, dict[str, dict]] = defaultdict(dict)
    ids_by_caller: dict[str, set[str]] = defaultdict(set)
    for path in paths:
        for row in read_tsv(path):
            caller = normalize_caller(first(row, ["CALLER", "caller"]))
            sv_id = first(row, ["SV_ID", "ID"])
            if caller == MISSING or sv_id == MISSING:
                continue
            if sv_id in by_caller[caller]:
                raise ValueError(
                    f"Duplicate SV_ID {sv_id!r} in caller evidence for {caller}. "
                    "Caller VCF IDs must be unique before Jasmine."
                )
            by_caller[caller][sv_id] = row
            ids_by_caller[caller].add(sv_id)
    return by_caller, ids_by_caller


def jasmine_provenance(sv: dict, caller_order: list[str]) -> dict[str, str]:
    """Map Jasmine SUPP_VEC/IDLIST back to caller -> original variant ID."""
    supp_vec = sv.get("INFO_SUPP_VEC", MISSING)
    id_list = sv.get("INFO_IDLIST", MISSING)
    if supp_vec in ("", MISSING) or id_list in ("", MISSING):
        return {}

    supporting = [i for i, bit in enumerate(str(supp_vec)) if bit == "1"]
    ids = [x for x in str(id_list).split(",") if x]

    if len(supporting) != len(ids):
        return {}

    out: dict[str, str] = {}
    for idx, sv_id in zip(supporting, ids):
        if idx < len(caller_order):
            out[caller_order[idx]] = sv_id
    return out


def annotsv_genes(row: dict) -> list[str]:
    genes: list[str] = []
    for name in [
        "Gene_name",
        "Gene",
        "GENE",
        "Genes",
        "gene",
        "SYMBOL",
        "AnnotSV_Gene",
    ]:
        genes.extend(split_genes(first(row, [name])))
    return sorted(set(genes))


def norm_chr(value: str) -> str:
    return value.strip()


def annotsv_coords(row: dict) -> tuple[str, int | None, int | None, str]:
    chrom = first(row, ["SV_chrom", "CHROM", "Chrom", "chrom", "Chr"])
    start = to_int(first(row, ["SV_start", "START", "POS", "Start", "Start_Pos"]))
    end = to_int(first(row, ["SV_end", "END", "End", "End_Pos"]))
    svtype = first(row, ["SV_type", "SVTYPE", "Type", "SV_Type"]).upper()
    return norm_chr(chrom), start, end, svtype


def build_annotsv_indexes(rows: list[dict]):
    by_id: dict[str, list[dict]] = defaultdict(list)
    by_coord: dict[tuple, list[dict]] = defaultdict(list)
    by_start: dict[tuple, list[dict]] = defaultdict(list)

    for row in rows:
        sv_id = first(row, ["SV_ID", "ID", "AnnotSV_ID", "Input_ID"])
        if sv_id != MISSING:
            by_id[sv_id].append(row)

        chrom, start, end, svtype = annotsv_coords(row)
        if chrom != MISSING and start is not None and svtype != MISSING:
            by_start[(chrom, start, svtype)].append(row)
            if end is not None:
                by_coord[(chrom, start, end, svtype)].append(row)

    return by_id, by_coord, by_start


def match_annotsv(sv: dict, indexes) -> tuple[list[dict], str]:
    by_id, by_coord, by_start = indexes
    sv_id = sv["SV_ID"]
    if sv_id in by_id:
        return by_id[sv_id], "SV_ID"

    chrom = norm_chr(sv["CHROM"])
    start = to_int(sv["START"])
    end = to_int(sv["END"])
    svtype = str(sv["SVTYPE"]).upper()

    if start is not None and end is not None:
        hit = by_coord.get((chrom, start, end, svtype), [])
        if hit:
            return hit, "COORD_EXACT"
        # A one-base representation shift can occur when tools translate
        # between VCF/BED-style interval conventions.
        for ds in (-1, 0, 1):
            for de in (-1, 0, 1):
                hit = by_coord.get((chrom, start + ds, end + de, svtype), [])
                if hit:
                    return hit, "COORD_PLUSMINUS1"

    if start is not None:
        hit = by_start.get((chrom, start, svtype), [])
        if hit:
            return hit, "START_SVTYPE"

    return [], "NO_MATCH"


def load_ranking(rows: list[dict]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for row in rows:
        gene = first(row, ["gene", "Gene", "GENE", "SYMBOL"])
        if gene != MISSING:
            index[gene.upper()] = row
    return index


def needlr_query_id(row: dict, known_sniffles_ids: set[str]) -> str | None:
    """Recover the original query Sniffles ID from a needLR v4.x TSV row."""
    for col in ["Query ID", "Query_ID", "QueryID", "Query_IDs"]:
        value = first(row, [col])
        if value != MISSING and value in known_sniffles_ids:
            return value

    raw_id = first(row, ["ID", "Query ID", "Query_ID", "QueryID"])
    if raw_id == MISSING:
        return None

    if raw_id in known_sniffles_ids:
        return raw_id

    # needLR merged IDs commonly contain a query component prefixed with
    # 'needLR.' followed by control IDs. Extract that query component first.
    for token in re.split(r"[,;:]", raw_id):
        token = token.strip()
        if token.startswith("needLR."):
            candidate = token[len("needLR.") :]
            if candidate in known_sniffles_ids:
                return candidate

    # Last-resort exact-token recovery for versions that retain the original ID
    # without a needLR prefix inside a composite ID.
    for token in re.split(r"[,;:]", raw_id):
        token = token.strip()
        if token in known_sniffles_ids:
            return token

    return None


def load_needlr(rows: list[dict], known_sniffles_ids: set[str]):
    index: dict[str, dict] = {}
    unresolved = 0
    for row in rows:
        query_id = needlr_query_id(row, known_sniffles_ids)
        if query_id is None:
            unresolved += 1
            continue
        # One output row per query SV is expected in RESULTS.tsv. If a version
        # emits duplicates, keep the first and make the condition visible.
        index.setdefault(query_id, row)
    return index, unresolved


def needlr_af(row: dict) -> str:
    return first(
        row,
        [
            "Allele_Freq_ALL",
            "Allele_Freq_ALL_Control",
            "AF",
            "MAX_AF",
            "AF_MAX",
            "SV_AF",
            "AF_1KGP",
            "1KGP_AF",
        ],
    )


def choose_annotsv_row(matches: list[dict], gene: str) -> dict:
    if not matches:
        return {}
    if gene != MISSING:
        for row in matches:
            if gene in annotsv_genes(row):
                return row
    return matches[0]


def add_caller_columns(row: OrderedDict, provenance: dict[str, str], caller_evidence):
    support_parts: list[str] = []
    for caller in ["Sniffles2", "cuteSV", "delly"]:
        prefix = caller_prefix(caller)
        sv_id = provenance.get(caller, MISSING)
        ev = caller_evidence.get(caller, {}).get(sv_id, {}) if sv_id != MISSING else {}

        row[f"{prefix}_ID"] = sv_id
        row[f"{prefix}_SUPPORT"] = first(ev, ["CALLER_SUPPORT"])
        row[f"{prefix}_GT"] = first(ev, ["CALLER_GT"])
        row[f"{prefix}_GQ"] = first(ev, ["CALLER_GQ"])
        row[f"{prefix}_DP"] = first(ev, ["CALLER_DP"])
        row[f"{prefix}_FILTER"] = first(ev, ["FILTER"])
        row[f"{prefix}_EVIDENCE_STATUS"] = first(ev, ["EVIDENCE_STATUS"])
        row[f"{prefix}_EVIDENCE_FAIL_REASONS"] = first(ev, ["EVIDENCE_FAIL_REASONS"])
        row[f"{prefix}_EVIDENCE_FLAGS"] = first(ev, ["EVIDENCE_FLAGS"])

        supp = row[f"{prefix}_SUPPORT"]
        if supp != MISSING:
            support_parts.append(f"{caller}:{supp}")

    row["SUPPORT"] = ";".join(support_parts) if support_parts else MISSING


def needlr_status(sv: dict, provenance: dict[str, str], needlr_by_query: dict[str, dict]):
    svtype = str(sv["SVTYPE"]).upper()
    svlen = to_int(sv["SVLEN"])
    sniffles_id = provenance.get("Sniffles2")

    if svtype == "BND":
        return "NOT_EVALUATED_BND", {}
    if svlen is not None and abs(svlen) >= NEEDLR_MAX_SIZE:
        return "NOT_EVALUATED_GE10MB", {}
    if not sniffles_id:
        return "NOT_EVALUATED_NO_SNIFFLES_CALL", {}

    row = needlr_by_query.get(sniffles_id)
    if row:
        return "ANNOTATED", row
    return "NO_NEEDLR_MATCH", {}


def main() -> int:
    ap = argparse.ArgumentParser(description="Create integrated SV-gene evidence TSV")
    ap.add_argument("--vcf", required=True, help="Jasmine master SV VCF")
    ap.add_argument("--annotsv", required=True, help="Genome-wide AnnotSV TSV")
    ap.add_argument("--needlr", help="needLR RESULTS TSV generated from the Sniffles2 query VCF")
    ap.add_argument("--ranking", help="Genome-wide candidate ranking TSV")
    ap.add_argument("--panel", help="Candidate gene panel list")
    ap.add_argument("--caller-tsv", action="append", default=[], help="Filtered per-caller evidence TSV; repeat once per caller")
    ap.add_argument(
        "--caller-order",
        default="Sniffles2,cuteSV,delly",
        help="Jasmine input-file order used to decode SUPP_VEC/IDLIST",
    )
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    caller_order = [normalize_caller(x) for x in args.caller_order.split(",") if x.strip()]
    sv_rows = read_vcf(args.vcf)
    annotsv_rows = read_tsv(args.annotsv)
    ranking = load_ranking(read_tsv(args.ranking))
    panel = read_panel(args.panel)
    caller_evidence, ids_by_caller = load_caller_evidence(args.caller_tsv)

    annotsv_indexes = build_annotsv_indexes(annotsv_rows)
    needlr_rows = read_tsv(args.needlr)
    needlr_by_query, unresolved_needlr = load_needlr(
        needlr_rows, ids_by_caller.get("Sniffles2", set())
    )

    if needlr_rows:
        print(
            f"[INFO] needLR query-ID mapping: {len(needlr_by_query)}/{len(needlr_rows)} "
            f"rows resolved to filtered Sniffles2 IDs; unresolved={unresolved_needlr}",
            file=sys.stderr,
        )

    output_rows: list[OrderedDict] = []
    annotsv_matched = 0
    provenance_matched = 0

    for sv in sv_rows:
        provenance = jasmine_provenance(sv, caller_order)
        if provenance:
            provenance_matched += 1

        annotsv_matches, annotsv_match_method = match_annotsv(sv, annotsv_indexes)
        if annotsv_matches:
            annotsv_matched += 1

        genes = sorted(
            {
                gene
                for match in annotsv_matches
                for gene in annotsv_genes(match)
            }
        )
        if not genes:
            genes = [MISSING]

        nl_status, nl_row = needlr_status(sv, provenance, needlr_by_query)
        callers = [c for c in caller_order if c in provenance]
        caller_count = to_int(sv.get("INFO_SUPP")) or len(callers)

        for gene in genes:
            ann = choose_annotsv_row(annotsv_matches, gene)
            rank = ranking.get(gene, {}) if gene != MISSING else {}

            row = OrderedDict(
                [
                    ("SV_ID", sv["SV_ID"]),
                    ("CHROM", sv["CHROM"]),
                    ("START", sv["START"]),
                    ("END", sv["END"]),
                    ("SVTYPE", sv["SVTYPE"]),
                    ("SVLEN", sv["SVLEN"]),
                    ("QUAL", sv["QUAL"]),
                    ("FILTER", sv["FILTER"]),
                    ("CALLERS", ";".join(callers) if callers else MISSING),
                    ("CALLER_COUNT", str(caller_count) if caller_count else MISSING),
                    ("JASMINE_SUPP", sv.get("INFO_SUPP", MISSING)),
                    ("JASMINE_SUPP_VEC", sv.get("INFO_SUPP_VEC", MISSING)),
                    ("JASMINE_IDLIST", sv.get("INFO_IDLIST", MISSING)),
                    ("GENES", gene),
                    ("PANEL_STATUS", "PANEL_GENE" if gene in panel else ("UNRESOLVED" if gene == MISSING else "NONPANEL_GENE")),
                    ("ANNOTSV_MATCH_METHOD", annotsv_match_method),
                    ("NEEDLR_STATUS", nl_status),
                    ("NEEDLR_QUERY_ID", provenance.get("Sniffles2", MISSING)),
                    ("NEEDLR_AF", needlr_af(nl_row) if nl_row else MISSING),
                    ("NEEDLR_CTRL_SUPPORT", first(nl_row, ["Ctrl_support", "Control_support"]) if nl_row else MISSING),
                    ("NEEDLR_POP_FREQ_ALL", first(nl_row, ["Pop_Freq_ALL", "Pop_Freq_ALL_Control"]) if nl_row else MISSING),
                    ("NEEDLR_AF_NOTE", "NONZERO_SEX_CHROMOSOME_AF_HAS_DENOMINATOR_LIMITATION" if sv["CHROM"] in {"chrX", "chrY", "X", "Y"} and nl_row and needlr_af(nl_row) not in {MISSING, "0", "0.0"} else MISSING),
                    ("OMIM", first(ann, ["OMIM", "AnnotSV_OMIM", "AnnotSV_OMIM_evidence"])),
                    ("GENCC", first(ann, ["GENCC", "GenCC", "AnnotSV_GENCC", "AnnotSV_GENCC_evidence"])),
                    ("ANNOTSV_CLASSIFICATION", first(ann, ["AnnotSV ranking", "AnnotSV_rank", "AnnotSV_Classification", "Ranking", "ACMG_class"])),
                    ("PHENOTYPE_SCORE", first(rank, ["phenotype_score", "PHENOTYPE_SCORE"])),
                    ("DISCOVERY_CLASS", first(rank, ["discovery_relevance_class", "classification", "CANDIDATE_CLASS"])),
                    ("DISCOVERY_SCORE", first(rank, ["integrated_discovery_score", "evidence_score"])),
                ]
            )

            add_caller_columns(row, provenance, caller_evidence)

            # Preserve every Jasmine INFO field.
            for key, value in sv.items():
                if key.startswith("INFO_"):
                    row[key] = value

            # Preserve the selected gene-specific AnnotSV row with an explicit
            # prefix so no annotation is silently discarded or confused with
            # master VCF fields.
            for key, value in ann.items():
                row[f"ANNotsv_{clean_colname(str(key))}"] = value if value not in ("", None) else MISSING

            output_rows.append(row)

    n_sv = len(sv_rows)
    if n_sv:
        print(
            f"[INFO] Jasmine provenance decoded: {provenance_matched}/{n_sv} "
            f"({100.0 * provenance_matched / n_sv:.1f}%)",
            file=sys.stderr,
        )
        print(
            f"[INFO] AnnotSV join: {annotsv_matched}/{n_sv} "
            f"({100.0 * annotsv_matched / n_sv:.1f}%)",
            file=sys.stderr,
        )

    if n_sv and provenance_matched / n_sv < 0.95:
        print(
            "[WARN] More than 5% of Jasmine SVs could not be decoded from SUPP_VEC/IDLIST. "
            "Check Jasmine caller order and header sanitation before interpreting caller concordance.",
            file=sys.stderr,
        )

    if annotsv_rows and n_sv and annotsv_matched / n_sv < 0.50:
        print(
            "[WARN] Less than half of Jasmine SVs matched AnnotSV by ID/coordinates. "
            "Inspect AnnotSV coordinate conventions before interpretation.",
            file=sys.stderr,
        )

    # Stable leading columns followed by dynamic caller/INFO/AnnotSV fields.
    fixed = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN", "QUAL", "FILTER",
        "CALLERS", "CALLER_COUNT", "JASMINE_SUPP", "JASMINE_SUPP_VEC", "JASMINE_IDLIST",
        "GENES", "PANEL_STATUS", "ANNOTSV_MATCH_METHOD",
        "NEEDLR_STATUS", "NEEDLR_QUERY_ID", "NEEDLR_AF", "NEEDLR_CTRL_SUPPORT",
        "NEEDLR_POP_FREQ_ALL", "NEEDLR_AF_NOTE",
        "OMIM", "GENCC", "ANNOTSV_CLASSIFICATION", "PHENOTYPE_SCORE",
        "DISCOVERY_CLASS", "DISCOVERY_SCORE",
        "SNIFFLES_ID", "SNIFFLES_SUPPORT", "SNIFFLES_GT", "SNIFFLES_GQ", "SNIFFLES_DP",
        "SNIFFLES_FILTER", "SNIFFLES_EVIDENCE_STATUS", "SNIFFLES_EVIDENCE_FAIL_REASONS", "SNIFFLES_EVIDENCE_FLAGS",
        "CUTESV_ID", "CUTESV_SUPPORT", "CUTESV_GT", "CUTESV_GQ", "CUTESV_DP",
        "CUTESV_FILTER", "CUTESV_EVIDENCE_STATUS", "CUTESV_EVIDENCE_FAIL_REASONS", "CUTESV_EVIDENCE_FLAGS",
        "DELLY_ID", "DELLY_SUPPORT", "DELLY_GT", "DELLY_GQ", "DELLY_DP",
        "DELLY_FILTER", "DELLY_EVIDENCE_STATUS", "DELLY_EVIDENCE_FAIL_REASONS", "DELLY_EVIDENCE_FLAGS",
        "SUPPORT",
    ]

    dynamic = sorted(
        {
            key
            for row in output_rows
            for key in row
            if key not in fixed
        }
    )
    columns = fixed + dynamic

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
            writer.writerow({col: row.get(col, MISSING) for col in columns})

    print(f"[OK] integrated_rows={len(output_rows)} master_svs={n_sv}")
    print(f"[OK] output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
