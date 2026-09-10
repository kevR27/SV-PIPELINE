#!/usr/bin/env python3
"""
build_integrated_sv_gene_tsv.py

Build a single integrated SV-gene analysis table by joining together:
  - the merged/filtered SV VCF (source of truth for SV_ID, coordinates, SVTYPE)
  - the genome-wide AnnotSV TSV (gene overlap, OMIM, GenCC, AnnotSV ranking)
  - the needLR TSV (ONT-native population allele frequency, HPO annotation)
  - the genome-wide candidate ranking TSV (phenotype score, candidate class)
  - the per-caller evidence TSVs (which callers/how much support backs each SV)
  - the candidate gene panel list (to flag PANEL_GENE vs NONPANEL_GENE)

One output row is written per (SV, gene) pair, so an SV overlapping several
genes appears on multiple rows. All joins are done by SV_ID (VCF <-> AnnotSV
<-> per-caller evidence) or by gene symbol (AnnotSV gene <-> needLR <-> ranking
<-> panel), using a small set of known column-name aliases per source, since
different tool versions name the same column slightly differently.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

# Placeholder used everywhere a value is absent, to keep the output table
# free of blank cells.
MISSING = "."

def open_text(path: str):
    """Open a plain-text or gzip-compressed file for reading, as text."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(raw: str) -> dict:
    """
    Parse a VCF INFO field string (e.g. "SVTYPE=DEL;SVLEN=-500;IMPRECISE")
    into a dict. Flag-style fields with no "=" are recorded as 'True'.
    """
    info = {}
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
    """
    Return the value of the first column name (in order) present in `row`
    with a non-empty value. Tries an exact match first, then falls back to
    a case-insensitive match — different tool versions capitalize column
    names inconsistently (e.g. "Gene_name" vs "gene_name").
    """
    for name in names:
        value = row.get(name, MISSING)
        if value not in ("", MISSING, None):
            return str(value)

    lowercase_row = {key.lower(): value for key, value in row.items()}
    for name in names:
        value = lowercase_row.get(name.lower(), MISSING)
        if value not in ("", MISSING, None):
            return str(value)

    return MISSING


def split_genes(value) -> list[str]:
    """Split a delimiter-separated gene list (",", ";", or "|") into a
    sorted, deduplicated, uppercased list of gene symbols."""
    if value in ("", MISSING, None):
        return []
    genes = {gene.strip().upper() for gene in re.split(r"[,;|]", str(value)) if gene.strip()}
    return sorted(genes)


def read_tsv(path: str) -> list[dict]:
    """Read a TSV into a list of row dicts. Returns [] if path is falsy."""
    if not path:
        return []
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_panel(path: str) -> set[str]:
    """Read a one-gene-per-line panel list into an uppercased set."""
    if not path:
        return set()
    with open_text(path) as fh:
        return {line.strip().upper() for line in fh if line.strip() and not line.startswith("#")}


def read_vcf(path: str):
    """
    Read a VCF into a list of row dicts (one per record) plus an ordered
    set of INFO keys seen across the file. Each row carries the core VCF
    fields plus one INFO_<key> column per INFO field found in that record,
    so no INFO annotation is lost even if it isn't one of the "known"
    fields referenced elsewhere in this script.

    # for input isn't it okay to use that of already parsed and add missing info?
    """
    rows = []
    info_keys = OrderedDict()

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

            svtype = first(info, ["SVTYPE"])
            if svtype == MISSING:
                if "[" in alt or "]" in alt:
                    svtype = "BND"
                elif alt.startswith("<") and alt.endswith(">"):
                    svtype = alt[1:-1]
                else:
                    svtype = MISSING

            row = OrderedDict([
                ("SV_ID", sv_id),
                ("CHROM", chrom),
                ("START", pos),
                ("END", first(info, ["END"])),
                ("SVTYPE", svtype),
                ("SVLEN", first(info, ["SVLEN"])),
                ("QUAL", qual),
                ("FILTER", filt),
                ("REF", ref),
                ("ALT", alt),
                ("INFO_RAW", info_raw),
            ])

            for key, value in info.items():
                clean_key = re.sub(r"[^A-Za-z0-9_.-]+", "_", key)
                row[f"INFO_{clean_key}"] = value

            row["CALLERS"] = first(info, ["CALLERS", "CALLER", "SOURCE", "SOURCES"])
            row["SUPPORT"] = first(info, ["SUPPORT", "SUPP", "RE", "SU"])

            rows.append(row)

    return rows, info_keys


def index_by_id(rows: list[dict]) -> dict:
    """Index a list of row dicts by SV_ID (trying SV_ID/ID/AnnotSV_ID),
    returning {sv_id: [matching rows]}."""
    index = defaultdict(list)
    for row in rows:
        sv_id = first(row, ["SV_ID", "ID", "AnnotSV_ID"])
        if sv_id != MISSING:
            index[sv_id].append(row)
    return index


def genes_from_annot(row: dict) -> list[str]:
    """Collect every gene symbol out of an AnnotSV row, trying each of the
    known gene-column aliases used across AnnotSV versions."""
    gene_columns = ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "GeneID", "AnnotSV_Gene"]
    genes = []
    for column in gene_columns:
        genes.extend(split_genes(first(row, [column])))
    return sorted(set(genes))


def load_ranking(rows: list[dict]) -> dict:
    """Index the candidate-ranking TSV by uppercased gene symbol."""
    index = {}
    for row in rows:
        gene = first(row, ["gene", "Gene", "GENE", "SYMBOL"])
        if gene != MISSING:
            index[gene.upper()] = row
    return index


def load_needlr(rows: list[dict]) -> dict:
    """Index the needLR TSV by uppercased gene symbol. A needLR row can
    list multiple genes, so it may be indexed under more than one key."""
    index = defaultdict(list)
    for row in rows:
        gene_field = first(row, ["Gene", "gene", "Gene_name", "GENE", "SYMBOL"])
        for gene in split_genes(gene_field):
            index[gene].append(row)
    return index


def main() -> int:
    parser = argparse.ArgumentParser(description="Create integrated SV-gene analysis TSV")
    parser.add_argument("--vcf", required=True, help="Merged/filtered SV VCF")
    parser.add_argument("--annotsv", required=True, help="Genome-wide AnnotSV TSV")
    parser.add_argument("--needlr", help="needLR RESULTS TSV")
    parser.add_argument("--ranking", help="Genome-wide candidate ranking TSV")
    parser.add_argument("--panel", help="Candidate gene panel list")
    parser.add_argument("--caller-tsv", action="append", default=[],
                         help="Per-caller evidence TSV; repeat once per caller")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sv_rows, _ = read_vcf(args.vcf)
    annotsv_rows = read_tsv(args.annotsv)
    needlr_rows = read_tsv(args.needlr)
    ranking = load_ranking(read_tsv(args.ranking))
    panel = read_panel(args.panel)
    annotsv_by_id = index_by_id(annotsv_rows)
    needlr_by_gene = load_needlr(needlr_rows)

    # Log the AnnotSV join quality up front. A silently-broken SV_ID join
    # (e.g. because AnnotSV's ID column format doesn't match the VCF's)
    # would otherwise just leave OMIM/GENCC/ranking columns looking
    # "genuinely empty" rather than "failed to match" — this makes that
    # failure mode visible instead of silent.
    n_sv = len(sv_rows)
    n_matched = sum(1 for sv in sv_rows if annotsv_by_id.get(sv["SV_ID"]))
    match_rate = (n_matched / n_sv * 100) if n_sv else 0.0
    print(f"[INFO] AnnotSV ID join: {n_matched}/{n_sv} SVs matched ({match_rate:.1f}%)", file=sys.stderr)
    if annotsv_rows and n_sv and match_rate < 50:
        print(
            "[WARN] Less than half of the SVs matched an AnnotSV row by ID. "
            "This usually means the AnnotSV ID column format does not match "
            "SV_ID from the VCF (check the AnnotSV/AnnotSV_ID column in "
            "--annotsv) — downstream OMIM/GENCC/ranking fields may be "
            "wrongly empty rather than genuinely absent.",
            file=sys.stderr,
        )

    # Index each caller's evidence TSV by SV_ID so we can report, per SV,
    # which caller(s) called it and with how much support.
    caller_by_id = defaultdict(list)
    for path in args.caller_tsv:
        for row in read_tsv(path):
            sv_id = first(row, ["SV_ID", "ID"])
            if sv_id != MISSING:
                caller_by_id[sv_id].append(row)

    output_rows = []

    for sv in sv_rows:
        sv_id = sv["SV_ID"]

        annotsv_matches = annotsv_by_id.get(sv_id, [])
        genes = sorted({gene for match in annotsv_matches for gene in genes_from_annot(match)})
        if not genes:
            genes = [MISSING]

        caller_matches = caller_by_id.get(sv_id, [])
        caller_names = sorted({match.get("CALLER", "") for match in caller_matches if match.get("CALLER", "")})
        callers = ";".join(caller_names) if caller_names else sv["CALLERS"]

        support_parts = []
        for match in caller_matches:
            support_value = match.get("CALLER_SUPPORT", MISSING)
            if support_value not in ("", MISSING):
                support_parts.append(f"{match.get('CALLER', 'UNKNOWN')}:{support_value}")
        support = ";".join(support_parts) if support_parts else sv["SUPPORT"]

        # One output row per gene the SV overlaps (or a single row with
        # GENES=MISSING if it overlaps none).
        for gene in genes:
            annotsv_row = annotsv_matches[0] if annotsv_matches else {}
            needlr_row = needlr_by_gene.get(gene, [{}])[0] if gene != MISSING else {}
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
                ("SUPPORT", support),
                ("GENES", gene),
                ("NEEDLR_AF", first(needlr_row, ["AF", "MAX_AF", "AF_MAX", "SV_AF", "AF_1KGP", "1KGP_AF"])),
                ("NEEDLR_HPO", first(needlr_row, ["HPO", "HPO_terms", "HPO_Terms", "HPO phenotypes"])),
                ("OMIM", first(annotsv_row, ["OMIM", "AnnotSV_OMIM", "AnnotSV_OMIM_evidence"])),
                ("GENCC", first(annotsv_row, ["GENCC", "GenCC", "AnnotSV_GENCC", "AnnotSV_GENCC_evidence"])),
                ("ANNotsv_Gene", first(annotsv_row, ["Gene_name", "Gene", "GENE", "Genes", "gene", "SYMBOL", "AnnotSV_Gene"])),
                ("ANNotsv_Classification", first(annotsv_row, ["AnnotSV ranking", "AnnotSV_rank", "AnnotSV_Classification"])),
                ("PANEL_STATUS", "PANEL_GENE" if gene in panel else ("UNRESOLVED" if gene == MISSING else "NONPANEL_GENE")),
                ("PHENOTYPE_SCORE", first(ranking_row, ["phenotype_score", "PHENOTYPE_SCORE"])),
                ("CANDIDATE_CLASS", first(ranking_row, ["classification", "CANDIDATE_CLASS"])),
                ("NEEDLR_STATUS", "ANNOTATED" if needlr_row else "NO_MATCH"),
            ])

            # If AnnotSV didn't carry OMIM/GENCC for this gene, needLR
            # sometimes does — fall back to it before giving up.
            if row["OMIM"] == MISSING:
                row["OMIM"] = first(needlr_row, ["OMIM", "OMIM_phenotypes", "OMIM phenotypes"])
            if row["GENCC"] == MISSING:
                row["GENCC"] = first(needlr_row, ["GENCC", "GenCC", "GenCC_phenotypes", "GenCC phenotypes"])

            # Carry through every raw INFO_* column collected from the VCF.
            for key, value in sv.items():
                if key.startswith("INFO_"):
                    row[key] = value

            output_rows.append(row)

    fixed_columns = [
        "SV_ID", "CHROM", "START", "END", "SVTYPE", "SVLEN", "QUAL", "FILTER",
        "CALLERS", "SUPPORT", "GENES", "NEEDLR_AF", "NEEDLR_HPO", "OMIM", "GENCC",
        "ANNotsv_Gene", "ANNotsv_Classification", "PANEL_STATUS", "PHENOTYPE_SCORE",
        "CANDIDATE_CLASS", "NEEDLR_STATUS",
    ]
    info_columns = sorted({key for row in output_rows for key in row if key.startswith("INFO_")})
    columns = fixed_columns + info_columns

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t",
            extrasaction="ignore", lineterminator="\n",
        )
        writer.writeheader()
        for row in output_rows:
            writer.writerow({column: row.get(column, MISSING) for column in columns})

    print(f"[OK] integrated_rows={len(output_rows)} info_columns={len(info_columns)}")
    print(f"[OK] output={output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
