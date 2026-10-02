#!/usr/bin/env python3
"""Summarize Ensembl VEP structural-variant consequences per master SV and gene.

The output is intentionally descriptive. It preserves transcript-level VEP
consequences, exon/intron numbering and whole-transcript effects without
converting them into pathogenicity evidence.

Expected VEP input is tab-delimited output containing at least:
Uploaded_variation, Location, Gene, Feature, Feature_type, Consequence, SYMBOL.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

MISSING = {None, "", ".", "-", "NA", "N/A", "NAN", "NONE", "NULL"}

csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(
        path, "r", encoding="utf-8", errors="replace"
    )


def clean(value) -> str:
    if value is None:
        return "."
    value = str(value).strip()
    return "." if value.upper() in MISSING else value


def split_terms(value) -> set[str]:
    value = clean(value)
    if value == ".":
        return set()
    return {x.strip() for x in re.split(r"[,;&|]", value) if x.strip()}


def parse_info(raw: str) -> dict[str, str]:
    out = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "True"
    return out


def normalize_chrom(value: str) -> str:
    value = clean(value)
    return value[3:] if value.lower().startswith("chr") else value


def normalize_svtype(value: str) -> str:
    value = clean(value).upper().strip("<>")
    aliases = {
        "DELETION": "DEL",
        "DUPLICATION": "DUP",
        "INSERTION": "INS",
        "INVERSION": "INV",
        "BREAKEND": "BND",
        "TRANSLOCATION": "TRA",
    }
    return aliases.get(value, value)


def parse_location(value: str):
    value = clean(value)
    match = re.fullmatch(r"([^:]+):(\d+)(?:-(\d+))?", value)
    if not match:
        return None
    chrom = normalize_chrom(match.group(1))
    start = int(match.group(2))
    end = int(match.group(3) or match.group(2))
    return chrom, start, end


def read_master_vcf(path: str):
    by_id = {}
    by_coord = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, sv_id, _ref, alt, _qual, _filter, info_raw = fields[:8]
            info = parse_info(info_raw)
            svtype = normalize_svtype(info.get("SVTYPE", "."))
            if svtype == ".":
                if "[" in alt or "]" in alt:
                    svtype = "BND"
                elif alt.startswith("<") and alt.endswith(">"):
                    svtype = normalize_svtype(alt)

            start = int(pos)
            try:
                end = int(info.get("END", start))
            except ValueError:
                end = start

            if sv_id in MISSING:
                raise ValueError("Master VCF requires non-missing IDs for VEP integration")

            record = {
                "SV_ID": sv_id,
                "CHROM": normalize_chrom(chrom),
                "START": start,
                "END": end,
                "SVTYPE": svtype,
            }
            by_id[sv_id] = record
            by_coord[(record["CHROM"], start, end)].append(record)
            if end != start:
                by_coord[(record["CHROM"], start, start)].append(record)
    return by_id, by_coord


def read_vep(path: str) -> list[dict]:
    header = None
    rows = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                fields = line.lstrip("#").rstrip("\n").split("\t")
                if "Uploaded_variation" in fields:
                    header = fields
                continue
            if header is None:
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                values.extend(["."] * (len(header) - len(values)))
            row = dict(zip(header, values))

            # Legacy default VEP output stores SYMBOL/BIOTYPE/CANONICAL/PICK
            # and other annotations inside the semicolon-delimited Extra field.
            extra = clean(row.get("Extra"))
            if extra != ".":
                for item in extra.split(";"):
                    if "=" not in item:
                        continue
                    key, value = item.split("=", 1)
                    row.setdefault(key, value)

            rows.append(row)
    if header is None:
        raise ValueError(f"No VEP tabular header found in {path}")
    return rows


def infer_svtype_from_vep(row: dict) -> str:
    allele = normalize_svtype(row.get("Allele", "."))
    if allele in {"DEL", "DUP", "INS", "INV", "BND", "TRA", "CNV", "TDUP"}:
        return allele
    variant_class = normalize_svtype(row.get("VARIANT_CLASS", "."))
    return variant_class


def resolve_sv(row: dict, by_id: dict, by_coord: dict):
    uploaded = clean(row.get("Uploaded_variation"))
    if uploaded in by_id:
        return by_id[uploaded], "SV_ID"

    location = parse_location(row.get("Location", "."))
    if location:
        chrom, start, end = location
        candidates = by_coord.get((chrom, start, end), [])
        if len(candidates) == 1:
            return candidates[0], "EXACT_LOCATION"
        if len(candidates) > 1:
            vep_type = infer_svtype_from_vep(row)
            typed = [r for r in candidates if r["SVTYPE"] == vep_type]
            if len(typed) == 1:
                return typed[0], "EXACT_LOCATION_AND_TYPE"
    return None, "NO_MATCH"


WHOLE_LOSS = {"transcript_ablation"}
WHOLE_GAIN = {"transcript_amplification"}
TRUNCATION = {"feature_truncation"}
ELONGATION = {"feature_elongation"}
SPLICE = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
    "splice_donor_region_variant",
    "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant",
}
EXONIC = {
    "exon_loss_variant",
    "coding_sequence_variant",
    "protein_altering_variant",
    "frameshift_variant",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "missense_variant",
    "synonymous_variant",
    "inframe_insertion",
    "inframe_deletion",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
}
INTRONIC = {"intron_variant"}
PROXIMAL = {"upstream_gene_variant", "downstream_gene_variant"}


def classify_row(row: dict) -> tuple[set[str], set[str]]:
    consequences = split_terms(row.get("Consequence"))
    effects = set()
    regions = set()

    if consequences & WHOLE_LOSS:
        effects.add("WHOLE_TRANSCRIPT_LOSS")
        regions.add("WHOLE_TRANSCRIPT")
    if consequences & WHOLE_GAIN:
        effects.add("WHOLE_TRANSCRIPT_GAIN")
        regions.add("WHOLE_TRANSCRIPT")
    if consequences & TRUNCATION:
        effects.add("TRANSCRIPT_TRUNCATION")
    if consequences & ELONGATION:
        effects.add("TRANSCRIPT_ELONGATION")
    if consequences & SPLICE:
        effects.add("SPLICE_RELEVANT")
        regions.add("EXONIC_OR_SPLICE")
    if consequences & EXONIC or clean(row.get("EXON")) != ".":
        effects.add("EXONIC_TRANSCRIPT_OVERLAP")
        regions.add("EXONIC_OR_SPLICE")
    if consequences & INTRONIC or clean(row.get("INTRON")) != ".":
        effects.add("INTRONIC_TRANSCRIPT_OVERLAP")
        regions.add("INTRONIC")
    if consequences & PROXIMAL:
        regions.add("PROXIMAL")
    if not regions:
        regions.add("OTHER_TRANSCRIPT_CONTEXT")
    if not effects:
        effects.add("NO_SPECIFIC_STRUCTURAL_EFFECT_RESOLVED")

    return effects, regions


def join_values(values) -> str:
    cleaned = sorted({clean(v) for v in values if clean(v) != "."})
    return ";".join(cleaned) if cleaned else "."


def aggregate_group(sv_id: str, gene: str, rows: list[dict]) -> dict:
    effects = set()
    regions = set()
    transcripts = set()
    whole_effect_transcripts = set()
    canonical = set()
    picked = set()
    overlap_bp = []
    overlap_pc = []

    for row in rows:
        row_effects, row_regions = classify_row(row)
        effects.update(row_effects)
        regions.update(row_regions)

        feature = clean(row.get("Feature"))
        if feature != ".":
            transcripts.add(feature)
            if row_effects & {"WHOLE_TRANSCRIPT_LOSS", "WHOLE_TRANSCRIPT_GAIN"}:
                whole_effect_transcripts.add(feature)
            if clean(row.get("CANONICAL")).upper() in {"YES", "1"}:
                canonical.add(feature)
            if clean(row.get("PICK")).upper() in {"YES", "1"}:
                picked.add(feature)

        for field, store in (("OverlapBP", overlap_bp), ("OverlapPC", overlap_pc)):
            value = clean(row.get(field))
            if value != ".":
                try:
                    store.append(float(value))
                except ValueError:
                    pass

    transcript_count = len(transcripts)
    whole_count = len(whole_effect_transcripts)
    if transcript_count and whole_count == transcript_count:
        gene_scope = "ALL_MATCHED_TRANSCRIPTS_WHOLE_EFFECT"
    elif whole_count:
        gene_scope = "SOME_MATCHED_TRANSCRIPTS_WHOLE_EFFECT"
    else:
        gene_scope = "NO_WHOLE_TRANSCRIPT_EFFECT"

    return {
        "SV_ID": sv_id,
        "GENE": gene,
        "VEP_MATCH_METHOD": join_values(r.get("_MATCH_METHOD") for r in rows),
        "VEP_GENE_IDS": join_values(r.get("Gene") for r in rows),
        "VEP_TRANSCRIPT_COUNT": str(transcript_count),
        "VEP_TRANSCRIPTS": join_values(transcripts),
        "VEP_WHOLE_TRANSCRIPT_COUNT": str(whole_count),
        "VEP_GENE_TRANSCRIPT_SCOPE": gene_scope,
        "VEP_CONSEQUENCES": join_values(
            term for r in rows for term in split_terms(r.get("Consequence"))
        ),
        "VEP_IMPACTS": join_values(r.get("IMPACT") for r in rows),
        "VEP_BIOTYPES": join_values(r.get("BIOTYPE") for r in rows),
        "VEP_EXON": join_values(r.get("EXON") for r in rows),
        "VEP_INTRON": join_values(r.get("INTRON") for r in rows),
        "VEP_CANONICAL_TRANSCRIPTS": join_values(canonical),
        "VEP_PICK_TRANSCRIPTS": join_values(picked),
        "VEP_OVERLAP_BP_MAX": str(int(max(overlap_bp))) if overlap_bp else ".",
        "VEP_OVERLAP_PC_MAX": f"{max(overlap_pc):.3f}" if overlap_pc else ".",
        "VEP_TRANSCRIPT_REGION_CLASS": ";".join(sorted(regions)),
        "VEP_STRUCTURAL_EFFECT": ";".join(sorted(effects)),
        "VEP_ROWS_JSON": json.dumps(
            [
                {
                    key: clean(r.get(key))
                    for key in (
                        "Uploaded_variation",
                        "Location",
                        "Allele",
                        "Gene",
                        "Feature",
                        "Feature_type",
                        "Consequence",
                        "IMPACT",
                        "SYMBOL",
                        "BIOTYPE",
                        "EXON",
                        "INTRON",
                        "CANONICAL",
                        "PICK",
                        "VARIANT_CLASS",
                        "OverlapBP",
                        "OverlapPC",
                        "DISTANCE",
                        "STRAND",
                    )
                }
                for r in rows
            ],
            separators=(",", ":"),
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vep", required=True)
    parser.add_argument("--vcf", required=True, help="Same Jasmine master VCF used as VEP input")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    by_id, by_coord = read_master_vcf(args.vcf)
    vep_rows = read_vep(args.vep)

    grouped = defaultdict(list)
    unmatched = 0
    non_transcript = 0
    no_symbol = 0

    for row in vep_rows:
        if clean(row.get("Feature_type")) != "Transcript":
            non_transcript += 1
            continue

        sv, method = resolve_sv(row, by_id, by_coord)
        if sv is None:
            unmatched += 1
            continue

        gene = clean(row.get("SYMBOL")).upper()
        if gene == ".":
            no_symbol += 1
            continue

        row["_MATCH_METHOD"] = method
        grouped[(sv["SV_ID"], gene)].append(row)

    columns = [
        "SV_ID",
        "GENE",
        "VEP_MATCH_METHOD",
        "VEP_GENE_IDS",
        "VEP_TRANSCRIPT_COUNT",
        "VEP_TRANSCRIPTS",
        "VEP_WHOLE_TRANSCRIPT_COUNT",
        "VEP_GENE_TRANSCRIPT_SCOPE",
        "VEP_CONSEQUENCES",
        "VEP_IMPACTS",
        "VEP_BIOTYPES",
        "VEP_EXON",
        "VEP_INTRON",
        "VEP_CANONICAL_TRANSCRIPTS",
        "VEP_PICK_TRANSCRIPTS",
        "VEP_OVERLAP_BP_MAX",
        "VEP_OVERLAP_PC_MAX",
        "VEP_TRANSCRIPT_REGION_CLASS",
        "VEP_STRUCTURAL_EFFECT",
        "VEP_ROWS_JSON",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(grouped):
            writer.writerow(aggregate_group(key[0], key[1], grouped[key]))

    print(
        f"[OK] VEP SV-gene summaries={len(grouped)} "
        f"unmatched_transcript_rows={unmatched} "
        f"rows_without_symbol={no_symbol} "
        f"non_transcript_rows_skipped={non_transcript} "
        f"output={output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
