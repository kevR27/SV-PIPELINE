#!/usr/bin/env python3
"""
Summarize all VEP transcript consequences for each SV/gene pair.

VEP is run with --flag_pick, not --pick.
Therefore all transcript consequences are retained. PICK and CANONICAL are
reported as useful annotations but are not used to discard other transcripts.

This prevents a large structural variant affecting several transcripts from
being represented by only one selected transcript.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

MISSING = "."


def parse_extra(text: str) -> dict[str, str]:
    """Parse VEP's semicolon-separated Extra column."""
    result: dict[str, str] = {}

    for item in str(text or "").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = "1"

    return result


def split_terms(value: str) -> list[str]:
    """Split comma/ampersand separated VEP terms."""
    if value in (None, "", MISSING):
        return []

    return [
        term
        for term in re.split(r"[,;&]", str(value))
        if term and term != MISSING
    ]


def numeric_values(items, keys: tuple[str, ...]) -> list[float]:
    """Collect simple numeric VEP Extra fields when present."""
    values: list[float] = []

    for _, extra in items:
        for key in keys:
            value = extra.get(key)

            if value in (None, "", MISSING):
                continue

            if re.fullmatch(r"[0-9.]+", value):
                values.append(float(value))

    return values


def transcript_region_class(consequences: set[str]) -> str:
    """Create a readable summary of the transcript regions affected."""
    regions: set[str] = set()

    for consequence in consequences:
        if "transcript_ablation" in consequence:
            regions.add("TRANSCRIPT_ABLATION")
        elif "exon_loss" in consequence:
            regions.add("EXON_LOSS")
        elif "coding_sequence_variant" in consequence:
            regions.add("CODING")
        elif "intron_variant" in consequence:
            regions.add("INTRONIC")
        elif (
            "upstream" in consequence
            or "downstream" in consequence
        ):
            regions.add("REGULATORY_PROXIMAL")

    return ";".join(sorted(regions)) or MISSING


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize all VEP SV transcript rows by SV and gene."
    )
    parser.add_argument("--vep", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    header = None
    groups = defaultdict(list)

    with open(
        args.vep,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#Uploaded_variation"):
                header = line[1:].rstrip("\n").split("\t")
                continue

            if line.startswith("#") or not line.strip():
                continue

            if header is None:
                continue

            fields = line.rstrip("\n").split("\t")
            row = dict(zip(header, fields))
            extra = parse_extra(row.get("Extra", ""))

            sv_id = row.get(
                "Uploaded_variation",
                MISSING,
            )
            gene = extra.get(
                "SYMBOL",
                MISSING,
            )

            if gene == MISSING:
                continue

            groups[(sv_id, gene.upper())].append(
                (row, extra)
            )

    output_rows: list[dict[str, str]] = []

    for (sv_id, gene), items in sorted(groups.items()):
        gene_ids = {
            row.get("Gene", MISSING)
            for row, _ in items
            if row.get("Gene", MISSING) != MISSING
        }

        transcripts = {
            row.get("Feature", MISSING)
            for row, _ in items
            if row.get("Feature", MISSING) != MISSING
        }

        consequences = {
            term
            for row, _ in items
            for term in split_terms(
                row.get("Consequence", MISSING)
            )
        }

        impacts = {
            extra.get("IMPACT", MISSING)
            for _, extra in items
            if extra.get("IMPACT", MISSING) != MISSING
        }

        biotypes = {
            extra.get("BIOTYPE", MISSING)
            for _, extra in items
            if extra.get("BIOTYPE", MISSING) != MISSING
        }

        exons = {
            row.get(
                "EXON",
                extra.get("EXON", MISSING),
            )
            for row, extra in items
            if row.get(
                "EXON",
                extra.get("EXON", MISSING),
            )
            != MISSING
        }

        introns = {
            row.get(
                "INTRON",
                extra.get("INTRON", MISSING),
            )
            for row, extra in items
            if row.get(
                "INTRON",
                extra.get("INTRON", MISSING),
            )
            != MISSING
        }

        canonical_transcripts = {
            row.get("Feature", MISSING)
            for row, extra in items
            if extra.get("CANONICAL") == "YES"
        }

        pick_transcripts = {
            row.get("Feature", MISSING)
            for row, extra in items
            if extra.get("PICK") in {"1", "YES"}
        }

        overlap_bp = numeric_values(
            items,
            ("BP_OVERLAP", "OVERLAP_BP"),
        )
        overlap_percent = numeric_values(
            items,
            ("PERCENT_OVERLAP", "OVERLAP_PC"),
        )

        output_rows.append(
            {
                "SV_ID": sv_id,
                "GENE": gene,
                "VEP_MATCH_METHOD": "SV_ID_AND_SYMBOL",
                "VEP_GENE_IDS": (
                    ";".join(sorted(gene_ids))
                    or MISSING
                ),
                "VEP_TRANSCRIPT_COUNT": str(
                    len(transcripts)
                ),
                "VEP_TRANSCRIPTS": (
                    ";".join(sorted(transcripts))
                    or MISSING
                ),
                "VEP_WHOLE_TRANSCRIPT_COUNT": str(
                    len(transcripts)
                ),
                "VEP_GENE_TRANSCRIPT_SCOPE": (
                    "ALL_TRANSCRIPTS_RETAINED_FLAG_PICK_ONLY"
                ),
                "VEP_CONSEQUENCES": (
                    ";".join(sorted(consequences))
                    or MISSING
                ),
                "VEP_IMPACTS": (
                    ";".join(sorted(impacts))
                    or MISSING
                ),
                "VEP_BIOTYPES": (
                    ";".join(sorted(biotypes))
                    or MISSING
                ),
                "VEP_EXON": (
                    ";".join(sorted(exons))
                    or MISSING
                ),
                "VEP_INTRON": (
                    ";".join(sorted(introns))
                    or MISSING
                ),
                "VEP_CANONICAL_TRANSCRIPTS": (
                    ";".join(sorted(canonical_transcripts))
                    or MISSING
                ),
                "VEP_PICK_TRANSCRIPTS": (
                    ";".join(sorted(pick_transcripts))
                    or MISSING
                ),
                "VEP_OVERLAP_BP_MAX": (
                    str(max(overlap_bp))
                    if overlap_bp
                    else MISSING
                ),
                "VEP_OVERLAP_PC_MAX": (
                    str(max(overlap_percent))
                    if overlap_percent
                    else MISSING
                ),
                "VEP_TRANSCRIPT_REGION_CLASS": (
                    transcript_region_class(consequences)
                ),
                "VEP_STRUCTURAL_EFFECT": (
                    ";".join(sorted(consequences))
                    or MISSING
                ),
            }
        )

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
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(
        f"[OK] SV_gene_groups={len(output_rows)} output={output}"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
