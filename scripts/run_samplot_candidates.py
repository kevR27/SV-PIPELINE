#!/usr/bin/env python3
"""Create Samplot images for the highest-priority SRS SV candidates."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

import pandas as pd


def find_column(table: pd.DataFrame, names: list[str]) -> str | None:
    """Find the first matching column, ignoring capitalization."""
    lower_names = {str(column).lower(): column for column in table.columns}
    for name in names:
        if name in table.columns:
            return name
        if name.lower() in lower_names:
            return lower_names[name.lower()]
    return None


def integer(value) -> int | None:
    """Convert a table value to an integer coordinate when possible."""
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def safe_filename_text(value, maximum_length: int = 55) -> str:
    """Remove characters that are awkward in image filenames."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")
    return (cleaned or "candidate")[:maximum_length]


def samplot_command(
    sample: str,
    bam: str,
    reference: str,
    image: Path,
    chromosome: str,
    start: int,
    end: int,
    svtype: str,
    minimum_mapq: int,
    window: int,
) -> list[str]:
    """Build the short-read Samplot command for one interval."""
    command = [
        "samplot", "plot",
        "-n", sample,
        "-b", bam,
        "-r", reference,
        "-o", str(image),
        "-c", chromosome,
        "-s", str(start),
        "-e", str(end),
        "-q", str(minimum_mapq),
        "--window", str(window),
    ]
    if svtype in {"DEL", "DUP", "INV"}:
        command.extend(["-t", svtype])
    return command


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--window", type=int, default=5000)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t", dtype=str, low_memory=False)
    id_column = find_column(table, ["SV_ID", "ID"])
    chromosome_column = find_column(table, ["CHROM", "Chr"])
    start_column = find_column(table, ["START", "POS"])
    end_column = find_column(table, ["END"])
    type_column = find_column(table, ["SVTYPE", "SV_type"])
    gene_column = find_column(table, ["GENE", "GENES", "gene"])

    if not all([id_column, chromosome_column, start_column, type_column]):
        raise ValueError("Input requires SV_ID, CHROM, START and SVTYPE.")

    selected = table.drop_duplicates(id_column).head(args.top_n)
    image_directory = Path(args.out_dir)
    image_directory.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for rank, (_, row) in enumerate(selected.iterrows(), start=1):
        chromosome = str(row[chromosome_column])
        start = integer(row[start_column])
        end = integer(row[end_column]) if end_column else None
        svtype = str(row[type_column]).upper()
        gene = str(row[gene_column]) if gene_column else "."
        sv_id = str(row[id_column])

        status = "SKIPPED"
        message = "."
        image_path = "."

        # Interchromosomal events require both breakpoints. Drawing them as a
        # same-chromosome interval would be misleading, so keep them for manual
        # two-breakpoint review instead.
        if svtype in {"BND", "TRA", "CTX"}:
            message = "Interchromosomal event retained for manual review"
        elif start is not None:
            end = end if end is not None and end > start else start + 1
            image = image_directory / (
                f"{rank:03d}_{safe_filename_text(svtype, 12)}_"
                f"{safe_filename_text(chromosome, 12)}_{start}_{end}_"
                f"{safe_filename_text(sv_id)}_{safe_filename_text(gene, 30)}.png"
            )

            command = samplot_command(
                args.sample,
                args.bam,
                args.reference,
                image,
                chromosome,
                start,
                end,
                svtype,
                args.min_mapq,
                args.window,
            )
            result = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            image_created = image.exists() and image.stat().st_size > 0
            if result.returncode == 0 and image_created:
                status = "OK"
                image_path = str(image)
            else:
                status = "FAILED"
                message = (result.stdout or "").strip().split("\n")[-1]

        manifest_rows.append(
            {
                "RANK": rank,
                "SV_ID": sv_id,
                "GENE": gene,
                "SVTYPE": svtype,
                "CHROM": chromosome,
                "START": start if start is not None else ".",
                "END": end if end is not None else ".",
                "STATUS": status,
                "IMAGE": image_path,
                "MESSAGE": message,
            }
        )

    manifest = Path(args.manifest)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "RANK", "SV_ID", "GENE", "SVTYPE", "CHROM", "START", "END",
        "STATUS", "IMAGE", "MESSAGE",
    ]
    pd.DataFrame(manifest_rows, columns=columns).to_csv(
        manifest,
        sep="\t",
        index=False,
    )

    rendered = sum(row["STATUS"] == "OK" for row in manifest_rows)
    print(
        f"[OK] requested={len(selected)} rendered={rendered} "
        f"manifest={manifest}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
