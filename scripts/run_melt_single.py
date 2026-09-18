#!/usr/bin/env python3
"""Run one MELTv2 Single analysis and expose its final comparison VCF.

MELT names its final output from BAM/transposon metadata rather than from an
explicit output argument.  This wrapper gives Snakemake a deterministic output
path and isolates each transposon-family run so stale VCFs cannot be selected.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one MELTv2 Single family")
    parser.add_argument("--jar", required=True, help="Path to MELT.jar")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--genes-bed", required=True)
    parser.add_argument("--transposon", required=True, help="MELT transposon ZIP")
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--memory-gb", type=int, default=8)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    workdir = Path(args.workdir).resolve()
    output = Path(args.output).resolve()
    run_dir = workdir / "run"

    # This directory belongs exclusively to one Snakemake output. Recreating
    # it avoids accidentally reusing an older *.final_comp.vcf.
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    command = [
        "java",
        f"-Xmx{args.memory_gb}G",
        "-jar",
        str(Path(args.jar).resolve()),
        "Single",
        "-bamfile",
        str(Path(args.bam).resolve()),
        "-h",
        str(Path(args.reference).resolve()),
        "-n",
        str(Path(args.genes_bed).resolve()),
        "-t",
        str(Path(args.transposon).resolve()),
        "-w",
        str(run_dir),
    ]
    subprocess.run(command, check=True)

    candidates = sorted(run_dir.rglob("*.final_comp.vcf"))
    if len(candidates) != 1:
        rendered = ", ".join(str(path) for path in candidates) or "none"
        raise RuntimeError(
            "Expected exactly one MELT *.final_comp.vcf in "
            f"{run_dir}, found {len(candidates)}: {rendered}"
        )

    shutil.copyfile(candidates[0], output)
    if output.stat().st_size == 0:
        raise RuntimeError(f"MELT output is empty: {output}")

    print(f"[OK] MELT final VCF: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
