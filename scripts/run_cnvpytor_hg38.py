#!/usr/bin/env python3
"""Run CNVpytor for GRCh38 without requiring other genomes' resource files.

CNVpytor 1.3.2 checks every built-in genome before parsing its arguments.
For this GRCh38-only WGS workflow, retain just hg38 in the registered
reference dictionary. Its GC/mask files must still exist.
"""

from pathlib import Path
import sys

from cnvpytor.genome import Genome
from cnvpytor.__main__ import main as cnvpytor_main


def main():
    hg38 = Genome.reference_genomes.get("hg38")
    if hg38 is None:
        raise RuntimeError("CNVpytor does not define the hg38 reference genome.")

    for key in ("gc_file", "mask_file"):
        path = hg38.get(key)
        if path is None or not Path(path).is_file() or Path(path).stat().st_size == 0:
            raise FileNotFoundError(
                f"CNVpytor hg38 {key} resource is missing or empty: {path}"
            )

    # Restrict only this process, not the installed CNVpytor environment.
    # Other genomes are irrelevant to the HG002 GRCh38 analysis.
    Genome.reference_genomes = {"hg38": hg38}

    return cnvpytor_main()


if __name__ == "__main__":
    sys.exit(main() or 0)
