#!/usr/bin/env python3
"""Check the prepared Ensembl gene BED before AnnotSV uses it."""
import argparse
from pathlib import Path
from sv_evidence_common import invalid_gene_labels


def check_resource(path):
    count, bad, examples = 0, 0, []
    with open(path) as f:
        for line_number, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            count += 1
            if len(fields) != 10 or invalid_gene_labels(fields[4]):
                bad += 1
                if len(examples) < 5:
                    examples.append(str(line_number))
    if not count or bad:
        raise ValueError(f"AnnotSV gene resource needs repair: {path}; records={count}, invalid={bad}, "
                         f"example lines={','.join(examples) or '.'}. Expected the 10-column Ensembl BED, "
                         "with a gene name in column 5. cmpl/incmpl are CDS status labels. "
                         "Rebuild this resource from the matching gene annotation release before rerunning AnnotSV.")
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotations-dir", required=True)
    p.add_argument("--genome-build", default="GRCh38")
    a = p.parse_args()
    root = Path(a.annotations_dir) / "Annotations_Human" / "Genes" / a.genome_build
    paths = [root / "genes.ENSEMBL.sorted.bed", root / "genes.ENSEMBL.sorted.tmp.bed"]
    path = next((p for p in paths if p.is_file()), None)
    if path is None:
        p.error(f"Prepared Ensembl gene resource not found in {root}")
    print(f"[OK] Checked {check_resource(path)} Ensembl gene records in {path}")


if __name__ == "__main__":
    main()
