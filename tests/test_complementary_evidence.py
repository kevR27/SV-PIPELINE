#!/usr/bin/env python3
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "plots")]

import build_cohort_sv_recurrence as recurrence


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ComplementaryEvidenceTests(unittest.TestCase):
    def test_straglr_inside_large_del_is_context_not_support(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            write_tsv(
                d / "integrated.tsv",
                [
                    {
                        "SV_ID": "large_del",
                        "GENES": "A",
                        "CHROM": "chr1",
                        "START": "100",
                        "END": "10000",
                        "SVTYPE": "DEL",
                    },
                    {
                        "SV_ID": "local_ins",
                        "GENES": "B",
                        "CHROM": "chr1",
                        "START": "500",
                        "END": "500",
                        "SVTYPE": "INS",
                    },
                ],
            )
            write_tsv(
                d / "straglr.tsv",
                [
                    {
                        "chrom": "chr1",
                        "start": "490",
                        "end": "520",
                        "locus": "repeat1",
                        "overlapping_genes": "B",
                        "copy_number": "25",
                        "supporting_reads": "8",
                        "genotype": "10/25",
                    }
                ],
            )

            proc = subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "plots" / "intersect_complementary_sv_evidence.py"),
                    "--integrated",
                    str(d / "integrated.tsv"),
                    "--straglr",
                    str(d / "straglr.tsv"),
                    "--output",
                    str(d / "out.tsv"),
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)

            rows = {row["SV_ID"]: row for row in read_tsv(d / "out.tsv")}
            self.assertEqual(rows["large_del"]["STRAGLR_MATCH"], "NO")
            self.assertEqual(rows["large_del"]["STRAGLR_CONTEXT"], "YES")
            self.assertIn("repeat1", rows["large_del"]["STRAGLR_CONTEXT_LOCI"])

            self.assertEqual(rows["local_ins"]["STRAGLR_MATCH"], "YES")
            self.assertEqual(rows["local_ins"]["STRAGLR_CONTEXT"], "NO")

    def test_swapped_bnd_orientation_is_reversed(self):
        a = recurrence.SV(
            sample="A",
            sv_id="a",
            svtype="BND",
            chrom="chr1",
            start=1000,
            end=1000,
            chrom2="chr2",
            pos2=2000,
            size=None,
            orientation="+-",
        )
        swapped_equivalent = recurrence.SV(
            sample="B",
            sv_id="b",
            svtype="BND",
            chrom="chr2",
            start=2000,
            end=2000,
            chrom2="chr1",
            pos2=1000,
            size=None,
            orientation="-+",
        )
        swapped_wrong_orientation = recurrence.SV(
            sample="C",
            sv_id="c",
            svtype="BND",
            chrom="chr2",
            start=2000,
            end=2000,
            chrom2="chr1",
            pos2=1000,
            size=None,
            orientation="+-",
        )

        self.assertTrue(recurrence.same_sv(a, swapped_equivalent, 1000, 500))
        self.assertFalse(recurrence.same_sv(a, swapped_wrong_orientation, 1000, 500))


if __name__ == "__main__":
    unittest.main()
