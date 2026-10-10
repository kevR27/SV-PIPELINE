#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from annotate_gnomad_sv_exact import exact_match, fetch_records


class FakeRecord:
    def __init__(self, pos, stop, svtype, svlen=None):
        self.pos = pos
        self.stop = stop
        self.info = {"SVTYPE": svtype}
        if svlen is not None:
            self.info["SVLEN"] = svlen
        self.alts = (f"<{svtype}>",)


class TestGnomadExactMatch(unittest.TestCase):
    def test_index_errors_are_not_converted_to_no_population_match(self):
        class BrokenVcf:
            header = SimpleNamespace(contigs={"chr1": 10000})

            def fetch(self, *args):
                raise ValueError("fetch requires an index")

        with self.assertRaisesRegex(ValueError, "index"):
            fetch_records(BrokenVcf(), "chr1", 100)
        self.assertEqual(fetch_records(BrokenVcf(), "chr2", 100), [])

    def test_exact_deletion_requires_same_pos_end_type(self):
        q = {"SVTYPE": "DEL", "START": "100", "END": "500", "SVLEN": "-400"}
        ok, scope = exact_match(q, FakeRecord(100, 500, "DEL", -400))
        self.assertTrue(ok)
        self.assertEqual(scope, "EXACT_POS_END_TYPE")

        ok, _ = exact_match(q, FakeRecord(101, 500, "DEL", -399))
        self.assertFalse(ok)

        ok, _ = exact_match(q, FakeRecord(100, 501, "DEL", -401))
        self.assertFalse(ok)

    def test_insertion_requires_same_position_and_size_when_available(self):
        q = {"SVTYPE": "INS", "START": "100", "END": "100", "SVLEN": "75"}
        ok, scope = exact_match(q, FakeRecord(100, 100, "INS", 75))
        self.assertTrue(ok)
        self.assertIn("NOT_SEQUENCE_CONFIRMED", scope)

        ok, _ = exact_match(q, FakeRecord(100, 100, "INS", 74))
        self.assertFalse(ok)

    def test_bnd_is_not_claimed_as_exact(self):
        q = {"SVTYPE": "BND", "START": "100", "END": "100", "SVLEN": "."}
        ok, scope = exact_match(q, FakeRecord(100, 100, "BND"))
        self.assertFalse(ok)
        self.assertEqual(scope, "SVTYPE_NOT_SUPPORTED_FOR_EXACT_MATCH")


if __name__ == "__main__":
    unittest.main()
