import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from annotate_genomiser_comparison import annotate_genomiser
from compare_genomiser_rankings import compare
from assess_sv_alleles import hypothesis_key


class UnifiedRankingTests(unittest.TestCase):
    def test_comparison_preserves_order_and_separates_gene_and_variant(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "genes.tsv").write_text("#RANK\tGENE_SYMBOL\tEXOMISER_GENE_COMBINED_SCORE\n1\tA\t0.9\n2\tA\t0.8\n3\tB\t0.7\n")
            (d / "variants.tsv").write_text("GENE_SYMBOL\tCONTIG\tSTART\tEND\tALT\tEXOMISER_VARIANT_SCORE\nA\t1\t100\t200\t<DEL>\t0.8\n")
            source = pd.DataFrame([
                {"GENE": "B", "SV_ID": "b", "CHROM": "chr1", "START": 50, "END": 60, "SVTYPE": "INV", "CLINICAL_REVIEW_RANK": 1},
                {"GENE": "A", "SV_ID": "a", "CHROM": "chr1", "START": 100, "END": 200, "SVTYPE": "DEL", "CLINICAL_REVIEW_RANK": 2},
                {"GENE": "A", "SV_ID": "c", "CHROM": "chr1", "START": 300, "END": 400, "SVTYPE": "DEL", "CLINICAL_REVIEW_RANK": 3},
            ])
            out = annotate_genomiser(source, d / "genes.tsv", d / "variants.tsv")
            pd.testing.assert_frame_equal(out[source.columns], source)
            self.assertEqual(out.GENOMISER_GENE_RANK.tolist(), [2, 1, 1])
            self.assertEqual(out.GENOMISER_SV_VARIANT_SCORE.tolist(), [".", "0.8", "."])
            self.assertEqual(out.GENOMISER_SV_MATCH_STATUS.tolist()[1], "EXACT_COORDINATE_CANDIDATE")

    def test_descriptive_overlap_not_accuracy(self):
        result = compare({"A": 1, "B": 2}, {"C": 1, "A": 2}, ks=(1, 2))
        self.assertEqual(result["top_k"]["1"]["shared_n"], 0)
        self.assertEqual(result["top_k"]["2"]["shared_n"], 1)
        self.assertAlmostEqual(result["top_k"]["2"]["jaccard"], 1/3)
        self.assertIn("no truth-based", result["interpretation"].lower())

    def test_hypothesis_selection_is_categorical(self):
        compatible = {"ALLELE_DISEASE_HON_RELEVANCE": "SUPPORTS_HON",
                      "ALLELE_DISEASE_MECHANISM_STATUS": "DISEASE_MECHANISM_COMPATIBLE",
                      "ALLELE_DISEASE_VALIDITY": "STRONG"}
        conflict = {**compatible, "ALLELE_DISEASE_MECHANISM_STATUS": "MECHANISM_MISMATCH_REVIEW"}
        self.assertLess(hypothesis_key(compatible), hypothesis_key(conflict))
        self.assertEqual(hypothesis_key(compatible), hypothesis_key({**compatible, "ALLELE_EVIDENCE_COMPLETENESS": 0}))


if __name__ == "__main__":
    unittest.main()
