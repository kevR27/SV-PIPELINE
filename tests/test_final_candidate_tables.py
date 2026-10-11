#!/usr/bin/env python3
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_candidate_tables import build_gene_table
from build_diagnostic_candidates import build_clinical_table
from build_integrated_sv_gene_tsv import annotsv_row_for_gene

ROOT = Path(__file__).resolve().parents[1]


def clinical_evidence(**changes):
    row = {
        "SV_ID": "positive", "GENES": "OPA1", "SVTYPE": "DEL", "CHROM": "chr3",
        "ACMG_CNV_CLASS": "LIKELY_PATHOGENIC", "ACMG_CNV_SCORE": "0.90",
        "ANNOTSV_RANKING_CRITERIA": "2C-1 (0.90)",
        "ANNOTSV_CLASSIFICATION_SOURCE": "FULL_EVENT_ROW",
        "ANNOTSV_CLASSIFICATION_CONFLICT": "NO",
        "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED": "0.8",
        "GENE_DISEASE_EVIDENCE_LEVEL": "DEFINITIVE", "GENCC_MOI": "AD",
        "ALLELE_EFFECT": "POSSIBLE_LOF", "ALLELE_GENOTYPE": "HET",
        "ALLELE_GENOTYPE_SOURCE": "synthetic_exact_id_qc_genotype",
        "ALLELE_DISEASE_ID": "SYNTHETIC_HON_MODEL", "ALLELE_DISEASE_SOURCE": "synthetic_test_fixture",
        "ALLELE_DISEASE_HON_RELEVANCE": "SUPPORTS_HON", "ALLELE_DISEASE_HON_SOURCE": "synthetic_test_fixture",
        "ALLELE_DISEASE_MOI": "AD", "ALLELE_DISEASE_MECHANISM_STATUS": "DISEASE_MECHANISM_POSSIBLE",
        "ALLELE_DISEASE_MECHANISM_DETAIL": "disease=SYNTHETIC_HON_MODEL;mechanism=HAPLOINSUFFICIENCY;validity=DEFINITIVE;source=synthetic_test_fixture",
        "ALLELE_TECHNICAL_STATUS": "MULTICALLER_READ_SUPPORTED",
        "CLINGEN_HI": "3", "CALLER_COUNT": "2", "CALLER_READ_SUPPORT": "Sniffles2:8;cuteSV:7",
        "NEEDLR_AF": ".", "GENE_RELEVANCE_TIER": "HIGH",
        "GENE_RELEVANCE_DISPLAY_SCORE": "9", "PANEL_STATUS": "PANEL_GENE",
    }
    return {**row, **changes}


class TestClinicalCandidateTables(unittest.TestCase):
    def test_joint_rank_retains_uncertainty_unrelated_findings_and_benign(self):
        positive = clinical_evidence()
        vus = clinical_evidence(SV_ID="vus", ACMG_CNV_CLASS="VUS", ACMG_CNV_SCORE="0.3", ANNOTSV_RANKING_CRITERIA="4A (0.30)")
        unrelated = clinical_evidence(SV_ID="other_disease", HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED="0.1")
        benign = clinical_evidence(SV_ID="benign", ACMG_CNV_CLASS="BENIGN", ACMG_CNV_SCORE="-1")
        incomplete = clinical_evidence(SV_ID="unknown_gt", ALLELE_GENOTYPE="UNKNOWN")
        out = build_clinical_table(pd.DataFrame([benign, unrelated, vus, incomplete, positive]))
        self.assertEqual(out.SV_ID.tolist(), ["positive", "unknown_gt", "vus", "other_disease", "benign"])
        indexed = out.set_index("SV_ID")
        self.assertEqual(indexed.loc["vus", "CLINICAL_CNV_CLASS"], "VUS")
        self.assertEqual(indexed.loc["benign", "CLINICAL_REVIEW_GROUP"], "BENIGN_OR_LIKELY_BENIGN")
        self.assertEqual(indexed.loc["other_disease", "CLINICAL_HON_RELEVANCE"], "LIMITED_QUERY_MATCH")
        self.assertEqual(indexed.loc["positive", "CLINICAL_POPULATION_STATUS"], "UNKNOWN_AF")
        self.assertIn("CASE_OR_PHENOTYPE", indexed.loc["vus", "CLINICAL_PHENOTYPE_IN_CNV_CRITERIA"])
        self.assertEqual(indexed.loc["positive", "CLINICAL_PHENOTYPE_IN_CNV_CRITERIA"], "NOT_SHOWN_IN_REPORTED_CODES")
        family_criteria = build_clinical_table(pd.DataFrame([clinical_evidence(ANNOTSV_RANKING_CRITERIA="2C-1;5B")])).iloc[0]
        self.assertEqual(family_criteria.CLINICAL_PHENOTYPE_IN_CNV_CRITERIA, "CASE_OR_PHENOTYPE_CRITERIA_REPORTED")

    def test_gene_level_associations_are_not_disease_specific_confirmation(self):
        row = clinical_evidence(ALLELE_DISEASE_ID=".", ALLELE_DISEASE_SOURCE=".",
                                ALLELE_DISEASE_MECHANISM_STATUS="GENE_DOSAGE_SUPPORT_ONLY")
        out = build_clinical_table(pd.DataFrame([row])).iloc[0]
        self.assertEqual(out.CLINICAL_REVIEW_GROUP, "P_LP_HON_INCOMPLETE_EVIDENCE")
        self.assertIn("DISEASE_SPECIFIC_HON_MODEL_NOT_AVAILABLE", out.CLINICAL_MISSING_OR_CONTEXT_NOTES)
        self.assertEqual(out.CLINICAL_DISEASE_EVIDENCE_SCOPE, "GENE_LEVEL_ASSOCIATION_NOT_DISEASE_SPECIFIC")

    def test_disease_model_requires_reviewed_hon_link_not_just_gene_match(self):
        supported = build_clinical_table(pd.DataFrame([clinical_evidence()])).iloc[0]
        self.assertEqual(supported.CLINICAL_REVIEW_GROUP, "P_LP_HON_MECHANISM_SUPPORTED")
        cases = [clinical_evidence(SV_ID="unknown_hon_link", ALLELE_DISEASE_HON_RELEVANCE="UNKNOWN"),
                 clinical_evidence(SV_ID="unverified_hon_link", ALLELE_DISEASE_HON_SOURCE="."),
                 clinical_evidence(SV_ID="other_disease", ALLELE_DISEASE_HON_RELEVANCE="OTHER_DISEASE",
                                   ALLELE_DISEASE_HYPOTHESIS_COUNT="2")]
        out = build_clinical_table(pd.DataFrame(cases)).set_index("SV_ID")
        self.assertEqual(out.loc["unknown_hon_link", "CLINICAL_REVIEW_GROUP"], "P_LP_HON_INCOMPLETE_EVIDENCE")
        self.assertEqual(out.loc["unverified_hon_link", "CLINICAL_DISEASE_HON_STATUS"], "UNKNOWN")
        self.assertEqual(out.loc["other_disease", "CLINICAL_REVIEW_GROUP"], "P_LP_OTHER_OR_UNRESOLVED_HON")
        self.assertEqual(out.loc["other_disease", "CLINICAL_HON_RELEVANCE"], "STRONG_QUERY_MATCH")
        self.assertIn("MULTIPLE_DISEASE_MODELS_REVIEW", out.loc["other_disease", "CLINICAL_MISSING_OR_CONTEXT_NOTES"])

    def test_reviewed_population_af_is_used_and_unverified_af_stays_unknown(self):
        common = clinical_evidence(SV_ID="reviewed_common", ALLELE_POPULATION_STATUS="ABOVE_RESEARCH_AF_THRESHOLD",
                                   ALLELE_POPULATION_AF="0.2", ALLELE_POPULATION_MATCH="EXACT_ALLELE",
                                   ALLELE_POPULATION_SOURCE="synthetic_population_release", NEEDLR_AF="0.001")
        unverified = clinical_evidence(SV_ID="unverified", ALLELE_POPULATION_STATUS="UNKNOWN_UNVERIFIED_EXTERNAL_AF",
                                       ALLELE_POPULATION_AF="0.2", ALLELE_POPULATION_SOURCE="synthetic_population_release")
        out = build_clinical_table(pd.DataFrame([common, unverified])).set_index("SV_ID")
        self.assertEqual(out.loc["reviewed_common", "CLINICAL_MAX_REPORTED_AF"], 0.2)
        self.assertEqual(out.loc["reviewed_common", "CLINICAL_REVIEW_GROUP"], "CONFLICT_OR_TECHNICAL_REVIEW")
        self.assertIn("POPULATION_AF_SOURCE_CONFLICT", out.loc["reviewed_common", "CLINICAL_REVIEW_FLAGS"])
        self.assertIn("reviewed_allele:synthetic_population_release", out.loc["reviewed_common", "CLINICAL_POPULATION_SOURCES"])
        self.assertEqual(out.loc["unverified", "CLINICAL_POPULATION_STATUS"], "UNKNOWN_AF")
        legacy = clinical_evidence(ALLELE_POPULATION_STATUS="ABOVE_RESEARCH_AF_THRESHOLD",
                                   ALLELE_POPULATION_DETAIL="AF=0.2;source=synthetic_population_release;threshold=0.01")
        out = build_clinical_table(pd.DataFrame([legacy])).iloc[0]
        self.assertEqual(out.CLINICAL_MAX_REPORTED_AF, 0.2)
        self.assertEqual(out.CLINICAL_REVIEW_GROUP, "CONFLICT_OR_TECHNICAL_REVIEW")
        out = build_clinical_table(pd.DataFrame([{**legacy, "ALLELE_POPULATION_DETAIL": "AF=.;source=."}])).iloc[0]
        self.assertIn("REVIEWED_POPULATION_AF_UNAVAILABLE", out.CLINICAL_REVIEW_FLAGS)

    def test_passing_depth_only_cnv_is_retained_without_rescuing_bad_reads(self):
        depth = dict(ALLELE_TECHNICAL_STATUS="UNKNOWN_READ_SUPPORT", CALLER_COUNT="0", CALLER_READ_SUPPORT=".",
                     CNVPYTOR_READ_DEPTH_MATCH="DEPTH_ONLY", CNVPYTOR_STATUS="PASS", CNVPYTOR_LEVEL="0.5",
                     CNVPYTOR_EVAL1="1e-10", CNVPYTOR_Q0="0.01", CNVPYTOR_PN="0", CNVPYTOR_FLAGS="PASS")
        cases = [clinical_evidence(SV_ID="depth_only", **depth),
                 clinical_evidence(SV_ID="bad_reads", **{**depth, "ALLELE_TECHNICAL_STATUS": "LOW_READ_SUPPORT"}),
                 clinical_evidence(SV_ID="missing_depth_qc", **{**depth, "CNVPYTOR_Q0": "."}),
                 clinical_evidence(SV_ID="wrong_copy_direction", **{**depth, "CNVPYTOR_LEVEL": "1.5"}),
                 clinical_evidence(SV_ID="dup_depth", SVTYPE="DUP", ALLELE_EFFECT="POSSIBLE_COPY_GAIN",
                                   **{**depth, "CNVPYTOR_LEVEL": "1.5"})]
        out = build_clinical_table(pd.DataFrame(cases)).set_index("SV_ID")
        self.assertEqual(out.loc["depth_only", "CLINICAL_TECHNICAL_STATUS"], "DEPTH_SUPPORTED_BREAKPOINT_REVIEW")
        self.assertEqual(out.loc["depth_only", "CLINICAL_REVIEW_GROUP"], "P_LP_HON_INCOMPLETE_EVIDENCE")
        self.assertEqual(out.loc["dup_depth", "CLINICAL_REVIEW_GROUP"], "P_LP_HON_INCOMPLETE_EVIDENCE")
        for key in ("bad_reads", "missing_depth_qc", "wrong_copy_direction"):
            self.assertEqual(out.loc[key, "CLINICAL_REVIEW_GROUP"], "CONFLICT_OR_TECHNICAL_REVIEW")

    def test_unknown_effect_does_not_keep_a_strong_clinical_mechanism_label(self):
        row = clinical_evidence(ALLELE_EFFECT="UNKNOWN", INHERITANCE_MECHANISM_CLASS="STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY")
        out = build_clinical_table(pd.DataFrame([row])).iloc[0]
        self.assertEqual(out.CLINICAL_MECHANISM_CLASS, "GENE_EFFECT_UNRESOLVED")
        self.assertEqual(out.INHERITANCE_MECHANISM_CLASS, "STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY")
        self.assertEqual(out.CLINICAL_MECHANISM_SUPPORTED, "NO_OR_UNRESOLVED")
        row = clinical_evidence(ALLELE_EFFECT="GOF", INHERITANCE_MECHANISM_CLASS="UNRESOLVED",
                                ALLELE_DISEASE_MECHANISM_STATUS="DISEASE_MECHANISM_COMPATIBLE")
        out = build_clinical_table(pd.DataFrame([row])).iloc[0]
        self.assertEqual(out.CLINICAL_MECHANISM_CLASS, "REPORTED_EFFECT_COMPATIBLE_WITH_REVIEWED_DISEASE_MODEL")
        self.assertEqual(out.CLINICAL_MECHANISM_SUPPORTED, "YES")
        row = clinical_evidence(INHERITANCE_MECHANISM_CLASS="STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY",
                                ALLELE_DISEASE_MECHANISM_STATUS="MECHANISM_MISMATCH_REVIEW")
        out = build_clinical_table(pd.DataFrame([row])).iloc[0]
        self.assertEqual(out.CLINICAL_MECHANISM_CLASS, "DISEASE_MECHANISM_CONFLICT_REVIEW")
        self.assertEqual(out.CLINICAL_MECHANISM_SUPPORTED, "NO_OR_UNRESOLVED")

    def test_supporting_hon_match_and_unknown_af_are_retained_without_promotion(self):
        cases = [clinical_evidence(SV_ID="supporting_vus", ACMG_CNV_CLASS="VUS", ACMG_CNV_SCORE="0.3",
                                  HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED="0.3"),
                 clinical_evidence(SV_ID="unclassified", ACMG_CNV_CLASS=".", ACMG_CNV_SCORE=".", ANNOTSV_RANKING_CRITERIA="."),
                 clinical_evidence(SV_ID="zero_af", NEEDLR_AF="0"),
                 clinical_evidence(SV_ID="rare_af", NEEDLR_AF="0.01")]
        out = build_clinical_table(pd.DataFrame(cases)).set_index("SV_ID")
        self.assertEqual(out.loc["supporting_vus", "CLINICAL_REVIEW_GROUP"], "HON_RELEVANT_VUS_OR_UNCLASSIFIED")
        self.assertEqual(out.loc["unclassified", "CLINICAL_CNV_CLASS"], "UNCLASSIFIED")
        self.assertEqual(out.loc["unclassified", "CLINICAL_POPULATION_STATUS"], "UNKNOWN_AF")
        self.assertEqual(out.loc["zero_af", "CLINICAL_MAX_REPORTED_AF"], 0)
        self.assertEqual(out.loc["rare_af", "CLINICAL_POPULATION_STATUS"], "REPORTED_AF_LE_1_PERCENT")

    def test_aliases_patient_zero_and_ar_partner_review(self):
        row = clinical_evidence(GENE=".", ALLELE_DISEASE_MOI="AR", GENCC_MOI="AR", ALLELE_TRANS_PARTNER_IDS="second_allele",
                                HPO_SEMANTIC_STATUS="NO_REFERENCE_HPO_FOR_GENE", HPO_QUERY_RESNIK_NORMALIZED="0")
        out = build_clinical_table(pd.DataFrame([row])).iloc[0]
        self.assertEqual(out.GENE, "OPA1")
        self.assertEqual(out.CLINICAL_HON_QUERY_SIMILARITY, 0)
        self.assertEqual(out.CLINICAL_HON_QUERY_SCOPE, "PATIENT_OBSERVED_HPO_QUERY")
        self.assertEqual(out.CLINICAL_INHERITANCE_STATUS, "TRANS_PARTNER_CANDIDATE_REVIEW")

    def test_no_class_bonus_from_score_criteria_or_mitocarta(self):
        left = clinical_evidence(SV_ID="a", ACMG_CNV_CLASS="PATHOGENIC", ACMG_CNV_SCORE="1.0", ANNOTSV_RANKING_CRITERIA="2A")
        right = clinical_evidence(SV_ID="z", ACMG_CNV_CLASS="PATHOGENIC", ACMG_CNV_SCORE="3.0", ANNOTSV_RANKING_CRITERIA="2A;3C;4A;5H", MITOCARTA_ENCODING="NUCLEAR_MITOCHONDRIAL_GENE")
        out = build_clinical_table(pd.DataFrame([right, left]))
        self.assertEqual(out.SV_ID.tolist(), ["a", "z"])

    def test_conflicts_and_poor_support_cannot_be_rescued_by_hon(self):
        cases = [
            clinical_evidence(SV_ID="class_score_conflict", ACMG_CNV_SCORE="-1"),
            clinical_evidence(SV_ID="caller_qc", ALLELE_TECHNICAL_STATUS="REVIEW_REQUIRED"),
            clinical_evidence(SV_ID="reference_gt", ALLELE_GENOTYPE="HOM_REF"),
            clinical_evidence(SV_ID="mechanism_conflict", ALLELE_DISEASE_MECHANISM_STATUS="MECHANISM_MISMATCH_REVIEW"),
            clinical_evidence(SV_ID="depth_conflict", DEPTH_SUPPORT_CLASS="DOES_NOT_SHOW_EXPECTED_COPY_CHANGE"),
            clinical_evidence(SV_ID="population_conflict", NEEDLR_AF="0.2"),
            clinical_evidence(SV_ID="split_only", ANNOTSV_CLASSIFICATION_SOURCE="GENE_ROW_FALLBACK_REVIEW"),
        ]
        out = build_clinical_table(pd.DataFrame(cases))
        self.assertTrue(out.CLINICAL_REVIEW_GROUP.eq("CONFLICT_OR_TECHNICAL_REVIEW").all())
        self.assertEqual(len(out), len(cases))

    def test_non_cnv_and_multi_gene_event_do_not_prove_causality(self):
        cases = [
            clinical_evidence(SV_ID="inversion", SVTYPE="INV", ACMG_CNV_CLASS="PATHOGENIC", ACMG_CNV_SCORE="1.0"),
            clinical_evidence(SV_ID="large", GENES="OPA1", ACMG_CNV_CLASS="PATHOGENIC", ACMG_CNV_SCORE="1.0"),
            clinical_evidence(SV_ID="large", GENES="OTHER", ACMG_CNV_CLASS="PATHOGENIC", ACMG_CNV_SCORE="1.0", ALLELE_EFFECT="UNKNOWN", CLINGEN_HI="."),
        ]
        out = build_clinical_table(pd.DataFrame(cases)).set_index(["SV_ID", "GENE"])
        self.assertEqual(out.loc[("inversion", "OPA1"), "CLINICAL_CNV_CLASS"], "NOT_APPLICABLE_NON_CNV")
        self.assertEqual(out.loc[("large", "OTHER"), "CLINICAL_REVIEW_GROUP"], "P_LP_HON_INCOMPLETE_EVIDENCE")
        self.assertEqual(out.loc[("large", "OTHER"), "CLINICAL_MECHANISM_SUPPORTED"], "NO_OR_UNRESOLVED")

    def test_full_record_is_authoritative_and_disagreement_is_visible(self):
        full = {"Annotation_mode": "full", "Gene_name": "OPA1;OTHER", "ACMG_class": "3", "AnnotSV ranking score": "0.3", "AnnotSV ranking criteria": "4E"}
        split = {"Annotation_mode": "split", "Gene_name": "OPA1", "ACMG_class": "5", "AnnotSV ranking score": "1.0"}
        view = annotsv_row_for_gene([full, split], "OPA1")
        self.assertEqual(view["ACMG_class"], "3")
        self.assertEqual(view["AnnotSV ranking score"], "0.3")
        self.assertEqual(view["_CLASSIFICATION_SOURCE"], "FULL_EVENT_ROW")
        conflict = annotsv_row_for_gene([full, {**full, "ACMG_class": "5"}, split], "OPA1")
        self.assertEqual(conflict["ACMG_class"], ".")
        self.assertEqual(conflict["_CLASSIFICATION_CONFLICT"], "YES")

    def test_patient_context_and_unannotated_rows_are_retained(self):
        evidence = pd.DataFrame([clinical_evidence(), clinical_evidence(SV_ID="unannotated", GENES=".", HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED=".")])
        context = pd.DataFrame([{"SV_ID": "positive", "GENE": "OPA1", "HPO_QUERY_RESNIK_NORMALIZED": "0.1", "HPO_SEMANTIC_STATUS": "EVALUATED", "ACMG_CNV_CLASS": "BENIGN"}])
        out = build_clinical_table(evidence, context).set_index("SV_ID")
        self.assertEqual(len(out), 2)
        self.assertEqual(out.loc["positive", "CLINICAL_HON_QUERY_SCOPE"], "PATIENT_OBSERVED_HPO_QUERY")
        self.assertEqual(out.loc["positive", "CLINICAL_CNV_CLASS"], "LIKELY_PATHOGENIC")
        self.assertEqual(out.loc["unannotated", "GENE"], ".")

    def test_integrator_carries_core_hon_and_event_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "input.vcf").write_text(
                "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr3\t1000\tevent\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=3000;SVLEN=-2000;SUPP=1;IDLIST=0_Sniffles2.DEL.1\n")
            pd.DataFrame([
                {"SV_ID": "event", "Annotation_mode": "full", "Gene_name": "OPA1",
                 "ACMG_class": "4", "AnnotSV ranking score": "0.90", "AnnotSV ranking criteria": "2C-1"},
                {"SV_ID": "event", "Annotation_mode": "split", "Gene_name": "OPA1", "ACMG_class": "5",
                 "Overlapped_CDS_percent": "25"},
            ]).to_csv(folder / "ann.tsv", sep="\t", index=False)
            pd.DataFrame([{"SV_ID": "Sniffles2.DEL.1", "CALLER": "Sniffles2", "CHROM": "chr3",
                           "START": "1000", "END": "3000", "SVTYPE": "DEL", "SVLEN": "-2000",
                           "CALLER_SUPPORT": "8", "CALLER_GT": "0/1", "CALLER_GQ": "30", "CALLER_DP": "20",
                           "EVIDENCE_STATUS": "PASS", "EVIDENCE_FLAGS": "."}]).to_csv(folder / "caller.tsv", sep="\t", index=False)
            pd.DataFrame([{"gene": "OPA1", "HON_CORE_SEMANTIC_SIMILARITY_NORMALIZED": "0.8",
                           "HON_MITO_SYNDROMIC_CONTEXT_SIMILARITY_NORMALIZED": "0.1",
                           "HON_SEMANTIC_SIMILARITY_NORMALIZED": "0.5"}]).to_csv(folder / "rank.tsv", sep="\t", index=False)
            subprocess.run([sys.executable, str(ROOT / "scripts/build_integrated_sv_gene_tsv.py"),
                            "--sequencing-type", "LRS", "--jasmine-id-prefixes",
                            "--vcf", str(folder / "input.vcf"), "--annotsv", str(folder / "ann.tsv"),
                            "--caller-tsv", str(folder / "caller.tsv"),
                            "--ranking", str(folder / "rank.tsv"), "--output", str(folder / "integrated.tsv")],
                           check=True, capture_output=True)
            evidence = pd.read_csv(folder / "integrated.tsv", sep="\t", dtype=str)
            out = build_clinical_table(evidence).iloc[0]
            self.assertEqual(out.CLINICAL_CNV_CLASS, "LIKELY_PATHOGENIC")
            self.assertEqual(out.ANNOTSV_CLASSIFICATION_SOURCE, "FULL_EVENT_ROW")
            self.assertEqual(out.CLINICAL_HON_QUERY_SIMILARITY, 0.8)
            self.assertEqual(out.ANNOTSV_RANKING_CRITERIA, "2C-1")
            pd.DataFrame([{"gene": "OPA1", "disease_id": "SYNTHETIC:HON", "moi": "AD",
                           "mechanism": "LOF", "validity": "STRONG", "source": "synthetic_test_fixture",
                           "hon_relevance": "SUPPORTS_HON", "hon_relevance_source": "synthetic_test_fixture"}]).to_csv(folder / "models.tsv", sep="\t", index=False)
            pd.DataFrame([{"sample_id": "SYNTHETIC", "SV_ID": "event", "gene": "OPA1", "disease_id": ".",
                           "population_af": "0.2", "population_match": "EXACT_ALLELE",
                           "population_source": "synthetic_population_release", "source": "synthetic_test_fixture"}]).to_csv(folder / "allele_evidence.tsv", sep="\t", index=False)
            subprocess.run([sys.executable, str(ROOT / "scripts/assess_sv_alleles.py"),
                            "--integrated", str(folder / "integrated.tsv"), "--sample", "SYNTHETIC",
                            "--disease-models", str(folder / "models.tsv"), "--allele-evidence", str(folder / "allele_evidence.tsv"),
                            "--output", str(folder / "assessed.tsv"), "--hypotheses-output", str(folder / "hypotheses.tsv"),
                            "--manifest", str(folder / "manifest.json")], check=True, capture_output=True)
            subprocess.run([sys.executable, str(ROOT / "scripts/build_diagnostic_candidates.py"),
                            "--input", str(folder / "assessed.tsv"), "--output", str(folder / "research.tsv"),
                            "--clinical-output", str(folder / "clinical.tsv")], check=True, capture_output=True)
            out = pd.read_csv(folder / "clinical.tsv", sep="\t", dtype=str).iloc[0]
            self.assertEqual(out.ALLELE_EFFECT, "POSSIBLE_LOF")
            self.assertEqual(out.ALLELE_GENOTYPE, "HET")
            self.assertEqual(out.CLINICAL_MECHANISM_SUPPORTED, "YES")
            self.assertEqual(out.CLINICAL_DISEASE_HON_STATUS, "SUPPORTS_HON")
            self.assertEqual(float(out.CLINICAL_MAX_REPORTED_AF), 0.2)
            self.assertEqual(out.CLINICAL_REVIEW_GROUP, "CONFLICT_OR_TECHNICAL_REVIEW")

    def test_empty_table_and_duplicate_keys(self):
        empty = build_clinical_table(pd.DataFrame(columns=["SV_ID", "GENES"]))
        self.assertIn("CLINICAL_REVIEW_RANK", empty)
        self.assertEqual(len(empty), 0)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            build_clinical_table(pd.DataFrame([clinical_evidence(), clinical_evidence()]))

    def test_cli_keeps_original_srs_review_and_clinical_output_is_unlimited(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            source = folder / "input.tsv"
            pd.DataFrame([clinical_evidence(), clinical_evidence(SV_ID="vus", ACMG_CNV_CLASS="VUS", ACMG_CNV_SCORE="0.3")]).to_csv(source, sep="\t", index=False)
            base = [sys.executable, str(ROOT / "scripts/build_diagnostic_candidates.py"), "--input", str(source), "--top-n", "1"]
            subprocess.run([*base, "--output", str(folder / "original.tsv")], check=True, capture_output=True)
            subprocess.run([*base, "--output", str(folder / "with_clinical.tsv"), "--clinical-output", str(folder / "clinical.tsv")], check=True, capture_output=True)
            self.assertEqual((folder / "original.tsv").read_bytes(), (folder / "with_clinical.tsv").read_bytes())
            self.assertEqual(len(pd.read_csv(folder / "clinical.tsv", sep="\t")), 2)

    def test_lrs_finalizer_keeps_original_ranks_and_retains_unannotated_events(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            sv = pd.DataFrame([clinical_evidence(), clinical_evidence(SV_ID="vus", ACMG_CNV_CLASS="VUS", ACMG_CNV_SCORE="0.3")])
            sv["GENE"] = sv["GENES"]
            sv.to_csv(folder / "sv.tsv", sep="\t", index=False)
            genes = pd.DataFrame([{"GENE": "OPA1", "PANEL_STATUS": "PANEL_GENE", "GENE_RELEVANCE_TIER": "HIGH", "GENE_RELEVANCE_DISPLAY_SCORE": "9"}])
            genes.to_csv(folder / "genes.tsv", sep="\t", index=False)
            pd.DataFrame(columns=["GENE", "SV_ID", "RECESSIVE_PAIR_STATUS"]).to_csv(folder / "pairs.tsv", sep="\t", index=False)
            evidence = pd.concat([sv, pd.DataFrame([clinical_evidence(SV_ID="unannotated", GENES=".")])], ignore_index=True)
            evidence.to_csv(folder / "evidence.tsv", sep="\t", index=False)
            base = [sys.executable, str(ROOT / "scripts/finalize_candidate_ranking.py"), "--sample", "demo", "--genes", str(folder / "genes.tsv"), "--sv", str(folder / "sv.tsv"), "--snv-sv", str(folder / "pairs.tsv")]
            subprocess.run([*base, "--gene-output", str(folder / "gene_orig.tsv"), "--sv-output", str(folder / "sv_orig.tsv")], check=True, capture_output=True)
            subprocess.run([*base, "--gene-output", str(folder / "gene_new.tsv"), "--sv-output", str(folder / "sv_new.tsv"), "--clinical-evidence", str(folder / "evidence.tsv"), "--clinical-output", str(folder / "clinical.tsv")], check=True, capture_output=True)
            self.assertEqual((folder / "gene_orig.tsv").read_bytes(), (folder / "gene_new.tsv").read_bytes())
            self.assertEqual((folder / "sv_orig.tsv").read_bytes(), (folder / "sv_new.tsv").read_bytes())
            self.assertEqual(len(pd.read_csv(folder / "clinical.tsv", sep="\t")), 3)


class TestFinalCandidateTables(unittest.TestCase):
    def test_gene_relevance_does_not_use_sv_count(self):
        ranking = pd.DataFrame([
            {
                "gene": "OPA1",
                "panel_gene": "YES",
                "human_HPO_count": "12",
                "optic_neuropathy_anchor_HPO_count": "2",
                "phenotype_score": "6.0",
                "gene_disease_evidence_score": "4.0",
                "gene_disease_evidence_level": "DEFINITIVE",
                "SV_count": "99",
                "SV_types": "DEL;INS",
                "SV_count_lt100kb": "50",
                "SV_count_100kb_to_1Mb": "20",
                "SV_count_1Mb_to_10Mb": "20",
                "SV_count_ge10Mb": "9",
                "breakpoint_defined_INV_BND_count": "0",
                "candidate_group": "PANEL_GENE",
            }
        ])
        out = build_gene_table(ranking, pd.DataFrame())
        self.assertEqual(float(out.iloc[0]["GENE_RELEVANCE_SCORE"]), 10.0)
        self.assertEqual(out.iloc[0]["SV_COUNT"], "99")


if __name__ == "__main__":
    unittest.main()
