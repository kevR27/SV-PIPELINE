"""Regressions found while checking saved outputs. All records here are synthetic."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import annotsv_reconciliation as reconciliation
import assess_sv_alleles as allele
import audit_annotsv_evidence as availability
import build_integrated_sv_gene_tsv as integration
from check_annotsv_gene_resource import check_resource
from sv_evidence_common import gene_symbols, invalid_gene_labels, sv_length


def table(path, rows, fields=None):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fields or list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


class SavedOutputTests(unittest.TestCase):
    def run_script(self, name, *args, success=True):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / name), *map(str, args)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result

    def test_length_fallback_and_insertion_endpoint(self):
        self.assertEqual(sv_length({"SVTYPE":"DEL", "REF":"A"*36, "ALT":"A"}), ("-35", "REF_ALT_LENGTH"))
        self.assertEqual(sv_length({"SVTYPE":"DEL", "REF":"N", "ALT":"<DEL>", "START":"100", "END":"135"}), ("-35", "END_MINUS_POS"))
        self.assertEqual(sv_length({"SVTYPE":"INS", "START":"100", "END":"100", "ALT":"<INS>"}), (".", "UNKNOWN"))
        self.assertEqual(sv_length({"SVTYPE":"INS", "SVLEN":"600", "REF":"N", "ALT":"N"}), ("600", "INFO_SVLEN"))
        self.assertEqual(sv_length({"SVTYPE":"BND", "SVLEN":"0", "START":"100", "END":"9999"}), (".", "NOT_APPLICABLE_BREAKEND"))

    def test_size_filter_catches_delly_missing_svlen_and_keeps_bnd(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            base = {"SV_ID":"small", "CHROM":"chr1", "START":"100", "END":"135", "SVTYPE":"DEL",
                    "SVLEN":".", "REF":"A"*36, "ALT":"A", "FILTER":"PASS", "CALLER_SUPPORT":"6"}
            table(d/"in.tsv", [base, {**base,"SV_ID":"bnd","SVTYPE":"BND","SVLEN":"0","ALT":"A[chr2:200["}])
            self.run_script("filter_sv_evidence.py", "--caller-tsv", d/"in.tsv", "--caller", "Delly",
                            "--output-tsv", d/"out.tsv", "--output-ids", d/"ids")
            a,b=read(d/"out.tsv")
            self.assertEqual(a["EVIDENCE_FAIL_REASONS"], "SIZE_BELOW_MIN")
            self.assertEqual(a["SVLEN_SOURCE"], "REF_ALT_LENGTH")
            self.assertEqual(b["EVIDENCE_STATUS"], "PASS")
            self.assertEqual((d/"ids").read_text(), "bnd\n")

    def test_blacklist_uses_both_bnd_points_not_cross_chromosome_span(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            r={"SV_ID":"b", "CHROM":"chr1", "START":"100", "END":"900", "ALT":"N[chr2:900[",
               "SVTYPE":"BND", "SVLEN":".", "FILTER":"PASS", "CALLER_SUPPORT":"8"}
            table(d/"in.tsv", [r])
            for bed, expected in [("chr1\t400\t500\n", "."), ("chr2\t899\t900\n", "BLACKLIST_REGION")]:
                (d/"blacklist.bed").write_text(bed)
                self.run_script("filter_sv_evidence.py", "--caller-tsv",d/"in.tsv","--caller","Delly",
                                "--blacklist-bed",d/"blacklist.bed","--output-tsv",d/"out.tsv","--output-ids",d/"ids")
                self.assertEqual("BLACKLIST_REGION" in read(d/"out.tsv")[0]["EVIDENCE_FLAGS"], expected != ".")

    def test_transcript_status_labels_are_not_genes_but_unk_is(self):
        self.assertEqual(gene_symbols("cmpl;incmpl;UNK;OPA1;NA"), ["OPA1","UNK"])
        self.assertEqual(invalid_gene_labels("cmpl;incmpl;UNK"), ["CMPL","INCMPL"])

    def test_gene_resource_check_rejects_status_labels_and_accepts_unk(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"genes.bed"
            p.write_text("1\t100\t200\t+\tUNK\tENST_TEST\t120\t180\t100,\t200,\n")
            self.assertEqual(check_resource(p),1)
            p.write_text(p.read_text().replace("UNK", "cmpl"))
            with self.assertRaisesRegex(ValueError,"needs repair"):
                check_resource(p)

    def test_skip_log_explains_missing_and_rejects_wrong_vcf(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"skip.tsv"
            p.write_text("1_100_135_DEL_AA_A: variantLength (35) < SVminSize (50) (line 3)\n")
            sv={"SV_ID":"s", "_VCF_LINE":3,"CHROM":"chr1","START":"100","END":"135","SVTYPE":"DEL","SVLEN":"-35","REF":"AA","ALT":"A"}
            notes=reconciliation.read_skip_log(p,[sv])
            r=reconciliation.record_evidence(sv,[],notes,True)
            self.assertEqual(r["ANNOTSV_RECORD_STATUS"],"BELOW_ANNOTSV_MIN_SIZE")
            with self.assertRaisesRegex(ValueError,"does not match"):
                reconciliation.read_skip_log(p,[{**sv,"START":"101"}])

    def test_reciprocal_breakend_notice_does_not_mean_record_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"skip.tsv";p.write_text("2_900_900_TRA_1: reciprocal breakend (line 3)\n")
            sv={"SV_ID":"b","_VCF_LINE":3,"CHROM":"chr1","START":"100","SVTYPE":"BND","REF":"N","ALT":"N[chr2:900["}
            notes=reconciliation.read_skip_log(p,[sv])
            self.assertEqual(reconciliation.record_evidence(sv,[{"Gene_name":"A"}],notes,True)["ANNOTSV_RECORD_STATUS"],"ANNOTATED_WITH_LOG_NOTICE")

    def test_saved_output_audit_does_not_claim_bundle_is_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"a.tsv";table(p,[{"P_loss_source":"CLN:1","B_loss_source":"DGV:1"}])
            audit={"sample":availability.sample_columns(p),"bundle":availability.bundle_inventory(None,"GRCh38")}
            statuses=availability.evidence_availability({},"DEL",audit)
            self.assertEqual(statuses["SV_PATHOGENIC_DB_STATUS"],"NO_OVERLAP_REPORTED")
            self.assertEqual(statuses["SV_DB_CLINVAR_AVAILABILITY"],"DATABASE_OBSERVED_AVAILABILITY_INCOMPLETE")
            self.assertEqual(availability.evidence_availability({},"INV",audit)["SV_PATHOGENIC_DB_STATUS"],"SOURCE_COLUMN_MISSING")

    def test_bnd_remote_gene_disruption_uses_remote_coordinate(self):
        row={"SVTYPE":"BND","CHROM":"chr1","START":"100","CHR2":"chr2","POS2":"900"}
        ann={"SV_chrom":"2","Tx_start":"800","Tx_end":"1000","Location":"exon2"}
        row["ANNOTSV_GENE_ROWS_JSON"]=json.dumps([ann])
        self.assertEqual(allele.disruption(row,{})[0][0],"BREAKPOINT_DISRUPTION_POSSIBLE")
        row["ANNOTSV_GENE_ROWS_JSON"]=json.dumps([{**ann,"Tx_start":"90","Tx_end":"110"}])
        self.assertEqual(allele.disruption(row,{})[0][0],"UNKNOWN")

    def test_contained_inversion_gene_has_no_disruption_points(self):
        row={"SVTYPE":"INV","START":"100","END":"1000", "ANNOTSV_GENE_ROWS_JSON":json.dumps([
             {"Tx_start":"200","Tx_end":"800","Overlapped_CDS_percent":"100","Location":"txStart-txEnd"}])}
        self.assertEqual(allele.disruption(row,{})[0][1],0)
        self.assertIn("NO_INTRAGENIC_BREAKPOINT",allele.functional_context(row))

    def test_current_filter_failure_is_visible_in_technical_assessment(self):
        row={"SVTYPE":"DEL","CALLER_COUNT":"2","CALLER_EVIDENCE_MATCH":"IDLIST","FILTER":"PASS",
             "CALLER_EVIDENCE_JSON":json.dumps([{"CALLER_SUPPORT":"9","EVIDENCE_STATUS":"FAIL","EVIDENCE_FAIL_REASONS":"SIZE_BELOW_MIN"}])}
        result=allele.technical(row,5)
        self.assertEqual(result[0],"REVIEW_REQUIRED")
        self.assertIn("SIZE_BELOW_MIN", result[2])

    def test_rebuild_saved_outputs_preserves_sv_and_quarantines_bad_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/"master.vcf").write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\tv\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=1000;SVLEN=900;SUPP_VEC=100\n")
            table(d/"ann.tsv",[{"ID":"v","AnnotSV_ID":"one","Annotation_mode":"full","Gene_name":"A;cmpl"},
                               {"ID":"v","AnnotSV_ID":"one","Annotation_mode":"split","Gene_name":"cmpl"}])
            (d/"panel").write_text("A\n")
            args=["--sample","TEST","--vcf",d/"master.vcf","--annotsv",d/"ann.tsv","--panel",d/"panel","--outdir",d/"out"]
            self.run_script("rebuild_sv_results.py",*args)
            rows=read(d/"out/TEST_allele_assessment.tsv")
            self.assertEqual({r["SV_ID"] for r in rows},{"v"})
            self.assertEqual({r["GENES"] for r in rows},{"A","."})
            self.assertTrue(all(r["ALLELE_PHENOTYPE_STATUS"]=="UNKNOWN_PATIENT_PHENOTYPES" for r in rows))
            self.assertTrue(all(r["NEEDLR_STATUS"]=="NOT_SUPPLIED" for r in rows))
            unresolved=next(r for r in rows if r["GENES"]==".")
            self.assertEqual(json.loads(unresolved["ANNOTSV_UNRESOLVED_GENE_ROWS_JSON"])[0]["Gene_name"],"cmpl")
            self.run_script("rebuild_sv_results.py",*args,success=False)

    def test_rank_counts_master_id_once_across_breakpoint_annotations(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/"genes").write_text("A\n");(d/"panel").write_text("A\n")
            table(d/"hpo",[],["gene_symbol","hpo_id","optic_neuropathy_anchor"])
            table(
                d/"hon_seeds.tsv",
                [{
                    "hpo_id":"HP:0000648",
                    "hpo_label":"Optic atrophy",
                    "hon_seed_role":"CORE_OCULAR_HON",
                }],
            )
            table(
                d/"edges.tsv",
                [],
                ["subject","object","predicate","category"],
            )
            table(d/"ann",[{"ID":"v","AnnotSV_ID":endpoint,"Gene_name":"A","Annotation_mode":"split"} for endpoint in ["left","right"]])
            self.run_script(
                "rank_sv_gene_candidates.py",
                "--annotsv",d/"ann","--genes",d/"genes",
                "--phenotypes",d/"hpo","--panel",d/"panel",
                "--hpo-seeds",d/"hon_seeds.tsv","--edges",d/"edges.tsv",
                "--output",d/"rank",
            )
            r=read(d/"rank")[0]
            self.assertEqual(r["SV_count"],"1")

    def test_vep_audit_separates_nearby_genes_from_unannotated_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            table(d/"ann",[{"ID":"a","Gene_name":"A"}])
            table(d/"records",[{"SV_ID":"a","ANNOTSV_RECORD_STATUS":"ANNOTATED"},
                               {"SV_ID":"b","ANNOTSV_RECORD_STATUS":"BELOW_ANNOTSV_MIN_SIZE"}])
            (d/"vep").write_text("## VEP command-line: vep --pick --symbol\n#Uploaded_variation\tFeature\tConsequence\tExtra\na\tTX1\tdownstream_gene_variant\tSYMBOL=B\nb\tTX2\tintron_variant\tSYMBOL=C\n")
            self.run_script("audit_vep_gene_coverage.py","--vep",d/"vep","--annotsv",d/"ann","--records",d/"records",
                            "--output-tsv",d/"out.tsv","--output-json",d/"out.json")
            self.assertEqual([r["RELATION"] for r in read(d/"out.tsv")],["NEARBY_GENE","ANNOTSV_RECORD_NOT_ANNOTATED"])
            self.assertEqual(json.loads((d/"out.json").read_text())["transcript_selection"],"PICK_RESTRICTED")


if __name__ == "__main__":
    unittest.main()
