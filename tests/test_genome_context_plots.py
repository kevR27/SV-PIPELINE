"""Coordinate/counting/phase checks using artificial SVs, not patient results."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "plots"))
sys.path.insert(0,str(ROOT / "scripts"))
import plot_genome_wide_sv_map as genome
import plot_event_review_summary as counts
from candidate_review import event_frequency_table
from plot_utils import abs_svlen
from build_snv_sv_candidates import pair_phase_relation, read_phase
from plot_candidate_locus import load_samplot_images, plot_assembly_alignment
import run_thesis_plots as launcher


def events():
    return pd.DataFrame([
        {"SV_ID":"del1","GENE":"TEST_A","CHROM":"1","START":101,"END":501,"SVTYPE":"DEL","SVLEN":-400,"CALLER_COUNT":2,"PANEL_STATUS":"PANEL_GENE","NEEDLR_AF":"0.01"},
        {"SV_ID":"del1","GENE":"TEST_B","CHROM":"1","START":101,"END":501,"SVTYPE":"DEL","SVLEN":-400,"CALLER_COUNT":2,"PANEL_STATUS":"NONPANEL_GENE","NEEDLR_AF":"0.01"},
        {"SV_ID":"ins1","GENE":"TEST_C","CHROM":"chr1","START":101,"END":101,"SVTYPE":"INS","SVLEN":70,"CALLER_COUNT":1,"PANEL_STATUS":"NONPANEL_GENE","NEEDLR_AF":"."},
        {"SV_ID":"bnd1","GENE":"TEST_D","CHROM":"chr2","START":1001,"END":201,"SVTYPE":"TRA","CHR2":"1","POS2":201,"CALLER_COUNT":2,"PANEL_STATUS":"PANEL_GENE","NEEDLR_AF":"0.02"},
    ])


class GenomeContextTests(unittest.TestCase):
    def test_breakend_coordinates_are_not_treated_as_sv_length(self):
        frame = pd.DataFrame([{"SVTYPE":"BND","START":100,"END":5000,"SVLEN":4900},
                              {"SVTYPE":"INS","START":100,"END":5000,"SVLEN":"."},
                              {"SVTYPE":"DEL","START":100,"END":5000,"SVLEN":"."}])
        sizes = abs_svlen(frame)
        self.assertTrue(pd.isna(sizes.iloc[0]))
        self.assertTrue(pd.isna(sizes.iloc[1]))
        self.assertEqual(sizes.iloc[2],4900)

    def test_depth_bed_is_not_discovered_as_methylation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample = root/"SYNTHETIC"
            (sample/"coverage/binned").mkdir(parents=True)
            (sample/"coverage/binned/SYNTHETIC.10kb.regions.bed.gz").touch()
            source = root/"master.tsv"
            events().to_csv(source,sep="\t",index=False)
            args = ["launcher","--root",str(root),"--sample","SYNTHETIC","--evidence-table",str(source),"--candidate-table",str(source),"--dry-run"]
            with patch.object(sys,"argv",args),patch.object(launcher,"run") as run:
                launcher.main()
            scripts = [Path(call.args[0][1]).name for call in run.call_args_list]
            self.assertNotIn("plot_methylation.py",scripts)

    def test_master_ids_once_and_breakend_uses_real_partner(self):
        frame = genome.prepare_events(events(),{"chr1":10000,"chr2":20000})
        self.assertEqual(len(frame),3)
        bnd = frame.set_index("SV_ID").loc["bnd1"]
        self.assertEqual(bnd["SVTYPE"],"BND")
        self.assertEqual(genome.endpoints(bnd),[("chr2",1000),("chr1",200)])

    def test_conflicting_propagation_is_not_arbitrarily_plotted(self):
        raw = events()
        raw.loc[1,"START"] = 102
        frame = genome.prepare_events(raw,{"chr1":10000,"chr2":20000})
        self.assertEqual(frame.set_index("SV_ID").at["del1","MAP_STATUS"],"CONFLICTING_COORDINATES")

    def test_reference_bounds_and_missing_partner_stay_explicit(self):
        raw = events()
        raw.loc[0:1,"START"] = 0
        raw.loc[3,"POS2"] = np.nan
        frame = genome.prepare_events(raw,{"chr1":10000,"chr2":20000}).set_index("SV_ID")
        self.assertEqual(frame.at["del1","MAP_STATUS"],"INVALID_START")
        self.assertEqual(frame.at["bnd1","PARTNER_STATUS"],"UNRESOLVED_OR_OFF_MAP")

    def test_rainfall_keeps_ties_and_never_crosses_chromosomes(self):
        frame = genome.prepare_events(events(),{"chr1":10000,"chr2":20000})
        rainfall = genome.rainfall_table(frame)
        self.assertEqual(len(rainfall),3)
        self.assertEqual(rainfall["DISTANCE_BP"].dropna().tolist(),[0.])
        self.assertTrue(rainfall[rainfall["CHROM"].eq("chr2")]["DISTANCE_BP"].isna().all())

    def test_bed_half_open_boundary(self):
        idx = genome.interval_index(pd.DataFrame([("chr1",100,200)],columns=["CHROM","START0","END0"]))
        self.assertFalse(genome.overlaps(idx,"chr1",99,100))
        self.assertTrue(genome.overlaps(idx,"chr1",100,101))
        self.assertFalse(genome.overlaps(idx,"chr1",200,201))

    def test_wrong_resource_length_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp)/"wrong.bed"
            bed.write_text("chr1\t0\t10001\n")
            with self.assertRaises(ValueError):
                genome.read_bed(bed,{"chr1":10000})

    def test_measured_depth_weighting_and_gap_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            depth = Path(tmp)/"depth.bed"
            depth.write_text("chr1\t0\t100\t0\nchr1\t100\t200\t10\nchr1\t200\t400\t20\nchrX\t0\t400\t5\n")
            mask = {"chr1":np.array([[0,100]])}
            frame = genome.load_depth(depth,{"chr1":400,"chrX":400},mask,display_bin=400)
            self.assertAlmostEqual(frame[frame["CHROM"].eq("chr1")].iloc[0]["MEAN_DEPTH"],50/3)
            self.assertTrue(np.allclose(frame["RELATIVE_DEPTH"],1))
            self.assertEqual(frame[frame["CHROM"].eq("chrX")].iloc[0]["BASELINE_SCOPE"],"OWN_CHROMOSOME_MEAN")

    def test_repeat_expectation_has_reproducible_evaluable_denominator(self):
        raw = events().iloc[:2]
        frame = genome.prepare_events(raw,{"chr1":10000})
        bed = pd.DataFrame([("chr1",0,10000,"repeat")],columns=["CHROM","START0","END0","LABEL"])
        with tempfile.TemporaryDirectory() as tmp, patch.object(genome,"save_figure"):
            genome.repeat_context(frame,{"Tandem":bed},{"chr1":10000},{},str(Path(tmp)/"repeat"),repeats=20,seed=5)
            result = pd.read_csv(Path(tmp)/"repeat_summary.tsv",sep="\t")
            self.assertEqual(result.at[0,"N_EVALUABLE"],1)
            self.assertEqual(result.at[0,"OBSERVED_PERCENT"],100)
            self.assertEqual(result.at[0,"RANDOM_MEAN_PERCENT"],100)

    def test_frequency_combines_unknown_with_boundary_but_not_common_or_invalid(self):
        raw = events()
        raw = pd.concat([raw,pd.DataFrame([{**raw.iloc[0].to_dict(),"SV_ID":"invalid","NEEDLR_AF":"105"}])],ignore_index=True)
        frequency = event_frequency_table(raw).set_index("SV_ID")
        self.assertEqual(frequency.at["ins1","AF_STATUS"],"AF_UNKNOWN")
        self.assertEqual(frequency.at["del1","AF_STATUS"],"MEASURED_AT_OR_BELOW_CUTOFF")
        membership = counts.evidence_membership(raw,pd.DataFrame({"SV_ID":["del1"]})).set_index("SV_ID")
        self.assertTrue(membership.at["ins1","AF ≤1% or unknown"])
        self.assertFalse(membership.at["bnd1","AF ≤1% or unknown"])
        self.assertFalse(membership.at["invalid","AF ≤1% or unknown"])
        self.assertTrue(membership["Patient HPO"].isna().all())

    def test_common_af_from_another_gene_row_is_not_lost(self):
        raw = events()
        raw.loc[1,"NEEDLR_AF"] = .1
        self.assertEqual(event_frequency_table(raw).set_index("SV_ID").at["del1","MAX_POPULATION_AF"],.1)

    def test_gene_count_uses_actual_coordinates_not_a_nearby_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp)/"genes.bed"
            bed.write_text("chr1\t300\t400\tTEST_A\nchr2\t5000\t6000\tTEST_D\n")
            membership = counts.evidence_membership(events(),pd.DataFrame({"SV_ID":["del1"]}),str(bed)).set_index("SV_ID")
            self.assertTrue(membership.at["del1","Gene overlap"])
            self.assertFalse(membership.at["ins1","Gene overlap"])
            self.assertFalse(membership.at["bnd1","Gene overlap"])

    def test_phase_labels_cannot_cross_chromosomes_or_blocks(self):
        small = {"CHROM":"1","SMALL_GT":"0|1","SMALL_PS":"100","SMALL_PHASED":"YES"}
        sv = {"CHROM":"chr1","LONGPHASE_GT":"1|0","LONGPHASE_PS":"100","LONGPHASE_PHASED":"YES"}
        self.assertEqual(pair_phase_relation(small,sv),"TRANS")
        self.assertEqual(pair_phase_relation(small,{**sv,"CHROM":"chr2"}),"UNRESOLVED")
        self.assertEqual(pair_phase_relation(small,{**sv,"LONGPHASE_PS":"101"}),"UNRESOLVED")
        self.assertEqual(pair_phase_relation(small,{**sv,"LONGPHASE_GT":"1|1"}),"UNRESOLVED")

    def test_multiallelic_record_does_not_assign_one_haplotype_to_every_alt(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = Path(tmp)/"phase.vcf"
            vcf.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\nchr1\t100\t.\tA\tC,G\t.\tPASS\t.\tGT:PS\t0|1:100\n")
            phase,_ = read_phase(str(vcf))
            self.assertEqual(phase[("1",100,"G")][2],"NO")

    def test_failed_samplot_is_not_reused_from_a_stale_image(self):
        with tempfile.TemporaryDirectory() as tmp:
            image = Path(tmp)/"stale.png"
            image.write_bytes(b"unused")
            manifest = Path(tmp)/"manifest.tsv"
            pd.DataFrame([{"SV_ID":"del1","STATUS":"FAILED","IMAGE":str(image)}]).to_csv(manifest,sep="\t",index=False)
            self.assertEqual(load_samplot_images([manifest]),{})

    def test_assembly_alignment_source_preserves_reverse_strand(self):
        with tempfile.TemporaryDirectory() as tmp, patch("plot_candidate_locus.save_figure"):
            paf = Path(tmp)/"assembly.paf"
            paf.write_text("ctg\t1000\t10\t400\t-\tchr1\t10000\t100\t500\t370\t400\t60\n")
            prefix = Path(tmp)/"alignment"
            plot_assembly_alignment(paf,prefix)
            source = pd.read_csv(str(prefix)+"_source.tsv",sep="\t")
            self.assertEqual(source.at[0,"STRAND"],"-")
            self.assertEqual(source.at[0,"REF_START0"],100)


if __name__ == "__main__":
    unittest.main()
