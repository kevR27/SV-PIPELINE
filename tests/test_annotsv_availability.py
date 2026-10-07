import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import audit_annotsv_evidence as audit


def tsv(path, fields, records):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(records)


class AvailabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.resources = audit.resource_paths(self.root)
        self.sample = self.root/'sample.tsv'

    def tearDown(self):
        self.temp.cleanup()

    def write_resource(self, field, content):
        p = self.resources[field]
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content)
        return p

    def make_audit(self, fields=('P_loss_source','B_loss_source'), records=None):
        tsv(self.sample,fields,records or [])
        return {'version':audit.VERSION, 'bundle':audit.bundle_inventory(self.root,'GRCh38'), 'sample':audit.sample_columns(self.sample)}

    def test_missing_column_distinct_from_no_reported_overlap(self):
        self.write_resource('P_loss_source','1\t1\t50\tDisease\tHP:0000001\tCLN:1\t1:1-50\n')
        self.write_resource('B_loss_source','1\t100\t150\tDGV:1\t1:100-150\t0.1\n')
        report=self.make_audit()
        result=audit.evidence_availability({},'DEL',report)
        self.assertEqual(result['SV_DB_CLINVAR_AVAILABILITY'],'NO_OVERLAP_REPORTED')
        self.assertEqual(result['SV_DB_GNOMAD_AVAILABILITY'],'DATABASE_NOT_REPRESENTED_IN_PROCESSED_RESOURCES')
        report=self.make_audit(fields=['B_loss_source'])
        self.assertEqual(audit.evidence_availability({},'DEL',report)['SV_PATHOGENIC_DB_STATUS'],'SOURCE_COLUMN_MISSING')

    def test_missing_empty_and_unrecognized_resources_distinct(self):
        report=self.make_audit()
        self.assertEqual(audit.evidence_availability({},'DEL',report)['SV_PATHOGENIC_DB_STATUS'],'PROCESSED_RESOURCE_NOT_FOUND')
        self.write_resource('P_loss_source','')
        report=self.make_audit()
        self.assertEqual(audit.evidence_availability({},'DEL',report)['SV_PATHOGENIC_DB_STATUS'],'PROCESSED_RESOURCE_EMPTY')
        self.write_resource('P_loss_source','unknown\tformat\n')
        report=self.make_audit()
        self.assertEqual(audit.evidence_availability({},'DEL',report)['SV_PATHOGENIC_DB_STATUS'],'RESOURCE_SCHEMA_UNRECOGNIZED')

    def test_output_positive_survives_missing_local_resource(self):
        report=self.make_audit()
        result=audit.evidence_availability({'P_loss_source':'CLN:123'},'DEL',report)
        self.assertEqual(result['SV_DB_CLINVAR_AVAILABILITY'],'OVERLAP_REPORTED')
        self.assertEqual(result['SV_PATHOGENIC_DB_STATUS'],'OVERLAP_REPORTED')

    def test_missing_annotation_and_bnd_not_negative(self):
        report=self.make_audit()
        self.assertEqual(audit.evidence_availability({},'DEL',report,False)['SV_DB_DGV_AVAILABILITY'],'NO_ANNOTSV_MATCH')
        self.assertEqual(audit.evidence_availability({},'BND',report)['SV_DB_DGV_AVAILABILITY'],'NOT_APPLICABLE')
        self.assertEqual(audit.evidence_availability({},'INV',None)['SV_PATHOGENIC_DB_STATUS'],'AVAILABILITY_NOT_AUDITED')

    def test_scan_correct_source_column_and_database_tokens(self):
        path=self.write_resource('P_inv_source','1\t1\t100\tDGV phenotype word\t.\tCLN:1\t1:1-100\n')
        result=audit.scan_resource('P_inv_source',path)
        self.assertEqual(result['database_records']['CLINVAR'],1)
        self.assertEqual(result['database_records']['DGV'],0)
        self.assertTrue(audit.mentions('HI3:GENE','CLINGEN'))
        self.assertFalse(audit.mentions('NOTCLN','CLINVAR'))

    def test_cache_reused_then_invalidated(self):
        path=self.write_resource('B_inv_source','1\t1\t2\tDGV:1\t1:1-2\t.1\n')
        cache=self.root/'cache.json'
        audit.bundle_inventory(self.root,'GRCh38',cache)
        with patch.object(audit,'scan_resource',side_effect=AssertionError('cache should avoid rescanning')):
            audit.bundle_inventory(self.root,'GRCh38',cache)
        path.write_text(path.read_text()+'1\t3\t4\tgnomAD:2\t1:3-4\t.2\n')
        result=audit.bundle_inventory(self.root,'GRCh38',cache)
        self.assertEqual(result['channels']['B_inv_source']['database_records']['GNOMAD'],1)

    def test_integration_keeps_legacy_flag_and_adds_reason(self):
        self.write_resource('P_loss_source','1\t1000\t2000\tDisease\t.\tCLN:1\t1:1000-2000\n')
        self.write_resource('B_loss_source','')
        fields=['SV_ID','Gene_name','Annotation_mode','P_loss_source','B_loss_source']
        report=self.make_audit(fields,[dict(zip(fields,['v','A','full','.','.']))])
        audit.atomic_json(self.root/'audit.json',report)
        vcf=self.root/'in.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\tv\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200;SVLEN=-100\n')
        subprocess.run([sys.executable,str(ROOT/'scripts/build_integrated_sv_gene_tsv.py'),
            '--vcf',str(vcf),'--annotsv',str(self.sample),'--annotsv-audit',str(self.root/'audit.json'),
            '--output',str(self.root/'integrated.tsv')],check=True,capture_output=True)
        with open(self.root/'integrated.tsv') as fh:
            row=next(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual(row['SV_DB_CLINVAR_OVERLAP'],'UNKNOWN')
        self.assertEqual(row['SV_DB_CLINVAR_AVAILABILITY'],'NO_OVERLAP_REPORTED')
        self.assertEqual(row['SV_BENIGN_DB_STATUS'],'PROCESSED_RESOURCE_EMPTY')

    def test_cli_48_rows_and_stale_sample_guard(self):
        report=self.make_audit()
        subprocess.run([sys.executable,str(ROOT/'scripts/audit_annotsv_evidence.py'),
                        '--annotations-dir',str(self.root),'--annotsv',str(self.sample),
                        '--output-json',str(self.root/'audit.json'),'--output-tsv',str(self.root/'audit.tsv')],check=True,capture_output=True)
        with open(self.root/'audit.tsv') as fh:
            records=list(csv.DictReader(fh,delimiter='\t'))
        self.assertEqual(len(records),48)
        self.assertEqual(audit.load_audit(self.root/'audit.json',self.sample)['version'],audit.VERSION)
        self.sample.write_text(self.sample.read_text()+'changed\tchanged\n')
        with self.assertRaisesRegex(ValueError,'stale'):
            audit.load_audit(self.root/'audit.json',self.sample)


class PanelFilterTests(unittest.TestCase):
    def test_exact_panel_filter_and_nonpanel_complement(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'panel.txt').write_text('OPA1\n')
            tsv(p/'ann.tsv',['Gene_name'],[{'Gene_name':g} for g in ['OPA1','OPA10','OPA1;NEWGENE','opa1|OTHER','.']])
            subprocess.run([sys.executable,str(ROOT/'scripts/filter_annotsv_panel.py'),'--annotsv',str(p/'ann.tsv'),'--genes',str(p/'panel.txt'),'--output',str(p/'panel.tsv')],check=True,capture_output=True)
            with open(p/'panel.tsv') as fh:
                rows=list(csv.DictReader(fh,delimiter='\t'))
            self.assertEqual(len(rows),3)
            self.assertNotIn('OPA10',[r['Gene_name'] for r in rows])
            subprocess.run([sys.executable,str(ROOT/'scripts/extract_annotsv_genes.py'),'--annotsv',str(p/'ann.tsv'),'--output',str(p/'all.txt')],check=True,capture_output=True)
            subprocess.run([sys.executable,str(ROOT/'scripts/extract_nonpanel_genes.py'),'--genes',str(p/'all.txt'),'--panel',str(p/'panel.txt'),'--output',str(p/'nonpanel.txt')],check=True,capture_output=True)
            self.assertEqual(set((p/'nonpanel.txt').read_text().split()),{'OPA10','NEWGENE','OTHER'})

    def test_legacy_workflow_has_no_pick_command(self):
        legacy=(ROOT/'snakemake_pipelines/lrs/Snakefile_LRS').read_text()
        self.assertIn('include: "Snakefile_LRS_update"',legacy)
        self.assertNotIn('--pick',legacy)
        for name in ['lrs/Snakefile_LRS_update','srs_wgs_pipeline/Snakefile_SRS_WGS']:
            text=(ROOT/'snakemake_pipelines'/name).read_text()
            self.assertIn('--flag_pick',text)
            self.assertNotIn('--pick ',text)


if __name__ == '__main__':
    unittest.main()
