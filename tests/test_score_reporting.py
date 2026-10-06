"""Regressions from the score/output audit; all fixtures are synthetic."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import assess_sv_alleles as allele
import rank_sv_gene_events as event_rank


def table(path, data, columns=None):
    columns = columns or list(dict.fromkeys(k for r in data for k in r))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, columns, delimiter='\t', restval='.')
        writer.writeheader(); writer.writerows(data)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


class ScoreReportingTests(unittest.TestCase):
    def ranking_resources(self, d):
        seeds = d / "hon_seeds.tsv"
        edges = d / "edges.tsv"
        table(
            seeds,
            [{
                "hpo_id": "HP:0000648",
                "hpo_label": "Optic atrophy",
                "hon_seed_role": "CORE_OCULAR_HON",
            }],
        )
        table(
            edges,
            [],
            ["subject", "object", "predicate", "category"],
        )
        return seeds, edges

    def command(self, script, *args):
        p = subprocess.run([sys.executable, str(ROOT / script), *map(str, args)],
                           text=True, capture_output=True, env={**os.environ, 'MPLBACKEND':'Agg'})
        self.assertEqual(p.returncode, 0, p.stderr)

    def test_large_field_and_null_inputs_preserve_original_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            original = [{'SV_ID':'a','GENES':'A','INFO_RAW':'X' * 350000,
                         'INTEGRATED_DISCOVERY_SCORE':'19','ACMG_CNV_SCORE':'0.99'},
                        {'SV_ID':'a','GENES':'B','INFO_RAW':'.',
                         'INTEGRATED_DISCOVERY_SCORE':'5','ACMG_CNV_SCORE':'0.99'},
                        {'SV_ID':'b','GENES':'.','INFO_RAW':'.',
                         'INTEGRATED_DISCOVERY_SCORE':'.','ACMG_CNV_SCORE':'.'}]
            table(d/'input.tsv', original)
            self.command('scripts/assess_sv_alleles.py', '--integrated', d/'input.tsv', '--sample','SYNTHETIC',
                         '--output',d/'out.tsv','--hypotheses-output',d/'hyp.tsv','--manifest',d/'manifest.json')
            output = read(d/'out.tsv')
            self.assertEqual(len(output), len(original))
            for a,b in zip(original,output):
                self.assertEqual(a,{k:b[k] for k in a})
                self.assertEqual(b['ALLELE_PHENOTYPE_STATUS'],'UNKNOWN_PATIENT_PHENOTYPES')

    def rank(self, d, annotations):
        (d/'genes').write_text('A\n');(d/'panel').write_text('A\n')
        table(d/'hpo.tsv',[],['gene_symbol','hpo_id','optic_neuropathy_anchor'])
        table(d/'ann.tsv',annotations)
        seeds, edges = self.ranking_resources(d)
        self.command(
            'scripts/rank_sv_gene_candidates.py',
            '--annotsv',d/'ann.tsv','--genes',d/'genes',
            '--panel',d/'panel','--phenotypes',d/'hpo.tsv',
            '--hpo-seeds',seeds,'--edges',edges,
            '--output',d/'rank.tsv'
        )
        return read(d/'rank.tsv')[0]

    def test_missing_omim_is_not_disease_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            for value in ['.','NA',' N/A ','None','NaN','null']:
                with self.subTest(value=value):
                    r=self.rank(Path(tmp),[{'SV_ID':'v','Gene_name':'A','Annotation_mode':'split','OMIM_phenotype':value}])
                    self.assertEqual(r['gene_disease_evidence_score'],'0.0')

    def test_clinvar_alone_does_not_define_human_disease_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d/'genes').write_text('A\n')
            (d/'panel').write_text('')
            table(d/'hpo.tsv',[],['gene_symbol','hpo_id','optic_neuropathy_anchor'])
            table(d/'ann.tsv',[{
                'SV_ID':'v','Gene_name':'A','Annotation_mode':'split',
                'ClinVar':'pathogenic_variant_overlap'
            }])
            seeds, edges = self.ranking_resources(d)
            self.command(
                'scripts/rank_sv_gene_candidates.py',
                '--annotsv',d/'ann.tsv','--genes',d/'genes',
                '--panel',d/'panel','--phenotypes',d/'hpo.tsv',
                '--hpo-seeds',seeds,'--edges',edges,
                '--output',d/'rank.tsv'
            )
            r = read(d/'rank.tsv')[0]
            self.assertEqual(r['gene_disease_evidence_score'],'0.0')
            self.assertEqual(r['candidate_group'],'OTHER_NONPANEL_CANDIDATE')

    def test_animal_only_record_retained_without_human_evidence_points(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=self.rank(Path(tmp),[{'SV_ID':'v','Gene_name':'A','Annotation_mode':'split','GenCC_classification':'Animal Model Only'}])
            self.assertEqual(r['gene_disease_evidence_score'],'0.0')
            self.assertEqual(r['GENCC'],'ANIMAL MODEL ONLY')

    def test_full_sv_class_retained_without_transferring_gene_disease(self):
        with tempfile.TemporaryDirectory() as tmp:
            r=self.rank(Path(tmp),[
                {'SV_ID':'v','Gene_name':'A;B','Annotation_mode':'full','ACMG_class':'5',
                 'AnnotSV_ranking_score':'0.99','OMIM_phenotype':'UNRELATED_GENE_DISEASE'},
                {'SV_ID':'v','Gene_name':'A','Annotation_mode':'split','ACMG_class':'.','AnnotSV_ranking_score':'.'}])
            self.assertEqual(r['AnnotSV_ACMG_classes'],'5')
            self.assertEqual(r['AnnotSV_ranking_scores'],'0.99')
            self.assertEqual(r['gene_disease_evidence_score'],'0.0')
            self.assertIn('NOT_GENE_PATHOGENICITY',r['AnnotSV_classification_scope'])

    def test_bnd_missing_length_is_informational(self):
        row={'SVTYPE':'BND','CALLER_COUNT':'2','FILTER':'PASS','CALLER_EVIDENCE_MATCH':'IDLIST',
             'CALLER_EVIDENCE_JSON':json.dumps([{'CALLER_SUPPORT':'20','EVIDENCE_FLAGS':'NO_SVLEN'}])}
        self.assertEqual(allele.technical(row,5)[1],2)
        self.assertIn('NO_SVLEN',allele.technical(row,5)[2])
        self.assertEqual(allele.technical({**row,'SVTYPE':'DEL'},5)[1],0)
        row['CALLER_EVIDENCE_JSON']=json.dumps([{'CALLER_SUPPORT':'20','EVIDENCE_FLAGS':'NO_SVLEN;LOW_GQ'}])
        self.assertEqual(allele.technical(row,5)[1],0)

    def test_absent_complementary_sources_are_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);table(d/'input.tsv',[{'SV_ID':'v','GENES':'A','CHROM':'chr1','START':'100','END':'200','SVTYPE':'DEL'}])
            self.command('plots/intersect_complementary_sv_evidence.py','--integrated',d/'input.tsv','--output',d/'out.tsv')
            r=read(d/'out.tsv')[0]
            for name in ['STRAGLR_MATCH','TLDR_MATCH','LONGPHASE_MATCH','LONGPHASE_PHASED']:
                self.assertEqual(r[name],'NOT_AVAILABLE')
            self.command('plots/build_gene_evidence_summary.py','--integrated',d/'out.tsv','--output',d/'summary.tsv')
            self.assertEqual(read(d/'summary.tsv')[0]['longphase_match_unavailable_SV_count'],'1')

    def test_event_technical_review_flags_qc_rescue(self):
        self.assertEqual(
            event_rank.technical_review_status({
                "CALLER_EVIDENCE_FLAGS": "Sniffles2:RESCUED_COV_VAR",
                "CALLER_EVIDENCE_MATCH": "IDLIST",
            }),
            "REVIEW_REQUIRED_CALLER_QC_FLAG",
        )
        self.assertEqual(
            event_rank.technical_review_status({
                "CALLER_EVIDENCE_FLAGS": ".",
                "CALLER_EVIDENCE_MATCH": "IDLIST",
            }),
            "NO_REVIEW_FLAG_FROM_CALLER_EVIDENCE",
        )
        self.assertEqual(
            event_rank.technical_review_status({
                "CALLER_EVIDENCE_FLAGS": ".",
                "CALLER_EVIDENCE_MATCH": "AMBIGUOUS_COORDINATE",
            }),
            "REVIEW_REQUIRED_AMBIGUOUS_CALLER_LINK",
        )

    def test_missing_population_evidence_is_neutral(self):
        priority = event_rank.POPULATION_PRIORITY
        self.assertGreater(priority["PROVISIONAL_LOW_AF_LE_0.01"], priority["NO_MATCH"])
        self.assertEqual(priority["NO_MATCH"], priority["UNKNOWN_OR_MISSING"])
        self.assertEqual(priority["NOT_EVALUABLE_BND"], priority["UNKNOWN_OR_MISSING"])
        self.assertEqual(priority["NOT_EVALUABLE_GE10MB"], priority["UNKNOWN_OR_MISSING"])
        self.assertGreater(priority["UNKNOWN_OR_MISSING"], priority["PROVISIONAL_HIGH_AF_GT_0.01"])

    def test_matrix_keeps_distinct_ids_and_follows_gene_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            base={'SV_ID':'a','GENES':'A','CHROM':'chr1','START':'100','END':'200','SVTYPE':'DEL',
                  'PHENOTYPE_SCORE':'1','INTEGRATED_DISCOVERY_SCORE':'9','PANEL_STATUS':'NONPANEL_GENE',
                  'LONGPHASE_MATCH':'NOT_AVAILABLE','CALLERS':'Manta;Delly'}
            table(d/'in.tsv',[base,{**base,'SV_ID':'b','END':'300'},
                              {**base,'SV_ID':'c','GENES':'C','PANEL_STATUS':'PANEL_GENE','INTEGRATED_DISCOVERY_SCORE':'1'}])
            self.command('plots/plot_candidate_evidence_matrix.py','--input',d/'in.tsv','--out-prefix',d/'matrix')
            result=read(d/'matrix_matrix.tsv')
            self.assertEqual([r['SV_ID'] for r in result],['a','b','c'])
            self.assertEqual(result[0]['PLOT_ORDER_BASIS'],'INTEGRATED_DISCOVERY_SCORE')
            self.assertEqual(result[0]['LongPhase match'],'')
            self.assertEqual(result[0]['Manta'],'1')

    def test_empty_context_and_gene_summary_are_valid_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);table(d/'in.tsv',[],['SV_ID','GENES','CHROM','START','END','SVTYPE'])
            self.command('plots/integrate_candidate_context.py','--integrated',d/'in.tsv','--output',d/'context.tsv')
            self.command('plots/build_gene_evidence_summary.py','--integrated',d/'context.tsv','--output',d/'summary.tsv')
            self.assertEqual(read(d/'summary.tsv'),[])
            self.assertIn('gene',(d/'summary.tsv').read_text())
