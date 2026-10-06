"""Regression fixtures for evidence misattribution and conservative assessment."""
import csv
import importlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'scripts'), str(ROOT / 'plots')]
import build_integrated_sv_gene_tsv as integrate
import sv_evidence_common as common
import assess_sv_alleles as allele
import parse_sv_caller_vcf as caller


def table(path, data, fields=None):
    fields = fields or list(data[0])
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(data)


class CoreTests(unittest.TestCase):
    def test_database_missingness(self):
        for scope, p, b, expected in [('LOSS','.','.','UNKNOWN'), ('LOSS','dbVar','.','NOT_REPORTED'), ('LOSS','ClinVar','.','YES'), ('NO_DEDICATED_ANNOTSV_SV_DATABASE_FIELDS','.','.','NOT_APPLICABLE')]:
            self.assertEqual(integrate.source_overlap_flag(p,b,'CLINVAR',scope=scope), expected)

    def test_population_shared_parser(self):
        self.assertEqual(common.population_frequency({'Allele_Freq_ALL_Control':'0.002','Allele_Freq_ALL':'0.5'}), ('0.002','Allele_Freq_ALL_Control','REPORTED'))
        self.assertEqual(common.population_frequency({'Pop_Freq_ALL':'nan'})[2], 'UNKNOWN')
        self.assertEqual(common.population_frequency({'AF':'inf'})[2], 'INVALID_FREQUENCY')
        self.assertEqual(common.population_frequency({'AF':'2'})[2], 'INVALID_FREQUENCY')
        self.assertEqual(common.population_frequency({'AF':'0'})[0], '0.0')

    def test_annotsv_gene_scope(self):
        full={'Annotation_mode':'full','Gene_name':'A;B','OMIM_phenotype':'WRONG','P_loss_source':'dbVar','ACMG_class':'3'}
        split={'Annotation_mode':'split','Gene_name':'A','OMIM_phenotype':'RIGHT','HI':'3','Tx':'t1'}
        result=integrate.annotsv_row_for_gene([full,split], 'A')
        self.assertEqual(result['OMIM_phenotype'],'RIGHT')
        self.assertEqual(result['P_loss_source'],'dbVar')
        self.assertEqual(result['HI'],'3')
        self.assertEqual(result['ACMG_class'],'3')
        self.assertEqual(integrate.annotsv_row_for_gene([full], 'B')['OMIM_phenotype'], '.')
        self.assertEqual(integrate.genes_from_annot({'Gene_name':'A','GeneID':'123'}), ['A'])

    def test_transcript_disagreement_not_arbitrary(self):
        result=integrate.annotsv_row_for_gene([{'Gene_name':'A','Annotation_mode':'split','Tx':'t1'}, {'Gene_name':'A','Annotation_mode':'split','Tx':'t2'}], 'A')
        self.assertEqual(result['Tx'],'.')
        self.assertEqual(len(json.loads(result['_GENE_ROWS'])),2)

    def test_breakend_partner_and_orientation(self):
        a={'CHROM':'chr1','START':'1000','SVTYPE':'BND','ALT':'N[chr2:2000['}
        self.assertIsNotNone(common.same_breakend(a,dict(a)))
        for alt in ['N[chr9:2000[','N[chr2:9000[','N]chr2:2000]']:
            b={**a,'ALT':alt}
            self.assertIsNone(integrate.evidence_match_score(a,b))
        self.assertIsNone(common.same_breakend(a,{'CHROM':'chr1','START':'1000','CHR2':'chr2','POS2':'2000'}))

    def test_needlr_ambiguity(self):
        sv={'SV_ID':'v','CHROM':'chr1','START':'1000','END':'1500','SVLEN':'-500','SVTYPE':'DEL'}
        idx=integrate.build_needlr_index([dict(sv),{**sv,'SV_ID':'v2','START':'1001'}])
        self.assertEqual(integrate.match_needlr(sv,idx)[1], 'AMBIGUOUS_MATCH')
        self.assertEqual(integrate.match_needlr({**sv,'SVTYPE':'BND'},idx)[1], 'NOT_EVALUABLE_BND')
        self.assertEqual(integrate.match_needlr({**sv,'SVLEN':'10000000'},idx)[1], 'NOT_EVALUABLE_GE_10MB')

    def test_caller_ambiguity_not_best_arbitrary(self):
        sv={'CHROM':'chr1','START':'1000','END':'1500','SVLEN':'-500','SVTYPE':'DEL'}
        candidates=[{**sv,'CALLER':'Sniffles2','SV_ID':'x'}, {**sv,'CALLER':'Sniffles2','SV_ID':'y'}]
        matches,method=integrate.caller_matches_for_sv(sv,{},integrate.build_spatial_index(candidates))
        self.assertEqual(matches,[])
        self.assertEqual(method,'AMBIGUOUS_COORDINATE')

    def test_jasmine_prefixed_source_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'caller.tsv'
            table(p,[{'SV_ID':'Sniffles2.DEL.1','CALLER':'Sniffles2','CHROM':'chr1','START':'100','END':'200','SVTYPE':'DEL'}])
            by_id,index=integrate.load_caller_evidence([str(p)],jasmine_prefixes=True)
            found,method=integrate.caller_matches_for_sv({'INFO_IDLIST':'0_Sniffles2.DEL.1'},by_id,index)
            self.assertEqual(method,'IDLIST')
            self.assertEqual(found[0]['SV_ID'],'Sniffles2.DEL.1')

    def test_partial_depth_is_not_total_depth(self):
        self.assertEqual(caller.parse_format('GT:GQ:DV','0/1:40:8')['CALLER_DP'],'.')
        self.assertEqual(caller.parse_format('GT:GQ:DR:DV','0/1:40:8:7')['CALLER_DP'],'15')


class AssessmentTests(unittest.TestCase):
    def setUp(self):
        self.args=SimpleNamespace(sample='P',rare_af=.01,min_support=5,min_gq=20,min_dp=5)
        self.row={'SV_ID':'v','GENES':'A','CHROM':'chr1','START':'100','END':'200','SVTYPE':'DEL','CALLER_COUNT':'2','FILTER':'PASS','CALLER_EVIDENCE_MATCH':'IDLIST','CALLER_EVIDENCE_JSON':json.dumps([{'CALLER':'Sniffles2','CALLER_SUPPORT':'8','CALLER_GT':'0/1','CALLER_GQ':'30','CALLER_DP':'20'}])}
        self.model={'disease_id':'D:1','moi':'AD','mechanism':'LOF','validity':'STRONG','hpo_terms':'HP:0000648','source':'curation'}

    def test_header_template_missing_optional_values(self):
        self.assertEqual(allele.disruption(self.row,{'effect':'.','source':'review'})[1], 'UNKNOWN')
        self.assertEqual(allele.population({}, {'population_af':'.'},.01)[0], 'UNKNOWN')

    def test_unknown_never_gets_absence_bonus(self):
        r=allele.assess({'SV_ID':'v','GENES':'A'}, {},{}, {},{}, {},self.args)
        self.assertEqual(r['ALLELE_RESEARCH_SCORE'],0)
        self.assertEqual(len(r['ALLELE_UNKNOWN_DOMAINS'].split(';')),6)

    def test_patient_phenotype_not_generic_anchor_score(self):
        r=allele.assess({**self.row,'PHENOTYPE_SCORE':'13'},self.model,{}, {},{}, {'A':{'HP:0000648'}},self.args)
        self.assertEqual(r['ALLELE_PHENOTYPE_STATUS'],'UNKNOWN_PATIENT_PHENOTYPES')

    def test_positive_and_negative_phenotypes(self):
        self.assertEqual(allele.phenotype('A',self.model,{'hpo_present':'HP:0000648'}, {})[1],2)
        self.assertEqual(allele.phenotype('A',self.model,{'hpo_present':'HP:0000001','hpo_absent':'HP:0000648'}, {})[0],'PATIENT_TERM_CONFLICT_REVIEW')

    def test_needlr_only_is_provisional(self):
        r=allele.population({'NEEDLR_STATUS':'MATCHED','NEEDLR_AF':'0'}, {},.01)
        self.assertEqual(r[0],'PROVISIONAL_LOW_FREQUENCY_MATCH')
        self.assertEqual(r[1],0)

    def test_exact_population_and_high_af(self):
        ev={'population_af':'.001','population_match':'EXACT_ALLELE','population_source':'resource'}
        self.assertEqual(allele.population({},ev,.01)[1],2)
        self.assertEqual(allele.population({},{**ev,'population_af':'.2'},.01)[0],'ABOVE_RESEARCH_AF_THRESHOLD')
        self.assertEqual(allele.population({},{**ev,'population_af':'0'},.01)[1],1)

    def test_gene_overlap_does_not_establish_lof(self):
        self.assertEqual(allele.disruption(self.row,{})[1],'UNKNOWN')
        row={**self.row,'ANNOTSV_GENE_ROWS_JSON':json.dumps([{'Overlapped_CDS_percent':'25'}])}
        self.assertEqual(allele.disruption(row,{})[1],'POSSIBLE_LOF')

    def test_disease_mismatch_and_refutation(self):
        self.assertEqual(allele.disease_mechanism(self.row,{**self.model,'mechanism':'GOF'},'LOF')[0],'MECHANISM_MISMATCH_REVIEW')
        self.assertEqual(allele.disease_mechanism(self.row,{**self.model,'validity':'REFUTED'},'LOF')[1],0)

    def test_coordinate_link_cannot_transfer_genotype(self):
        gt,_=allele.patient_genotype({**self.row,'CALLER_EVIDENCE_MATCH':'COORDINATE_FALLBACK'}, {},{},'P',20,5)
        self.assertEqual(gt,'UNKNOWN')

    def test_genotype_conflict(self):
        records=json.loads(self.row['CALLER_EVIDENCE_JSON'])
        records.append({**records[0],'CALLER':'Delly','CALLER_GT':'1/1'})
        gt,_=allele.patient_genotype({**self.row,'CALLER_EVIDENCE_JSON':json.dumps(records)}, {},{},'P',20,5)
        self.assertEqual(gt,'CONFLICT')

    def test_ar_requires_second_allele(self):
        r=allele.inheritance(self.row,{**self.model,'moi':'AR'}, {},'HET',{},20,5)
        self.assertEqual(r[0],'SECOND_ALLELE_OR_TRANS_PHASE_REQUIRED')
        self.assertEqual(r[1],0)
        self.assertEqual(allele.inheritance(self.row,{**self.model,'moi':'AR'}, {},'HOM_ALT',{},20,5)[1],2)

    def test_no_de_novo_from_missing_parent_calls(self):
        patient={'mother_id':'M','father_id':'F'}
        self.assertEqual(allele.inheritance(self.row,self.model,patient,'HET',{},20,5)[0],'DOMINANT_GENOTYPE_COMPATIBLE')
        family={(x,'v'):{'gt':'0/0','gq':'30','dp':'20','source':'joint_genotyping'} for x in ['M','F']}
        self.assertEqual(allele.inheritance(self.row,self.model,patient,'HET',family,20,5)[0],'CANDIDATE_DE_NOVO')

    def test_six_domains_no_clinical_label(self):
        evidence={'effect':'LOF','effect_status':'VALIDATED','source':'assay','population_af':'.001','population_match':'EXACT_ALLELE','population_source':'resource'}
        r=allele.assess(self.row,self.model,{'hpo_present':'HP:0000648'},evidence,{}, {},self.args)
        self.assertEqual(r['ALLELE_RESEARCH_SCORE'],11)
        self.assertEqual(r['ALLELE_INTERPRETATION'],'RESEARCH_PRIORITIZATION_ONLY_NOT_ACMG')

    def test_family_trans_candidate_requires_distinct_compatible_allele(self):
        model={**self.model,'moi':'AR'}
        other={**self.row,'SV_ID':'w','START':'300','END':'400','ANNOTSV_GENE_ROWS_JSON':'[{"Overlapped_CDS_percent":"50"}]'}
        family={}
        for sid,mg,fg in [('v','0/1','0/0'),('w','0/0','0/1')]:
            for sample,gt in [('M',mg),('F',fg),('P','0/1')]:
                family[(sample,sid)]={'gt':gt,'gq':'30','dp':'20','source':'joint'}
        self.args.gene_rows={'A':[self.row,other]}
        found=allele.trans_candidates(self.row,model,{'mother_id':'M','father_id':'F'},family,self.args)
        self.assertEqual(found,['w'])
        other['START']='150'
        self.assertEqual(allele.trans_candidates(self.row,model,{'mother_id':'M','father_id':'F'},family,self.args),[])


class CommandTests(unittest.TestCase):
    def test_cli_preserves_multi_gene_and_intergenic_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            data=[{'SV_ID':'x','GENES':'A'},{'SV_ID':'x','GENES':'B'},{'SV_ID':'y','GENES':'.'}]
            table(d/'input.tsv',data)
            models=[{'gene':'A','disease_id':'D:1','moi':'AD','mechanism':'LOF','validity':'STRONG','source':'curation'}, {'gene':'A','disease_id':'D:2','moi':'AR','mechanism':'GOF','validity':'LIMITED','source':'curation'}]
            table(d/'models.tsv',models)
            cmd=[sys.executable,str(ROOT/'scripts/assess_sv_alleles.py'),'--integrated',str(d/'input.tsv'),'--sample','P','--disease-models',str(d/'models.tsv'),'--output',str(d/'out.tsv'),'--hypotheses-output',str(d/'hyp.tsv'),'--manifest',str(d/'manifest.json')]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            out=allele.rows(d/'out.tsv');hyp=allele.rows(d/'hyp.tsv')
            self.assertEqual([(r['SV_ID'],r['GENES']) for r in out],[('x','A'),('x','B'),('y','.')])
            self.assertEqual(len(hyp),4)
            self.assertTrue(json.loads((d/'manifest.json').read_text())['input_sha256']['integrated'])
            # Empty inputs still have all output columns.
            table(d/'input.tsv',[],['SV_ID','GENES'])
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            self.assertEqual(allele.rows(d/'out.tsv'),[])
            self.assertIn('ALLELE_RESEARCH_SCORE',(d/'out.tsv').read_text())

    def test_hpo_duplicates_do_not_inflate_ranking(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'genes').write_text('A\n');(d/'panel').write_text('A\n')
            table(d/'ann.tsv',[{'Gene_name':'A','Annotation_mode':'split','AnnotSV_ID':'v'}])
            pheno={'gene_symbol':'A','hpo_id':'HP:0000648','optic_neuropathy_anchor':'1','source':'human'}
            table(
                d/'hon_seeds.tsv',
                [{
                    'hpo_id':'HP:0000648',
                    'hpo_label':'Optic atrophy',
                    'hon_seed_role':'CORE_OCULAR_HON',
                }],
            )
            table(
                d/'edges.tsv',
                [],
                ['subject','object','predicate','category'],
            )
            scores=[]
            for n in [1,2]:
                table(d/'hpo.tsv',[pheno]*n)
                subprocess.run(
                    [
                        sys.executable,
                        str(ROOT/'scripts/rank_sv_gene_candidates.py'),
                        '--annotsv',str(d/'ann.tsv'),
                        '--genes',str(d/'genes'),
                        '--panel',str(d/'panel'),
                        '--phenotypes',str(d/'hpo.tsv'),
                        '--hpo-seeds',str(d/'hon_seeds.tsv'),
                        '--edges',str(d/'edges.tsv'),
                        '--output',str(d/'out.tsv'),
                    ],
                    check=True,
                    capture_output=True,
                )
                scores.append(allele.rows(d/'out.tsv')[0]['phenotype_score'])
            self.assertEqual(scores[0],scores[1])

    def test_integrator_end_to_end_retains_unannotated_and_gene_specific(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            (d/'in.vcf').write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\tv\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200;SVLEN=-100\nchr2\t100\tb\tN\tN[chr9:300[\t.\tPASS\tSVTYPE=BND\n')
            table(d/'ann.tsv',[{'SV_ID':'v','Gene_name':'A;B','Annotation_mode':'full','HI':'.','P_loss_source':'dbVar'}, {'SV_ID':'v','Gene_name':'A','Annotation_mode':'split','HI':'3','P_loss_source':'.'}])
            subprocess.run([sys.executable,str(ROOT/'scripts/build_integrated_sv_gene_tsv.py'),'--vcf',str(d/'in.vcf'),'--annotsv',str(d/'ann.tsv'),'--output',str(d/'out.tsv')],check=True,capture_output=True)
            out=allele.rows(d/'out.tsv')
            self.assertEqual(len(out),3)
            self.assertEqual(out[0]['CLINGEN_HI'],'3')
            self.assertEqual(out[1]['CLINGEN_HI'],'.')
            self.assertEqual(out[2]['CHR2'],'chr9')
            self.assertEqual(out[2]['POS2'],'300')
            self.assertEqual(out[2]['SV_DB_CLINVAR_OVERLAP'],'NOT_APPLICABLE')
            self.assertTrue(json.loads(out[0]['ANNOTSV_GENE_ROWS_JSON']))

    def test_needlr_reuses_complete_native_results_only(self):
        snakefile = (ROOT / "snakemake_pipelines/lrs/Snakefile_LRS_update").read_text()
        self.assertIn("rule needlr_annotation:", snakefile)
        self.assertIn("query_hash=", snakefile)
        self.assertIn("sha256sum", snakefile)
        self.assertIn(".needlr_query.sha256", snakefile)
        self.assertIn("Native result matches current query; reusing it.", snakefile)
        self.assertIn("Query changed or no validated native result; running annotation.", snakefile)
        self.assertIn('if [ -s "$results_tsv" ] && [ -s "$results_vcf" ]', snakefile)
        self.assertIn("needLR annotate", snakefile)

    def test_methylation_units(self):
        context=importlib.import_module('integrate_candidate_context')
        plot=importlib.import_module('plot_methylation')
        line='chr1\t1\t2\th\t0\t+\t1\t2\t0\t10\t0.5'
        self.assertEqual(context.parse_bedmethyl_record(line,context.methylation_scale('percent'),5)['percent'],.5)
        self.assertEqual(plot.detect_scale('percent'),1)
        self.assertEqual(context.parse_bedmethyl_record(line,context.methylation_scale('fraction'),5)['percent'],50)

    def test_panel_and_unknown_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            table(d/'in.tsv',[{'SV_ID':'v','GENES':'A','PANEL_STATUS':'NONPANEL_GENE','SV_DB_CLINVAR_OVERLAP':'UNKNOWN'}, {'SV_ID':'w','GENES':'B','PANEL_STATUS':'PANEL_GENE','SV_DB_CLINVAR_OVERLAP':'NOT_REPORTED'}])
            subprocess.run([sys.executable,str(ROOT/'plots/build_gene_evidence_summary.py'),'--integrated',str(d/'in.tsv'),'--output',str(d/'out.tsv')],check=True,capture_output=True)
            out={r['gene']:r for r in allele.rows(d/'out.tsv')}
            self.assertEqual(out['A']['panel_gene'],'NO')
            self.assertEqual(out['B']['panel_gene'],'YES')
            self.assertEqual(out['A']['clinvar_SV_unknown_count'],'1')
            self.assertEqual(out['B']['clinvar_SV_not_reported_count'],'1')


if __name__ == '__main__':
    unittest.main()
