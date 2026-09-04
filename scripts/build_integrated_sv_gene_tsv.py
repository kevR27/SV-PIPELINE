#!/usr/bin/env python3
"""Build an integrated SV-gene TSV from canonical VCF + AnnotSV + needLR + ranking."""
from __future__ import annotations
import argparse,csv,gzip,re
from collections import OrderedDict,defaultdict
from pathlib import Path
MISSING='.'
def open_text(path):
    if path.endswith('.gz'): return gzip.open(path,'rt',encoding='utf-8',errors='replace')
    return open(path,'r',encoding='utf-8',errors='replace')
def parse_info(raw):
    out={}
    if not raw or raw==MISSING: return out
    for item in raw.split(';'):
        if not item: continue
        if '=' in item: k,v=item.split('=',1); out[k]=v
        else: out[item]='True'
    return out
def first(row,names):
    for n in names:
        v=row.get(n,MISSING)
        if v not in ('',MISSING,None): return str(v)
    low={k.lower():v for k,v in row.items()}
    for n in names:
        v=low.get(n.lower(),MISSING)
        if v not in ('',MISSING,None): return str(v)
    return MISSING
def split_genes(value):
    if value in ('',MISSING,None): return []
    return sorted({x.strip().upper() for x in re.split(r'[,;|]',str(value)) if x.strip()})
def read_tsv(path):
    if not path: return []
    with open_text(path) as fh: return list(csv.DictReader(fh,delimiter='\t'))
def read_panel(path):
    if not path: return set()
    with open_text(path) as fh: return {x.strip().upper() for x in fh if x.strip() and not x.startswith('#')}
def read_vcf(path):
    rows=[]; info_keys=OrderedDict()
    with open_text(path) as fh:
        for line in fh:
            if line.startswith('##INFO=<'):
                m=re.search(r'ID=([^,>]+)',line)
                if m: info_keys.setdefault(m.group(1),None)
                continue
            if line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<8: raise ValueError(f'Malformed VCF record: {line.rstrip()}')
            chrom,pos,sid,ref,alt,qual,filt,info_raw=f[:8]; info=parse_info(info_raw)
            for k in info: info_keys.setdefault(k,None)
            svtype=first(info,['SVTYPE'])
            if svtype==MISSING: svtype='BND' if '[' in alt or ']' in alt else (alt[1:-1] if alt.startswith('<') and alt.endswith('>') else MISSING)
            row=OrderedDict([('SV_ID',sid),('CHROM',chrom),('START',pos),('END',first(info,['END'])),('SVTYPE',svtype),('SVLEN',first(info,['SVLEN'])),('QUAL',qual),('FILTER',filt),('REF',ref),('ALT',alt),('INFO_RAW',info_raw)])
            for k,v in info.items(): row[f'INFO_{re.sub(r"[^A-Za-z0-9_.-]+","_",k)}']=v
            row['CALLERS']=first(info,['CALLERS','CALLER','SOURCE','SOURCES']); row['SUPPORT']=first(info,['SUPPORT','SUPP','RE','SU']); rows.append(row)
    return rows,info_keys
def index_by_id(rows):
    d=defaultdict(list)
    for r in rows:
        sid=first(r,['SV_ID','ID','AnnotSV_ID'])
        if sid!=MISSING: d[sid].append(r)
    return d
def genes_from_annot(row):
    genes=[]
    for c in ['Gene_name','Gene','GENE','Genes','gene','SYMBOL','GeneID','AnnotSV_Gene']: genes.extend(split_genes(first(row,[c])))
    return sorted(set(genes))
def load_ranking(rows):
    d={}
    for r in rows:
        g=first(r,['gene','Gene','GENE','SYMBOL'])
        if g!=MISSING: d[g.upper()]=r
    return d
def load_needlr(rows):
    d=defaultdict(list)
    for r in rows:
        g=first(r,['Gene','gene','Gene_name','GENE','SYMBOL'])
        for gene in split_genes(g): d[gene].append(r)
    return d
def main()->int:
    ap=argparse.ArgumentParser(description='Create integrated SV-gene analysis TSV.'); ap.add_argument('--vcf',required=True); ap.add_argument('--annotsv',required=True); ap.add_argument('--needlr'); ap.add_argument('--ranking'); ap.add_argument('--panel'); ap.add_argument('--caller-tsv',action='append',default=[]); ap.add_argument('--output',required=True); a=ap.parse_args()
    sv_rows,_=read_vcf(a.vcf); ann=read_tsv(a.annotsv); need=read_tsv(a.needlr); rank=load_ranking(read_tsv(a.ranking)); panel=read_panel(a.panel); ann_by=index_by_id(ann); need_by=load_needlr(need)
    caller_by=defaultdict(list)
    for path in a.caller_tsv:
        for r in read_tsv(path):
            sid=first(r,['SV_ID','ID'])
            if sid!=MISSING: caller_by[sid].append(r)
    final=[]
    for sv in sv_rows:
        sid=sv['SV_ID']; am=ann_by.get(sid,[]); genes=sorted({g for x in am for g in genes_from_annot(x)}) or [MISSING]; cm=caller_by.get(sid,[]); names=sorted({x.get('CALLER','') for x in cm if x.get('CALLER','')}); callers=';'.join(names) if names else sv['CALLERS']; sp=[]
        for x in cm:
            s=x.get('CALLER_SUPPORT',MISSING)
            if s not in ('',MISSING): sp.append(f"{x.get('CALLER','UNKNOWN')}:{s}")
        support=';'.join(sp) if sp else sv['SUPPORT']
        for gene in genes:
            ar=am[0] if am else {}; nr=need_by.get(gene,[{}])[0] if gene!=MISSING else {}; rr=rank.get(gene,{}) if gene!=MISSING else {}
            row=OrderedDict([('SV_ID',sv['SV_ID']),('CHROM',sv['CHROM']),('START',sv['START']),('END',sv['END']),('SVTYPE',sv['SVTYPE']),('SVLEN',sv['SVLEN']),('QUAL',sv['QUAL']),('FILTER',sv['FILTER']),('CALLERS',callers),('SUPPORT',support),('GENES',gene),('NEEDLR_AF',first(nr,['AF','MAX_AF','AF_MAX','SV_AF','AF_1KGP','1KGP_AF'])),('NEEDLR_HPO',first(nr,['HPO','HPO_terms','HPO_Terms','HPO phenotypes'])),('OMIM',first(ar,['OMIM','AnnotSV_OMIM','AnnotSV_OMIM_evidence'])),('GENCC',first(ar,['GENCC','GenCC','AnnotSV_GENCC','AnnotSV_GENCC_evidence'])),('ANNotsv_Gene',first(ar,['Gene_name','Gene','GENE','Genes','gene','SYMBOL','AnnotSV_Gene'])),('ANNotsv_Classification',first(ar,['AnnotSV ranking','AnnotSV_rank','AnnotSV_Classification'])),('PANEL_STATUS','PANEL_GENE' if gene in panel else ('UNRESOLVED' if gene==MISSING else 'NONPANEL_GENE')),('PHENOTYPE_SCORE',first(rr,['phenotype_score','PHENOTYPE_SCORE'])),('CANDIDATE_CLASS',first(rr,['classification','CANDIDATE_CLASS'])),('NEEDLR_STATUS','ANNOTATED' if nr else 'NO_MATCH')])
            if row['OMIM']==MISSING: row['OMIM']=first(nr,['OMIM','OMIM_phenotypes','OMIM phenotypes'])
            if row['GENCC']==MISSING: row['GENCC']=first(nr,['GENCC','GenCC','GenCC_phenotypes','GenCC phenotypes'])
            for k,v in sv.items():
                if k.startswith('INFO_'): row[k]=v
            final.append(row)
    fixed=['SV_ID','CHROM','START','END','SVTYPE','SVLEN','QUAL','FILTER','CALLERS','SUPPORT','GENES','NEEDLR_AF','NEEDLR_HPO','OMIM','GENCC','ANNotsv_Gene','ANNotsv_Classification','PANEL_STATUS','PHENOTYPE_SCORE','CANDIDATE_CLASS','NEEDLR_STATUS']; info_cols=sorted({k for r in final for k in r if k.startswith('INFO_')}); cols=fixed+info_cols
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); [w.writerow({c:r.get(c,MISSING) for c in cols}) for r in final]
    print(f'[OK] integrated_rows={len(final)} info_columns={len(info_cols)}'); print(f'[OK] output={out}'); return 0
if __name__=='__main__': raise SystemExit(main())
