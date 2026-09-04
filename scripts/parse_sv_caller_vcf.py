#!/usr/bin/env python3
"""Create a compact caller-aware SV evidence table from a VCF/VCF.GZ."""
from __future__ import annotations
import argparse,csv,gzip
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
def first(info,*keys):
    for k in keys:
        v=info.get(k,MISSING)
        if v not in ('',MISSING): return v
    return MISSING
def infer_svtype(info,alt):
    v=first(info,'SVTYPE')
    if v!=MISSING: return v
    if '[' in alt or ']' in alt: return 'BND'
    if alt.startswith('<') and alt.endswith('>'): return alt[1:-1]
    return MISSING
def evidence(caller,info):
    c=caller.lower()
    if c=='sniffles2': return {'CALLER_SUPPORT':first(info,'SUPPORT','RE','SUPP'),'CALLER_RNAMES':first(info,'RNAMES'),'CALLER_STRANDS':first(info,'STRANDS'),'CALLER_IMPRECISE':first(info,'IMPRECISE'),'CALLER_MOSAIC':first(info,'MOSAIC')}
    if c=='cutesv': return {'CALLER_SUPPORT':first(info,'RE','SUPPORT','SUPP'),'CALLER_RNAMES':first(info,'RNAMES'),'CALLER_STRANDS':first(info,'STRANDS')}
    if c=='delly': return {'CALLER_SUPPORT':first(info,'SU','SUPPORT'),'CALLER_PE':first(info,'PE'),'CALLER_SR':first(info,'SR'),'CALLER_PRECISE':first(info,'PRECISE')}
    if c=='manta': return {'CALLER_SUPPORT':first(info,'SU','PR','SR','SUPPORT'),'CALLER_PR':first(info,'PR'),'CALLER_SR':first(info,'SR')}
    if c in {'jasmine','survivor'}: return {'CALLER_SUPPORT':first(info,'SUPP','SUPPORT'),'CALLER_SUPPORT_VECTOR':first(info,'SUPP_VEC')}
    return {'CALLER_SUPPORT':first(info,'SUPPORT','SUPP','RE','SU')}
def main()->int:
    ap=argparse.ArgumentParser(description='Parse caller-specific SV evidence.'); ap.add_argument('--vcf',required=True); ap.add_argument('--caller',required=True); ap.add_argument('--output',required=True); a=ap.parse_args()
    rows=[]
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<8: raise ValueError(f'Malformed VCF record: {line.rstrip()}')
            chrom,pos,sid,ref,alt,qual,filt,info_raw=f[:8]; info=parse_info(info_raw)
            row={'CALLER':a.caller,'SV_ID':sid,'CHROM':chrom,'START':pos,'END':first(info,'END'),'SVTYPE':infer_svtype(info,alt),'SVLEN':first(info,'SVLEN'),'QUAL':qual,'FILTER':filt,'REF':ref,'ALT':alt,'INFO_RAW':info_raw}
            row.update(evidence(a.caller,info)); rows.append(row)
    cols=['CALLER','SV_ID','CHROM','START','END','SVTYPE','SVLEN','QUAL','FILTER','REF','ALT','INFO_RAW','CALLER_SUPPORT','CALLER_RNAMES','CALLER_STRANDS','CALLER_IMPRECISE','CALLER_MOSAIC','CALLER_PE','CALLER_SR','CALLER_PR','CALLER_PRECISE','CALLER_SUPPORT_VECTOR']
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter='\t',extrasaction='ignore',lineterminator='\n'); w.writeheader(); w.writerows(rows)
    print(f'[OK] caller={a.caller} records={len(rows)} output={out}'); return 0
if __name__=='__main__': raise SystemExit(main())
