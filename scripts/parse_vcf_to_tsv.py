#!/usr/bin/env python3
"""Expand a VCF/VCF.GZ into a TSV without losing INFO or FORMAT data."""
from __future__ import annotations
import argparse,gzip,re,sys
from collections import OrderedDict
from pathlib import Path
from typing import TextIO
MISSING='.'
def open_text(path:str)->TextIO:
    if path=='-': return sys.stdin
    if path.endswith('.gz'): return gzip.open(path,'rt',encoding='utf-8',errors='replace')
    return open(path,'r',encoding='utf-8',errors='replace')
def clean_key(key:str)->str: return re.sub(r'[^A-Za-z0-9_.-]+','_',key).strip('_') or 'UNKNOWN'
def parse_info(raw:str)->OrderedDict[str,str]:
    out=OrderedDict()
    if not raw or raw==MISSING: return out
    for item in raw.split(';'):
        if not item: continue
        if '=' in item:
            k,v=item.split('=',1); out[k]=v
        else: out[item]='True'
    return out
def parse_format(raw:str,samples:list[str],values:list[str])->OrderedDict[str,str]:
    out=OrderedDict()
    if not raw or raw==MISSING: return out
    keys=raw.split(':')
    for i,sample in enumerate(samples):
        if i>=len(values): break
        vals=values[i].split(':'); sk=clean_key(sample)
        for j,key in enumerate(keys): out[f'FORMAT_{sk}_{clean_key(key)}']=vals[j] if j<len(vals) else MISSING
    return out
def first(info:dict[str,str],*keys:str)->str:
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
def infer_end(pos:int,ref:str,info)->str:
    v=first(info,'END')
    if v!=MISSING: return v
    return str(pos+max(len(ref),1)-1) if ref not in ('',MISSING) else MISSING
def main()->int:
    ap=argparse.ArgumentParser(description='Expand every VCF INFO/FORMAT key into TSV columns.')
    ap.add_argument('--vcf',required=True); ap.add_argument('--output',required=True); ap.add_argument('--caller',default='')
    a=ap.parse_args(); info_keys=OrderedDict(); format_keys=OrderedDict(); samples=[]; rows=[]
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith('##INFO=<'):
                m=re.search(r'ID=([^,>]+)',line)
                if m: info_keys.setdefault(m.group(1),None)
            elif line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:]; break
    with open_text(a.vcf) as fh:
        for line in fh:
            if line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<8: raise ValueError(f'Malformed VCF record: {line.rstrip()}')
            chrom,pos,sid,ref,alt,qual,filt,info_raw=f[:8]
            format_raw=f[8] if len(f)>8 else MISSING; sample_values=f[9:] if len(f)>9 else []
            try: pos_int=int(pos)
            except ValueError as e: raise ValueError(f'Invalid POS {pos!r}') from e
            info=parse_info(info_raw)
            for k in info: info_keys.setdefault(k,None)
            fmt=parse_format(format_raw,samples,sample_values)
            for k in fmt: format_keys.setdefault(k,None)
            row=OrderedDict()
            if a.caller: row['CALLER']=a.caller
            row.update({'SV_ID':sid,'CHROM':chrom,'START':pos,'END':infer_end(pos_int,ref,info),'SVTYPE':infer_svtype(info,alt),'SVLEN':first(info,'SVLEN'),'QUAL':qual,'FILTER':filt,'REF':ref,'ALT':alt,'INFO_RAW':info_raw,'FORMAT_RAW':format_raw})
            for k in info_keys: row[f'INFO_{clean_key(k)}']=info.get(k,MISSING)
            for k in format_keys: row[k]=fmt.get(k,MISSING)
            rows.append(row)
    cols=(['CALLER'] if a.caller else [])+['SV_ID','CHROM','START','END','SVTYPE','SVLEN','QUAL','FILTER','REF','ALT','INFO_RAW','FORMAT_RAW']+[f'INFO_{clean_key(k)}' for k in info_keys]+list(format_keys)
    out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',encoding='utf-8',newline='') as fh:
        fh.write('\t'.join(cols)+'\n')
        for row in rows: fh.write('\t'.join(row.get(c,MISSING) for c in cols)+'\n')
    print(f'[OK] records={len(rows)} info_fields={len(info_keys)} format_fields={len(format_keys)}')
    print(f'[OK] output={out}'); return 0
if __name__=='__main__': raise SystemExit(main())
