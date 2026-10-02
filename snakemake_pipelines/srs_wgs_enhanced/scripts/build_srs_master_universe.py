#!/usr/bin/env python3
"""Build an SRS master SV universe by adding independent CNVpytor depth calls.

The breakpoint universe (Manta+DELLY merged by SURVIVOR) remains the primary
breakpoint representation. CNVpytor is an orthogonal read-depth caller. Same-type
DEL/DUP calls with sufficient reciprocal overlap are attached as RD evidence.
Quality-passing CNVpytor calls without a breakpoint match are retained as
CNVPYTOR_ONLY symbolic DEL/DUP records so depth-only CNVs are not silently lost.

This script does NOT reinterpret caller agreement as pathogenicity or clinical
confidence.
"""
from __future__ import annotations
import argparse, gzip, math, re
from pathlib import Path

MISSING="."

def opentext(path):
    return gzip.open(path,"rt",encoding="utf-8",errors="replace") if str(path).endswith(".gz") else open(path,encoding="utf-8",errors="replace")

def info_dict(raw):
    d={}
    for item in str(raw).split(";"):
        if not item: continue
        if "=" in item:
            k,v=item.split("=",1); d[k]=v
        else: d[item]="True"
    return d

def num(x):
    try:
        y=float(x); return y if math.isfinite(y) else None
    except Exception: return None

def norm_chrom(x):
    x=str(x)
    if x=="MT": return "chrM"
    return x if x.startswith("chr") else "chr"+x

def norm_type(x):
    x=str(x).upper()
    if "DEL" in x: return "DEL"
    if "DUP" in x or "GAIN" in x: return "DUP"
    return x

def overlap(a0,a1,b0,b1):
    lo=max(min(a0,a1),min(b0,b1)); hi=min(max(a0,a1),max(b0,b1))
    ov=max(0,hi-lo+1)
    if ov==0: return 0.0
    la=max(1,abs(a1-a0)+1); lb=max(1,abs(b1-b0)+1)
    return min(ov/la,ov/lb)

def parse_cnvpytor(path,max_q0,max_pn,max_eval):
    rows=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for i,line in enumerate(fh,1):
            if not line.strip() or line.startswith("#"): continue
            f=line.strip().split()
            if len(f)<4: continue
            typ=norm_type(f[0])
            if typ not in {"DEL","DUP"}: continue
            m=re.match(r"([^:]+):(\d+)-(\d+)",f[1].replace(",",""))
            if not m: continue
            chrom,start,end=norm_chrom(m.group(1)),int(m.group(2)),int(m.group(3))
            size=num(f[2]); level=num(f[3])
            eval1=num(f[4]) if len(f)>4 else None
            eval2=num(f[5]) if len(f)>5 else None
            eval3=num(f[6]) if len(f)>6 else None
            eval4=num(f[7]) if len(f)>7 else None
            q0=num(f[8]) if len(f)>8 else None
            pn=num(f[9]) if len(f)>9 else None
            dg=num(f[10]) if len(f)>10 else None
            flags=[]
            if q0 is not None and q0>max_q0: flags.append("HIGH_Q0")
            if pn is not None and pn>max_pn: flags.append("HIGH_N_CONTENT")
            if eval1 is not None and eval1>max_eval: flags.append("WEAK_RD_EVALUE")
            status="PASS" if not flags else "REVIEW"
            rows.append(dict(index=i,chrom=chrom,start=start,end=end,svtype=typ,size=size,
                             level=level,eval1=eval1,eval2=eval2,eval3=eval3,eval4=eval4,
                             q0=q0,pn=pn,dg=dg,status=status,flags=";".join(flags) if flags else "."))
    return rows

def parse_vcf(path):
    meta=[]; header=None; records=[]
    with opentext(path) as fh:
        for line in fh:
            if line.startswith("##"): meta.append(line); continue
            if line.startswith("#CHROM"): header=line; continue
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<8: continue
            inf=info_dict(f[7])
            svt=norm_type(inf.get("SVTYPE", f[4].strip("<>")))
            end=num(inf.get("END"))
            records.append({"fields":f,"info":inf,"chrom":norm_chrom(f[0]),"start":int(f[1]),
                            "end":int(end) if end is not None else int(f[1]),"svtype":svt})
    if header is None: raise ValueError("VCF lacks #CHROM header")
    return meta,header,records

def fmt(x):
    if x is None: return "."
    if isinstance(x,float):
        return f"{x:.8g}"
    return str(x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--breakpoint-vcf",required=True)
    ap.add_argument("--cnvpytor-calls",required=True)
    ap.add_argument("--reference-fai",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--audit-output",required=True)
    ap.add_argument("--min-reciprocal-overlap",type=float,default=0.5)
    ap.add_argument("--max-q0",type=float,default=0.5)
    ap.add_argument("--max-pn",type=float,default=0.5)
    ap.add_argument("--max-eval1",type=float,default=1e-4)
    args=ap.parse_args()

    meta,header,recs=parse_vcf(args.breakpoint_vcf)
    cnvs=parse_cnvpytor(args.cnvpytor_calls,args.max_q0,args.max_pn,args.max_eval1)

    matched=set()
    for r in recs:
        if r["svtype"] not in {"DEL","DUP"}: continue
        best=None; bestov=0.0
        for j,c in enumerate(cnvs):
            if c["chrom"]!=r["chrom"] or c["svtype"]!=r["svtype"]: continue
            ov=overlap(r["start"],r["end"],c["start"],c["end"])
            if ov>=args.min_reciprocal_overlap and ov>bestov:
                best,bestov=j,ov
        if best is not None:
            c=cnvs[best]; matched.add(best)
            add={
                "CNVPYTOR_MATCH":"YES",
                "CNVPYTOR_STATUS":c["status"],
                "CNVPYTOR_RECIPROCAL_OVERLAP":fmt(bestov),
                "CNVPYTOR_LEVEL":fmt(c["level"]),
                "CNVPYTOR_EVAL1":fmt(c["eval1"]),
                "CNVPYTOR_Q0":fmt(c["q0"]),
                "CNVPYTOR_PN":fmt(c["pn"]),
                "CNVPYTOR_FLAGS":c["flags"],
            }
            r["info"].update(add)
            r["fields"][7]=";".join([k if v=="True" else f"{k}={v}" for k,v in r["info"].items()])

    # Add quality-passing depth-only CNVs.
    sample_cols=max(0,len(header.rstrip("\n").split("\t"))-9)
    for j,c in enumerate(cnvs):
        if j in matched or c["status"]!="PASS": continue
        svlen=(c["end"]-c["start"]+1) * (-1 if c["svtype"]=="DEL" else 1)
        info={
            "SVTYPE":c["svtype"],"END":str(c["end"]),"SVLEN":str(svlen),
            "SRS_SOURCE":"CNVPytor","CNVPYTOR_ONLY":"1","CNVPYTOR_STATUS":c["status"],
            "CNVPYTOR_LEVEL":fmt(c["level"]),"CNVPYTOR_EVAL1":fmt(c["eval1"]),
            "CNVPYTOR_Q0":fmt(c["q0"]),"CNVPYTOR_PN":fmt(c["pn"]),"CNVPYTOR_FLAGS":c["flags"],
        }
        fields=[c["chrom"],str(c["start"]),f"CNVPYTOR_ONLY_{c['index']}","N",f"<{c['svtype']}>",".","PASS",
                ";".join(f"{k}={v}" for k,v in info.items())]
        if sample_cols:
            fields+=["GT"]+["./."]*sample_cols
        recs.append({"fields":fields,"info":info,"chrom":c["chrom"],"start":c["start"],"end":c["end"],"svtype":c["svtype"]})

    contig_order={}
    with open(args.reference_fai) as fh:
        for i,line in enumerate(fh):
            contig_order[line.split("\t",1)[0]]=i
    recs.sort(key=lambda r:(contig_order.get(r["chrom"],10**9),r["start"],r["end"],r["fields"][2]))

    defs=[
        '##INFO=<ID=SRS_SOURCE,Number=1,Type=String,Description="Additional SRS evidence source">\n',
        '##INFO=<ID=CNVPYTOR_ONLY,Number=0,Type=Flag,Description="Depth-only CNV added from CNVpytor">\n',
        '##INFO=<ID=CNVPYTOR_MATCH,Number=1,Type=String,Description="CNVpytor RD call overlaps breakpoint event">\n',
        '##INFO=<ID=CNVPYTOR_STATUS,Number=1,Type=String,Description="CNVpytor RD evidence QC status">\n',
        '##INFO=<ID=CNVPYTOR_RECIPROCAL_OVERLAP,Number=1,Type=Float,Description="Minimum reciprocal overlap with CNVpytor call">\n',
        '##INFO=<ID=CNVPYTOR_LEVEL,Number=1,Type=Float,Description="CNVpytor normalized read-depth level">\n',
        '##INFO=<ID=CNVPYTOR_EVAL1,Number=1,Type=Float,Description="CNVpytor primary RD e-value">\n',
        '##INFO=<ID=CNVPYTOR_Q0,Number=1,Type=Float,Description="CNVpytor q0 fraction">\n',
        '##INFO=<ID=CNVPYTOR_PN,Number=1,Type=Float,Description="CNVpytor reference N fraction">\n',
        '##INFO=<ID=CNVPYTOR_FLAGS,Number=1,Type=String,Description="CNVpytor review flags">\n',
    ]
    existing="".join(meta)
    meta += [x for x in defs if x.split("ID=",1)[1].split(",",1)[0] not in existing]

    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True)
    with open(out,"w",encoding="utf-8") as fh:
        fh.writelines(meta); fh.write(header if header.endswith("\n") else header+"\n")
        for r in recs: fh.write("\t".join(r["fields"])+"\n")

    audit=Path(args.audit_output); audit.parent.mkdir(parents=True,exist_ok=True)
    with open(audit,"w",encoding="utf-8") as fh:
        fh.write("CNVPYTOR_INDEX\tCHROM\tSTART\tEND\tSVTYPE\tSTATUS\tFLAGS\tMATCHED_TO_BREAKPOINT_UNIVERSE\n")
        for j,c in enumerate(cnvs):
            fh.write(f"{c['index']}\t{c['chrom']}\t{c['start']}\t{c['end']}\t{c['svtype']}\t{c['status']}\t{c['flags']}\t{'YES' if j in matched else 'NO'}\n")
    print(f"[OK] breakpoint_records={len(recs)-sum(1 for j,c in enumerate(cnvs) if j not in matched and c['status']=='PASS')} cnvpytor_calls={len(cnvs)} matched={len(matched)} output={out}")

if __name__=="__main__":
    main()
