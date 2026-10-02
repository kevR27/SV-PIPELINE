#!/usr/bin/env python3
"""Attach SRS-specific orthogonal evidence without changing the master SV rows.

Adds:
- CNVpytor RD evidence copied from the enhanced master VCF.
- GRIDSS breakpoint/assembly evidence by breakpoint proximity.
- MELT insertion evidence (when enabled).
- ExpansionHunter locus overlap (not counted as an SV caller).
- A descriptive TECHNICAL_EVIDENCE_CLASS. This is not ACMG classification.
"""
from __future__ import annotations
import argparse,csv,gzip,re
from collections import defaultdict
from pathlib import Path
MISSING="."

def opentext(p):
    return gzip.open(p,"rt",encoding="utf-8",errors="replace") if str(p).endswith(".gz") else open(p,encoding="utf-8",errors="replace")
def info(raw):
    d={}
    for x in str(raw).split(";"):
        if not x: continue
        if "=" in x:
            k,v=x.split("=",1); d[k]=v
        else:d[x]="True"
    return d
def ni(x):
    try:return int(float(str(x).replace(",","")))
    except:return None
def nf(x):
    try:return float(x)
    except:return None
def chrom(x):
    x=str(x)
    return "chrM" if x=="MT" else (x if x.startswith("chr") else "chr"+x)
def svtype(x):
    x=str(x or ".").upper()
    if "BND" in x or x=="TRA":return "BND"
    for t in ("DEL","DUP","INS","INV","CNV"):
        if t in x:return t
    return x
def parse_vcf(p):
    out=[]
    with opentext(p) as h:
        for line in h:
            if line.startswith("#"):continue
            f=line.rstrip().split("\t")
            if len(f)<8:continue
            d=info(f[7])
            out.append(dict(chrom=chrom(f[0]),pos=ni(f[1]),id=f[2],alt=f[4],filter=f[6],info=d,
                            sample=(dict(zip(f[8].split(":"),f[9].split(":"))) if len(f)>9 else {})))
    return out
def breakend(rec):
    m=re.search(r"([\[\]])([^\[\]]+):(\d+)[\[\]]",rec["alt"])
    if m:return rec["chrom"],rec["pos"],chrom(m.group(2)),int(m.group(3))
    return rec["chrom"],rec["pos"],chrom(rec["info"].get("CHR2",".")),ni(rec["info"].get("POS2",rec["info"].get("END")))
def interval_overlap(a,b,c,d):
    if None in (a,b,c,d):return 0
    return max(0,min(max(a,b),max(c,d))-max(min(a,b),min(c,d))+1)
def load_master(p):
    d={}
    for r in parse_vcf(p): d[r["id"]]=r
    return d
def load_gridss(p):
    by=defaultdict(list)
    if not p:return by
    for r in parse_vcf(p):
        a,p1,b,p2=breakend(r)
        if p1 is not None:by[a].append((p1,b,p2,r))
    return by
def load_melt(p):
    by=defaultdict(list)
    if not p:return by
    for r in parse_vcf(p):by[r["chrom"]].append(r)
    return by
def load_eh(p):
    by=defaultdict(list)
    for r in parse_vcf(p):
        r["end"]=ni(r["info"].get("END")) or r["pos"]
        by[r["chrom"]].append(r)
    return by
def first(row,*names):
    lower={k.lower():k for k in row}
    for n in names:
        k=n if n in row else lower.get(n.lower())
        if k and str(row.get(k,"")).strip() not in ("","."):return str(row[k])
    return "."
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--integrated",required=True);ap.add_argument("--master-vcf",required=True)
    ap.add_argument("--gridss");ap.add_argument("--melt");ap.add_argument("--expansionhunter",required=True)
    ap.add_argument("--breakpoint-tolerance",type=int,default=500);ap.add_argument("--output",required=True)
    a=ap.parse_args()
    with opentext(a.integrated) as h:
        rd=csv.DictReader(h,delimiter="\t");rows=list(rd);cols=list(rd.fieldnames or [])
    master=load_master(a.master_vcf);gridss=load_gridss(a.gridss);melt=load_melt(a.melt);eh=load_eh(a.expansionhunter)
    for row in rows:
        sid=first(row,"SV_ID","ID");m=master.get(sid)
        inf=m["info"] if m else {}
        row["CNVPYTOR_RD_MATCH"]=inf.get("CNVPYTOR_MATCH","NO" if "CNVPYTOR_ONLY" not in inf else "DEPTH_ONLY")
        row["CNVPYTOR_RD_STATUS"]=inf.get("CNVPYTOR_STATUS",".")
        row["CNVPYTOR_LEVEL"]=inf.get("CNVPYTOR_LEVEL",".")
        row["CNVPYTOR_EVAL1"]=inf.get("CNVPYTOR_EVAL1",".")
        row["CNVPYTOR_Q0"]=inf.get("CNVPYTOR_Q0",".")
        row["CNVPYTOR_PN"]=inf.get("CNVPYTOR_PN",".")
        row["CNVPYTOR_FLAGS"]=inf.get("CNVPYTOR_FLAGS",".")
        c=chrom(first(row,"CHROM","Chr"));s=ni(first(row,"START","POS","SV_start"));e=ni(first(row,"END","SV_end"))
        t=svtype(first(row,"SVTYPE","SV_type"))
        c2=chrom(first(row,"CHR2")) if first(row,"CHR2")!="." else c
        p2=ni(first(row,"POS2"))
        expected=[(c,s)]
        if t=="BND" and p2 is not None:expected.append((c2,p2))
        elif e is not None and e!=s:expected.append((c,e))
        hits=[]
        for ec,ep in expected:
            if ep is None:continue
            for gp,gc2,gp2,gr in gridss.get(ec,[]):
                if gp is not None and abs(gp-ep)<=a.breakpoint_tolerance:
                    hits.append(gr["id"] if gr["id"] not in ("",".") else f"{ec}:{gp}")
        row["GRIDSS_BREAKPOINT_MATCH"]="YES" if hits else ("NO" if a.gridss else "NOT_RUN")
        row["GRIDSS_BREAKPOINT_RECORDS"]=";".join(sorted(set(hits))) if hits else "."
        mm=[]
        if a.melt and t in {"INS","BND"} and s is not None:
            for r in melt.get(c,[]):
                if r["pos"] is not None and abs(r["pos"]-s)<=a.breakpoint_tolerance:mm.append(r["id"])
        row["MELT_MATCH"]="YES" if mm else ("NO" if a.melt else "NOT_RUN")
        row["MELT_INSERTIONS"]=";".join(sorted(set(mm))) if mm else "."
        hm=[]
        if s is not None:
            ee=e if e is not None and t!="BND" else s
            for r in eh.get(c,[]):
                if interval_overlap(s,ee,r["pos"],r["end"])>0:hm.append(r["info"].get("REPID",r["id"]))
        row["EXPANSIONHUNTER_LOCUS_OVERLAP"]="YES" if hm else "NO"
        row["EXPANSIONHUNTER_LOCI"]=";".join(sorted(set(hm))) if hm else "."
        cc=nf(first(row,"CALLER_COUNT","SUPP")) or 0
        rdpass=row["CNVPYTOR_RD_STATUS"]=="PASS" or row["CNVPYTOR_RD_MATCH"]=="DEPTH_ONLY"
        grid=row["GRIDSS_BREAKPOINT_MATCH"]=="YES"
        if "CNVPYTOR_ONLY" in inf:
            tech="RD_ONLY_CNV"
        elif cc>=2 and rdpass:
            tech="BREAKPOINT_CONCORDANT_PLUS_RD"
        elif cc>=2 and grid:
            tech="BREAKPOINT_CONCORDANT_PLUS_ASSEMBLY"
        elif cc>=2:
            tech="BREAKPOINT_CALLER_CONCORDANT"
        elif cc>=1 and (rdpass or grid):
            tech="SINGLE_BREAKPOINT_CALLER_PLUS_ORTHOGONAL"
        elif cc>=1:
            tech="SINGLE_BREAKPOINT_CALLER"
        else:
            tech="REVIEW_REQUIRED"
        row["TECHNICAL_EVIDENCE_CLASS"]=tech
        row["TECHNICAL_EVIDENCE_SCOPE"]="DESCRIPTIVE_RESEARCH_EVIDENCE_NOT_ACMG_CLASSIFICATION"
    new=["CNVPYTOR_RD_MATCH","CNVPYTOR_RD_STATUS","CNVPYTOR_LEVEL","CNVPYTOR_EVAL1","CNVPYTOR_Q0","CNVPYTOR_PN","CNVPYTOR_FLAGS",
         "GRIDSS_BREAKPOINT_MATCH","GRIDSS_BREAKPOINT_RECORDS","MELT_MATCH","MELT_INSERTIONS",
         "EXPANSIONHUNTER_LOCUS_OVERLAP","EXPANSIONHUNTER_LOCI","TECHNICAL_EVIDENCE_CLASS","TECHNICAL_EVIDENCE_SCOPE"]
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    with open(out,"w",newline="",encoding="utf-8") as h:
        w=csv.DictWriter(h,fieldnames=cols+[x for x in new if x not in cols],delimiter="\t",extrasaction="ignore",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print(f"[OK] rows={len(rows)} output={out}")
if __name__=="__main__":main()
