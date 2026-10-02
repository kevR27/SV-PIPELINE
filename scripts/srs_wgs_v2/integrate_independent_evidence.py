#!/usr/bin/env python3
"""Attach independent SRS evidence without redefining the master SV universe.

Inputs
------
1. The existing SRS allele-assessment table derived from Manta+DELLY+SURVIVOR.
2. CNVpytor calls (read-depth evidence).
3. GRIDSS VCF (breakpoint-assembly evidence).
4. Optional AnnotSV outputs for unmatched CNVpytor and GRIDSS calls.

Outputs
-------
- enhanced master table with CNVpytor and GRIDSS support columns;
- CNVpytor-only candidate table;
- GRIDSS-only candidate table.

GRIDSS and CNVpytor are NOT counted as equivalent callers. They represent
different evidence modalities.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

MISSING="."

def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path,"rt",encoding="utf-8",errors="replace")
    return open(path,"r",encoding="utf-8",errors="replace")

def to_int(v):
    try: return int(float(str(v).replace(",","")))
    except Exception: return None

def normalize_chrom(v):
    s=str(v or "")
    if not s: return s
    if s.startswith("chr"): return s
    if s=="MT": return "chrM"
    return "chr"+s

def normalize_type(v):
    s=str(v or MISSING).upper()
    if "BND" in s or s=="TRA": return "BND"
    for t in ("DEL","DUP","INV","INS","CNV"):
        if t in s: return t
    return s

def reciprocal_overlap(a1,a2,b1,b2):
    if None in (a1,a2,b1,b2): return 0.0
    alo,ahi=sorted((a1,a2)); blo,bhi=sorted((b1,b2))
    ov=max(0,min(ahi,bhi)-max(alo,blo)+1)
    if ov<=0: return 0.0
    return min(ov/max(1,ahi-alo+1), ov/max(1,bhi-blo+1))

def parse_info(raw):
    d={}
    for item in str(raw or "").split(";"):
        if not item: continue
        if "=" in item:
            k,v=item.split("=",1); d[k]=v
        else: d[item]="True"
    return d

def parse_bnd_alt(chrom,pos,alt):
    m=re.search(r"([\[\]])([^\[\]]+):(\d+)[\[\]]",str(alt))
    if not m: return None
    return normalize_chrom(chrom),to_int(pos),normalize_chrom(m.group(2)),int(m.group(3))

def read_table(path):
    if not path: return [],[]
    with open_text(path) as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        return list(rd), list(rd.fieldnames or [])

def find_col(cols,*names):
    low={c.lower():c for c in cols}
    for n in names:
        if n in cols: return n
        if n.lower() in low: return low[n.lower()]
    return None

def load_cnvpytor(path):
    rows,cols=read_table(path)
    for r in rows:
        r["_CHROM"]=normalize_chrom(r.get("CHROM"))
        r["_START"]=to_int(r.get("START"))
        r["_END"]=to_int(r.get("END"))
        r["_TYPE"]=normalize_type(r.get("SVTYPE"))
    return rows

def load_gridss(path):
    events=defaultdict(list)
    if not path: return {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            f=line.rstrip("\n").split("\t")
            if len(f)<8: continue
            chrom,pos,vid,ref,alt,qual,filt,info_raw=f[:8]
            info=parse_info(info_raw)
            event=info.get("EVENT") or info.get("MATEID") or vid
            endpoint=parse_bnd_alt(chrom,pos,alt)
            events[event].append({
                "ID":vid,"CHROM":normalize_chrom(chrom),"POS":to_int(pos),
                "ALT":alt,"QUAL":qual,"FILTER":filt,"INFO":info,"ENDPOINT":endpoint
            })
    out={}
    for event,recs in events.items():
        local=[]
        remote=[]
        for r in recs:
            local.append((r["CHROM"],r["POS"]))
            if r["ENDPOINT"]:
                remote.append((r["ENDPOINT"][2],r["ENDPOINT"][3]))
        pts=[]
        seen=set()
        for p in local+remote:
            if p[0] and p[1] is not None and p not in seen:
                seen.add(p); pts.append(p)
        quals=[]
        for r in recs:
            try: quals.append(float(r["QUAL"]))
            except Exception: pass
        simple_type=MISSING; simple_len=MISSING
        for r in recs:
            if r["INFO"].get("SIMPLE_SVTYPE"):
                simple_type=r["INFO"]["SIMPLE_SVTYPE"]
            if r["INFO"].get("SIMPLE_SVLEN"):
                simple_len=r["INFO"]["SIMPLE_SVLEN"]
        out[event]={
            "records":recs,"points":pts,
            "max_qual":max(quals) if quals else None,
            "simple_type":simple_type,"simple_len":simple_len,
        }
    return out

def annotsv_gene_index(path):
    if not path: return defaultdict(set), {}
    rows,cols=read_table(path)
    id_col=find_col(cols,"SV_ID","ID","AnnotSV_ID")
    gene_col=find_col(cols,"Gene_name","Gene","GENE","SYMBOL")
    class_col=find_col(cols,"ACMG_class","AnnotSV_class","Class")
    genes=defaultdict(set); classes={}
    if not id_col: return genes,classes
    for r in rows:
        vid=r.get(id_col,MISSING)
        if gene_col:
            for g in re.split(r"[,;/|]",str(r.get(gene_col,""))):
                g=g.strip()
                if g and g!=".": genes[vid].add(g)
        if class_col and r.get(class_col) not in (None,"","."):
            classes[vid]=r.get(class_col)
    return genes,classes

def gridss_matches_interval(ev,chrom,start,end,tol):
    if start is None or end is None: return False
    pts=[p for p in ev["points"] if p[0]==chrom and p[1] is not None]
    if len(pts)<2: return False
    poss=[p[1] for p in pts]
    return any(abs(a-start)<=tol and abs(b-end)<=tol for a in poss for b in poss if a!=b)

def gridss_matches_bnd(ev,chrom,start,chrom2,end2,tol):
    if None in (start,end2) or not chrom2: return False
    pts=ev["points"]
    has1=any(c==chrom and abs(p-start)<=tol for c,p in pts)
    has2=any(c==chrom2 and abs(p-end2)<=tol for c,p in pts)
    return has1 and has2

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--master",required=True)
    ap.add_argument("--cnvpytor")
    ap.add_argument("--cnvpytor-annotsv")
    ap.add_argument("--gridss")
    ap.add_argument("--gridss-annotsv")
    ap.add_argument("--reciprocal-overlap",type=float,default=0.5)
    ap.add_argument("--breakpoint-tolerance",type=int,default=500)
    ap.add_argument("--output",required=True)
    ap.add_argument("--cnvpytor-only",required=True)
    ap.add_argument("--gridss-only",required=True)
    args=ap.parse_args()

    rows,cols=read_table(args.master)
    c_chrom=find_col(cols,"CHROM","SV_chrom","Chr")
    c_start=find_col(cols,"START","SV_start","POS")
    c_end=find_col(cols,"END","SV_end")
    c_type=find_col(cols,"SVTYPE","SV_type","Type")
    c_chr2=find_col(cols,"CHR2")
    c_pos2=find_col(cols,"POS2")
    if not all((c_chrom,c_start,c_type)):
        raise ValueError("Master table requires CHROM, START/POS and SVTYPE")

    cnvs=load_cnvpytor(args.cnvpytor) if args.cnvpytor else []
    gridss=load_gridss(args.gridss) if args.gridss else {}
    used_cnv=set(); used_gridss=set()

    extra_cols=[
        "CNVPYTOR_RD_MATCH","CNVPYTOR_RD_MATCH_ID","CNVPYTOR_BIN",
        "CNVPYTOR_RD_LEVEL","CNVPYTOR_RECIPROCAL_OVERLAP","CNVPYTOR_EVAL1",
        "GRIDSS_BREAKPOINT_MATCH","GRIDSS_EVENT","GRIDSS_MAX_QUAL",
        "GRIDSS_SIMPLE_SVTYPE","GRIDSS_SIMPLE_SVLEN",
    ]

    for row in rows:
        chrom=normalize_chrom(row.get(c_chrom))
        start=to_int(row.get(c_start))
        end=to_int(row.get(c_end)) if c_end else start
        if end is None: end=start
        svtype=normalize_type(row.get(c_type))

        best=None
        cnv_matches=[]
        if svtype in {"DEL","DUP","CNV"} and start is not None and end is not None:
            for i,c in enumerate(cnvs):
                if c["_CHROM"]!=chrom: continue
                if svtype in {"DEL","DUP"} and c["_TYPE"]!=svtype: continue
                score=reciprocal_overlap(start,end,c["_START"],c["_END"])
                if score>=args.reciprocal_overlap:
                    cnv_matches.append((score,i,c))
                    if best is None or score>best[0]:
                        best=(score,i,c)
        if best:
            # A biological CNV may be emitted at more than one CNVpytor bin.
            # Mark every concordant multi-scale call as represented by this
            # master event; report the best-overlap row as the attached evidence.
            for _score,j,_call in cnv_matches:
                used_cnv.add(j)
            score,i,c=best
            row.update({
                "CNVPYTOR_RD_MATCH":"YES",
                "CNVPYTOR_RD_MATCH_ID":c.get("CNVPYTOR_ID",MISSING),
                "CNVPYTOR_BIN":c.get("BIN_SIZE",MISSING),
                "CNVPYTOR_RD_LEVEL":c.get("RD_LEVEL",MISSING),
                "CNVPYTOR_RECIPROCAL_OVERLAP":f"{score:.4f}",
                "CNVPYTOR_EVAL1":c.get("EVAL1",MISSING),
            })
        else:
            row.update({
                "CNVPYTOR_RD_MATCH":"NO" if args.cnvpytor else "NOT_RUN",
                "CNVPYTOR_RD_MATCH_ID":MISSING,"CNVPYTOR_BIN":MISSING,
                "CNVPYTOR_RD_LEVEL":MISSING,"CNVPYTOR_RECIPROCAL_OVERLAP":MISSING,
                "CNVPYTOR_EVAL1":MISSING,
            })

        matched=[]
        if args.gridss and start is not None:
            for event,ev in gridss.items():
                ok=False
                if svtype=="BND":
                    chr2=normalize_chrom(row.get(c_chr2)) if c_chr2 else ""
                    pos2=to_int(row.get(c_pos2)) if c_pos2 else end
                    ok=gridss_matches_bnd(ev,chrom,start,chr2,pos2,args.breakpoint_tolerance)
                elif svtype in {"DEL","DUP","INV","CNV"}:
                    ok=gridss_matches_interval(ev,chrom,start,end,args.breakpoint_tolerance)
                if ok:
                    matched.append((event,ev)); used_gridss.add(event)
        if matched:
            row["GRIDSS_BREAKPOINT_MATCH"]="YES"
            row["GRIDSS_EVENT"]=";".join(x[0] for x in matched)
            q=[x[1]["max_qual"] for x in matched if x[1]["max_qual"] is not None]
            row["GRIDSS_MAX_QUAL"]=str(max(q)) if q else MISSING
            row["GRIDSS_SIMPLE_SVTYPE"]=";".join(sorted({str(x[1]["simple_type"]) for x in matched if x[1]["simple_type"]!=MISSING})) or MISSING
            row["GRIDSS_SIMPLE_SVLEN"]=";".join(sorted({str(x[1]["simple_len"]) for x in matched if x[1]["simple_len"]!=MISSING})) or MISSING
        else:
            row["GRIDSS_BREAKPOINT_MATCH"]="NO" if args.gridss else "NOT_RUN"
            row["GRIDSS_EVENT"]=row["GRIDSS_MAX_QUAL"]=MISSING
            row["GRIDSS_SIMPLE_SVTYPE"]=row["GRIDSS_SIMPLE_SVLEN"]=MISSING

    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cols+extra_cols,delimiter="\t",extrasaction="ignore",lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    cnv_genes,cnv_classes=annotsv_gene_index(args.cnvpytor_annotsv)
    cnv_only_cols=["CNVPYTOR_ID","BIN_SIZE","CHROM","START","END","SVTYPE","SVLEN","RD_LEVEL","EVAL1","EVAL2","EVAL3","EVAL4","Q0","GENES","ANNOTSV_CLASS"]
    with open(args.cnvpytor_only,"w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=cnv_only_cols,delimiter="\t",lineterminator="\n"); w.writeheader()
        for i,c in enumerate(cnvs):
            if i in used_cnv: continue
            o={k:c.get(k,MISSING) for k in cnv_only_cols}
            vid=c.get("CNVPYTOR_ID",MISSING)
            o["GENES"]=";".join(sorted(cnv_genes.get(vid,set()))) or MISSING
            o["ANNOTSV_CLASS"]=cnv_classes.get(vid,MISSING)
            w.writerow(o)

    grid_genes,grid_classes=annotsv_gene_index(args.gridss_annotsv)
    grid_cols=["GRIDSS_EVENT","ENDPOINTS","MAX_QUAL","SIMPLE_SVTYPE","SIMPLE_SVLEN","RECORD_IDS","GENES","ANNOTSV_CLASS"]
    with open(args.gridss_only,"w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=grid_cols,delimiter="\t",lineterminator="\n"); w.writeheader()
        for event,ev in gridss.items():
            if event in used_gridss: continue
            ids=[r["ID"] for r in ev["records"]]
            genes=set()
            classes=[]
            for vid in ids:
                genes.update(grid_genes.get(vid,set()))
                if vid in grid_classes: classes.append(grid_classes[vid])
            w.writerow({
                "GRIDSS_EVENT":event,
                "ENDPOINTS":";".join(f"{c}:{p}" for c,p in ev["points"]),
                "MAX_QUAL":ev["max_qual"] if ev["max_qual"] is not None else MISSING,
                "SIMPLE_SVTYPE":ev["simple_type"],"SIMPLE_SVLEN":ev["simple_len"],
                "RECORD_IDS":";".join(ids),"GENES":";".join(sorted(genes)) or MISSING,
                "ANNOTSV_CLASS":";".join(sorted(set(classes))) or MISSING,
            })

    print(f"[OK] master_rows={len(rows)} cnvpytor_only={len(cnvs)-len(used_cnv)} gridss_only={len(gridss)-len(used_gridss)}")

if __name__=="__main__":
    main()
