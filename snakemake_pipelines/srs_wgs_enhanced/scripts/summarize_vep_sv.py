#!/usr/bin/env python3
"""Summarize all VEP SV transcript consequences per (SV_ID, gene).

Designed for VEP output produced with --flag_pick (not --pick), so all transcript
rows remain available while PICK/CANONICAL are retained as annotations.
"""
import argparse,csv,re
from collections import defaultdict
from pathlib import Path

MISSING="."

def parse_extra(text):
    out={}
    for item in str(text or "").split(";"):
        if "=" in item:
            k,v=item.split("=",1);out[k]=v
        elif item:
            out[item]="1"
    return out

def split_items(value):
    if value in (None,"",MISSING): return []
    return [x for x in re.split(r"[,;&]",str(value)) if x and x!="."]

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--vep",required=True)
    ap.add_argument("--output",required=True)
    a=ap.parse_args()

    header=None
    groups=defaultdict(list)
    with open(a.vep,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            if line.startswith("##"): continue
            if line.startswith("#Uploaded_variation"):
                header=line[1:].rstrip("\n").split("\t");continue
            if line.startswith("#") or not line.strip() or header is None: continue
            f=line.rstrip("\n").split("\t")
            row=dict(zip(header,f))
            extra=parse_extra(row.get("Extra",""))
            sid=row.get("Uploaded_variation",MISSING)
            gene=extra.get("SYMBOL",MISSING)
            if gene==MISSING:
                continue
            groups[(sid,gene.upper())].append((row,extra))

    cols=[
        "SV_ID","GENE","VEP_MATCH_METHOD","VEP_GENE_IDS","VEP_TRANSCRIPT_COUNT",
        "VEP_TRANSCRIPTS","VEP_WHOLE_TRANSCRIPT_COUNT","VEP_GENE_TRANSCRIPT_SCOPE",
        "VEP_CONSEQUENCES","VEP_IMPACTS","VEP_BIOTYPES","VEP_EXON","VEP_INTRON",
        "VEP_CANONICAL_TRANSCRIPTS","VEP_PICK_TRANSCRIPTS","VEP_OVERLAP_BP_MAX",
        "VEP_OVERLAP_PC_MAX","VEP_TRANSCRIPT_REGION_CLASS","VEP_STRUCTURAL_EFFECT"
    ]
    rows=[]
    for (sid,gene),items in sorted(groups.items()):
        genes=sorted({r.get("Gene",MISSING) for r,e in items if r.get("Gene",MISSING)!=MISSING})
        tx=sorted({r.get("Feature",MISSING) for r,e in items if r.get("Feature",MISSING)!=MISSING})
        cons=sorted({x for r,e in items for x in split_items(r.get("Consequence",MISSING))})
        impacts=sorted({e.get("IMPACT",MISSING) for r,e in items if e.get("IMPACT",MISSING)!=MISSING})
        biotypes=sorted({e.get("BIOTYPE",MISSING) for r,e in items if e.get("BIOTYPE",MISSING)!=MISSING})
        exons=sorted({r.get("EXON",e.get("EXON",MISSING)) for r,e in items if r.get("EXON",e.get("EXON",MISSING))!=MISSING})
        introns=sorted({r.get("INTRON",e.get("INTRON",MISSING)) for r,e in items if r.get("INTRON",e.get("INTRON",MISSING))!=MISSING})
        canonical=sorted({r.get("Feature",MISSING) for r,e in items if e.get("CANONICAL")=="YES"})
        picked=sorted({r.get("Feature",MISSING) for r,e in items if e.get("PICK") in {"1","YES"}})
        bp=[float(e[k]) for r,e in items for k in ("BP_OVERLAP","OVERLAP_BP") if e.get(k) not in (None,"",".") and re.fullmatch(r"[0-9.]+",e[k])]
        pc=[float(e[k]) for r,e in items for k in ("PERCENT_OVERLAP","OVERLAP_PC") if e.get(k) not in (None,"",".") and re.fullmatch(r"[0-9.]+",e[k])]
        region=[]
        for c in cons:
            if "transcript_ablation" in c: region.append("TRANSCRIPT_ABLATION")
            elif "exon_loss" in c: region.append("EXON_LOSS")
            elif "coding_sequence_variant" in c: region.append("CODING")
            elif "intron_variant" in c: region.append("INTRONIC")
            elif "upstream" in c or "downstream" in c: region.append("REGULATORY_PROXIMAL")
        rows.append({
            "SV_ID":sid,"GENE":gene,"VEP_MATCH_METHOD":"SV_ID_AND_SYMBOL",
            "VEP_GENE_IDS":";".join(genes) or MISSING,
            "VEP_TRANSCRIPT_COUNT":str(len(tx)),
            "VEP_TRANSCRIPTS":";".join(tx) or MISSING,
            "VEP_WHOLE_TRANSCRIPT_COUNT":str(len(tx)),
            "VEP_GENE_TRANSCRIPT_SCOPE":"ALL_VEP_TRANSCRIPT_ROWS_FLAG_PICK_NOT_FILTERED",
            "VEP_CONSEQUENCES":";".join(cons) or MISSING,
            "VEP_IMPACTS":";".join(impacts) or MISSING,
            "VEP_BIOTYPES":";".join(biotypes) or MISSING,
            "VEP_EXON":";".join(exons) or MISSING,
            "VEP_INTRON":";".join(introns) or MISSING,
            "VEP_CANONICAL_TRANSCRIPTS":";".join(canonical) or MISSING,
            "VEP_PICK_TRANSCRIPTS":";".join(picked) or MISSING,
            "VEP_OVERLAP_BP_MAX":str(max(bp)) if bp else MISSING,
            "VEP_OVERLAP_PC_MAX":str(max(pc)) if pc else MISSING,
            "VEP_TRANSCRIPT_REGION_CLASS":";".join(sorted(set(region))) or MISSING,
            "VEP_STRUCTURAL_EFFECT":";".join(cons) or MISSING,
        })
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    with open(out,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=cols,delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
    print(f"[OK] sv_gene_groups={len(rows)} output={out}")

if __name__=="__main__":
    main()
