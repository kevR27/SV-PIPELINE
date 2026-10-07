#!/usr/bin/env python3
"""Render the highest-priority Illumina SV calls with Samplot."""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path
import pandas as pd


def first(df, names):
    lower = {str(c).lower(): c for c in df.columns}
    return next((n if n in df.columns else lower[n.lower()] for n in names if n in df.columns or n.lower() in lower), None)


def number(value):
    try: return int(float(value))
    except (TypeError, ValueError): return None


def safe(value, length=55):
    return (re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "candidate")[:length]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("--input", required=True); p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True); p.add_argument("--sample", required=True); p.add_argument("--out-dir", required=True)
    p.add_argument("--top-n", type=int, default=20); p.add_argument("--min-mapq", type=int, default=20); p.add_argument("--window", type=int, default=5000)
    p.add_argument("--manifest", required=True); a = p.parse_args()
    df = pd.read_csv(a.input, sep="\t", dtype=str, low_memory=False)
    idc, chromc, startc = first(df, ["SV_ID", "ID"]), first(df, ["CHROM", "Chr"]), first(df, ["START", "POS"])
    endc, typec, genec = first(df, ["END"]), first(df, ["SVTYPE", "SV_type"]), first(df, ["GENE", "GENES", "gene"])
    if not all((idc, chromc, startc, typec)): raise ValueError("Input requires SV_ID, CHROM, START and SVTYPE")
    work = df.drop_duplicates(idc).head(a.top_n); outdir = Path(a.out_dir); outdir.mkdir(parents=True, exist_ok=True); manifest = []
    for rank, (_, row) in enumerate(work.iterrows(), 1):
        chrom, start, end = str(row[chromc]), number(row[startc]), number(row[endc]) if endc else None
        svtype, gene, sid = str(row[typec]).upper(), str(row[genec]) if genec else ".", str(row[idc])
        status, message, image = "SKIPPED", ".", "."
        if start is not None and svtype not in {"BND", "TRA", "CTX"}:
            end = end if end is not None and end > start else start + 1
            output = outdir / f"{rank:03d}_{safe(svtype,12)}_{safe(chrom,12)}_{start}_{end}_{safe(sid)}_{safe(gene,30)}.png"
            cmd = ["samplot", "plot", "-n", a.sample, "-b", a.bam, "-r", a.reference, "-o", str(output), "-c", chrom, "-s", str(start), "-e", str(end), "-q", str(a.min_mapq), "--window", str(a.window)]
            if svtype in {"DEL", "DUP", "INV"}: cmd += ["-t", svtype]
            result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            status = "OK" if result.returncode == 0 and output.exists() and output.stat().st_size else "FAILED"
            image = str(output) if status == "OK" else "."; message = "." if status == "OK" else (result.stdout or "").strip().split("\n")[-1]
        elif svtype in {"BND", "TRA", "CTX"}:
            message = "Interchromosomal event retained for manual two-breakpoint review"
        manifest.append({"RANK": rank, "SV_ID": sid, "GENE": gene, "SVTYPE": svtype, "CHROM": chrom, "START": start or ".", "END": end or ".", "STATUS": status, "IMAGE": image, "MESSAGE": message})
    target = Path(a.manifest); target.parent.mkdir(parents=True, exist_ok=True); pd.DataFrame(manifest).to_csv(target, sep="\t", index=False)
    print(f"[OK] requested={len(work)} rendered={sum(x['STATUS'] == 'OK' for x in manifest)} manifest={target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
