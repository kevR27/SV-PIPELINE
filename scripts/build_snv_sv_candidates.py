#!/usr/bin/env python3
"""Report coding/splice small variants in genes that also contain candidate SVs.

This is a secondary recessive-gene review layer. It does not classify variants
or infer inheritance without appropriate evidence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path

import pandas as pd


KEEP_IMPACT = {"HIGH", "MODERATE"}
KEEP_CONSEQUENCE = {"splice_acceptor_variant", "splice_donor_variant"}


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")


def read_phase(path: str) -> dict[tuple[str, int, str], tuple[str, str, str]]:
    phase = {}
    with open_text(path) as handle:
        sample_index = None
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip().split("\t")
                sample_index = 9 if len(fields) > 9 else None
                continue
            if line.startswith("#") or sample_index is None:
                continue

            fields = line.rstrip().split("\t")
            chrom = fields[0].removeprefix("chr")
            pos = int(fields[1])
            alts = fields[4].split(",")
            fmt = fields[8].split(":")
            values = fields[sample_index].split(":")
            data = dict(zip(fmt, values))
            gt = data.get("GT", ".")
            ps = data.get("PS", ".")
            phased = "YES" if "|" in gt and ps not in {"", "."} else "NO"

            for alt in alts:
                phase[(chrom, pos, alt)] = (gt, ps, phased)
    return phase


def alt_haplotype(gt: str) -> str:
    if "|" not in gt:
        return "."
    alleles = gt.split("|")
    if alleles == ["1", "0"]:
        return "H1"
    if alleles == ["0", "1"]:
        return "H2"
    return "."


def read_vep(path: str) -> pd.DataFrame:
    header = None
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#Uploaded_variation"):
                header = line[1:].rstrip().split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            rows.append(dict(zip(header, line.rstrip().split("\t"))))
    return pd.DataFrame(rows)


def parse_location(value: str):
    text = str(value).split("-")[0]
    chrom, pos = text.split(":", 1)
    return chrom.removeprefix("chr"), int(pos)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vep", required=True)
    parser.add_argument("--phased-vcf", required=True)
    parser.add_argument("--sv-candidates", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    sv = pd.read_csv(args.sv_candidates, sep="\t", dtype=str, low_memory=False)
    genes = set(sv["GENE"].dropna().astype(str))
    phase = read_phase(args.phased_vcf)
    vep = read_vep(args.vep)

    if vep.empty:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        return

    small_rows = []
    for _, row in vep.iterrows():
        gene = str(row.get("SYMBOL", "."))
        if gene not in genes:
            continue

        consequence = str(row.get("Consequence", "."))
        impact = str(row.get("IMPACT", ".")).upper()
        consequences = set(consequence.split(","))
        if impact not in KEEP_IMPACT and not (consequences & KEEP_CONSEQUENCE):
            continue

        try:
            chrom, pos = parse_location(row.get("Location", "."))
        except Exception:
            continue
        alt = str(row.get("Allele", "."))
        gt, ps, phased = phase.get((chrom, pos, alt), (".", ".", "NO"))
        small_rows.append({
            "GENE": gene,
            "SMALL_VARIANT": row.get("Uploaded_variation", "."),
            "CHROM": chrom,
            "POS": pos,
            "ALT": alt,
            "CONSEQUENCE": consequence,
            "IMPACT": impact,
            "CLINVAR_OR_EXISTING_ID": row.get("Existing_variation", "."),
            "SMALL_GT": gt,
            "SMALL_PS": ps,
            "SMALL_PHASED": phased,
            "SMALL_ALT_HAPLOTYPE": alt_haplotype(gt),
        })

    small = pd.DataFrame(small_rows)
    output_rows = []

    for _, small_var in small.iterrows():
        for _, sv_row in sv[sv["GENE"].eq(small_var["GENE"])].iterrows():
            sv_gt = str(sv_row.get("LONGPHASE_GT", "."))
            sv_ps = str(sv_row.get("LONGPHASE_PS", "."))
            sv_haplotype = alt_haplotype(sv_gt)

            relation = "UNRESOLVED"
            if (
                small_var["SMALL_PHASED"] == "YES"
                and str(sv_row.get("LONGPHASE_PHASED", "")).upper() == "YES"
                and small_var["SMALL_PS"] not in {"", "."}
                and small_var["SMALL_PS"] == sv_ps
                and small_var["SMALL_ALT_HAPLOTYPE"] != "."
                and sv_haplotype != "."
            ):
                relation = "CIS" if small_var["SMALL_ALT_HAPLOTYPE"] == sv_haplotype else "TRANS"

            output_rows.append({
                **small_var.to_dict(),
                "SV_ID": sv_row["SV_ID"],
                "SVTYPE": sv_row["SVTYPE"],
                "SV_GENE_EFFECT": sv_row["SV_GENE_EFFECT"],
                "SV_POPULATION_STATUS": sv_row["POPULATION_STATUS"],
                "SV_CALL_SUPPORT": sv_row["CALL_SUPPORT"],
                "SV_GT": sv_gt,
                "SV_PS": sv_ps,
                "SV_ALT_HAPLOTYPE": sv_haplotype,
                "PHASE_RELATION": relation,
                "INTERPRETATION": (
                    "Same-gene small variant and SV. TRANS requires matching phase-set evidence; "
                    "UNRESOLVED does not imply cis or trans."
                ),
            })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(output_rows).to_csv(output, sep="\t", index=False)
    print(f"[OK] snv_sv_pairs={len(output_rows)} output={output}")


if __name__ == "__main__":
    main()
