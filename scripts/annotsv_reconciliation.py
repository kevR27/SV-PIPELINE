"""Explain missing AnnotSV records using the log from the same master VCF."""
import re
from collections import Counter, defaultdict
from pathlib import Path

from sv_evidence_common import breakend, gene_symbols, invalid_gene_labels, number, sv_length


def chrom(value):
    return str(value).removeprefix("chr")


def read_skip_log(path, svs):
    """Validate line number, locus and type before attaching a log message.

    AnnotSV's reciprocal-breakend notices describe the remote endpoint. Never
    join a log by line number alone: a log from a different VCF could look valid.
    """
    if not path:
        return {}
    by_line = {r["_VCF_LINE"]: r for r in svs}
    notes = defaultdict(list)
    pattern = re.compile(r"(.+)_(-?\d+)_(-?\d+)_(DEL|DUP|INV|INS|TRA|BND|CNV)_(.*): (.+) \(line (\d+)\)$")
    for log_line, text in enumerate(Path(path).read_text().splitlines(), 1):
        if not text.strip():
            continue
        match = pattern.fullmatch(text)
        if not match:
            raise ValueError(f"Unrecognized AnnotSV skip-log format at line {log_line}")
        contig, start, end, kind, alleles, reason, input_line = match.groups()
        sv = by_line.get(int(input_line))
        normalized_kind = "BND" if kind == "TRA" else kind
        valid = sv is not None and normalized_kind == sv["SVTYPE"]
        if valid:
            if reason == "reciprocal breakend" and kind in {"BND", "TRA"}:
                _, _, remote, position, _ = breakend(sv)
                valid = chrom(remote) == contig and position == int(start)
            else:
                valid = chrom(sv["CHROM"]) == contig and number(sv["START"]) == int(start)
                if kind not in {"BND", "TRA"}:
                    literal = alleles.upper() == f"{sv['REF']}_{sv['ALT']}".upper()
                    # AnnotSV abbreviates long sequence-resolved deletions.
                    abbreviated = (kind == "DEL" and len(sv["REF"]) > 10000
                                   and alleles.upper() == f"{sv['REF'][0]}_<DEL>".upper()
                                   and number(sv.get("END")) == int(end))
                    valid = valid and (literal or abbreviated)
        if not valid:
            raise ValueError(f"AnnotSV log line {log_line} does not match master VCF line {input_line}; "
                             "supply the unchanged VCF used by that AnnotSV run")
        notes[sv["SV_ID"]].append(reason)
    return dict(notes)


def record_evidence(sv, annotations, notes, log_supplied=False):
    reasons = notes.get(sv["SV_ID"], [])
    if annotations:
        status = "ANNOTATED_WITH_LOG_NOTICE" if reasons else "ANNOTATED"
    elif any("< SVminSize" in r for r in reasons):
        status = "BELOW_ANNOTSV_MIN_SIZE"
    elif any(r.startswith('chromosome "') and r.endswith('unknown') for r in reasons):
        status = "UNRECOGNIZED_CONTIG"
    elif reasons:
        status = "SKIPPED_BY_ANNOTSV"
    else:
        status = "NO_ANNOTATION_NO_LOG_ENTRY" if log_supplied else "NO_ANNOTATION_LOG_NOT_SUPPLIED"
    invalid = sorted({g for r in annotations for g in invalid_gene_labels(r.get("Gene_name"))})
    has_genes = any(gene_symbols(r.get("Gene_name")) for r in annotations)
    mapping = ("UNRESOLVED_GENE_NAMES_PRESENT" if invalid else "ANNOTATED_GENE_NAMES" if has_genes
               else "NO_GENE_REPORTED" if annotations else "NOT_ANNOTATED")
    length, source = sv_length(sv)
    size = number(length)
    scope = ("BREAKEND_NO_SINGLE_LENGTH" if sv["SVTYPE"] in {"BND", "TRA"} else
             "UNKNOWN_LENGTH" if size is None else "BELOW_50_BP" if abs(size) < 50 else "AT_LEAST_50_BP")
    return {"SV_ID": sv["SV_ID"], "SV_SIZE_SCOPE": scope,
            "ANNOTSV_RECORD_STATUS": status, "ANNOTSV_UNANNOTATED_REASON": ";".join(reasons) or ".",
            "ANNOTSV_GENE_MAPPING_STATUS": mapping, "ANNOTSV_INVALID_GENE_LABELS": ";".join(invalid) or "."}


def summarize(records):
    return {"master_records": len(records),
            **{key: dict(Counter(r[key] for r in records)) for key in
               ("ANNOTSV_RECORD_STATUS", "SV_SIZE_SCOPE", "ANNOTSV_GENE_MAPPING_STATUS")}}
