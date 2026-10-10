"""Read event evidence without turning missing information into negative evidence.

These helpers do not recalculate ranks or infer pathogenicity. Compact final
tables and older detailed tables are both supported. Joins always use SV ID
and gene within one sample.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pandas as pd

from plot_utils import first_existing, numeric

MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}
PANEL_GROUPS = ["PANEL_GENE", "NONPANEL_GENE", "UNRESOLVED_PANEL"]
GROUP_LABELS = {"PANEL_GENE": "Panel", "NONPANEL_GENE": "Non-panel", "UNRESOLVED_PANEL": "Panel status unknown"}
GROUP_SUFFIXES = {"PANEL_GENE": "panel", "NONPANEL_GENE": "nonpanel", "UNRESOLVED_PANEL": "unresolved_panel"}
STATE_COLORS = {
    "REPORTED": "#D9EAF6",
    "SUPPORTING": "#B8DFD2",
    "REVIEW": "#F5D19B",
    "NO_MATCH": "#F5F5F5",
    "UNRESOLVED": "#E2D7EF",
    "UNAVAILABLE": "#D1D5DB",
    "NOT_APPLICABLE": "#FFF8E5",
}
STATE_LABELS = {
    "REPORTED": "Reported context",
    "SUPPORTING": "Computational support",
    "REVIEW": "Review flag",
    "NO_MATCH": "Evaluated: no compatible finding",
    "UNRESOLVED": "Unresolved",
    "UNAVAILABLE": "Not evaluated / unavailable",
    "NOT_APPLICABLE": "Not applicable",
}
DOMAINS = ["Gene relevance", "Gene effect", "Mechanism", "Call support", "Population", "Local depth", "SV phase", "SNV + SV", "Straglr", "TLDR", "Methylation"]
SRS_DOMAINS = ["Gene relevance", "Gene effect", "Mechanism", "Call support", "Population", "CNVpytor depth", "GRIDSS", "Repeat locus", "MELT"]


def domains_for(platform="lrs"):
    return SRS_DOMAINS if platform == "srs" else DOMAINS

# Only these small evidence columns are read from the large integrated table.
EVIDENCE_FIELDS = [
    "SV_ID", "ID", "GENE", "GENES", "Gene", "ANNotsv_Gene", "SVTYPE",
    "FINAL_GENE_RELEVANCE_TIER", "GENE_RELEVANCE_TIER", "GENE_RELEVANCE_DISPLAY_SCORE",
    "GENE_RELEVANCE_PHENOTYPE_SCOPE", "FINAL_PHENOTYPE_RANKING_SCOPE", "PHENOTYPE_SCORE_SCOPE",
    "SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP", "VEP_TRANSCRIPT_REGION_CLASS", "VEP_EXON", "VEP_TRANSCRIPTS", "VEP_CANONICAL_TRANSCRIPTS",
    "FINAL_INHERITANCE_MECHANISM_CLASS", "INHERITANCE_MECHANISM_CLASS", "INHERITANCE_MECHANISM_DETAIL",
    "CALLERS", "CALLER_COUNT", "SUPP", "READ_SUPPORT", "CALLER_READ_SUPPORT", "EVENT_MAX_CALLER_READ_SUPPORT",
    "CALLER_EVIDENCE_FLAGS", "CALLER_EVIDENCE_MATCH", "EVENT_TECHNICAL_REVIEW", "EVENT_TECHNICAL_TIER_V2",
    "NEEDLR_AF", "NEEDLR_STATUS", "POPULATION_STATUS", "GNOMAD_SV_AF", "GNOMAD_SV_EXACT_MATCH",
    "GNOMAD_SV_MATCH_SCOPE", "EVENT_POPULATION_AF_CONFLICT",
    "DEPTH_EVALUATION", "DEPTH_RATIO", "DEPTH_PATTERN", "DEPTH_SUPPORT_CLASS",
    "LONGPHASE", "LONGPHASE_MATCH", "LONGPHASE_GT", "LONGPHASE_PS", "LONGPHASE_PHASED",
    "RECESSIVE_PAIR_STATUS", "GENE_INHERITANCE_CLASS", "GENE_MOI_SET",
    "STRAGLR", "STRAGLR_MATCH", "STRAGLR_CONTEXT", "TLDR", "TLDR_MATCH",
    "METHYLATION_CONTEXT", "METHYLATION_5MC_MEAN_PERCENT", "ALLELE_GENOTYPE", "GENOTYPE",
    "DIAGNOSTIC_FOCUS", "TECHNICAL_SUPPORT", "CNVPYTOR_READ_DEPTH_MATCH",
    "CNVPYTOR_STATUS", "CNVPYTOR_LEVEL", "CNVPYTOR_FLAGS", "GRIDSS_SUPPORT",
    "MELT_SUPPORT", "EXPANSIONHUNTER_LOCUS_OVERLAP",
]


def text(row, names, default="."):
    for name in names:
        value = row.get(name)
        if value is not None and not pd.isna(value):
            value = str(value).strip()
            if value.upper() not in MISSING:
                return value
    return default


def number(row, names):
    value = text(row, names)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def normalise_events(df):
    work = df.copy()
    id_col = first_existing(work, ["SV_ID", "ID"])
    gene_col = first_existing(work, ["GENE", "GENES", "ANNotsv_Gene", "Gene"])
    if id_col is None or gene_col is None:
        raise ValueError("Expected an SV-gene event table with SV_ID and GENE/GENES. Do not pass the gene-level ranking.")
    work["SV_ID"] = work[id_col].fillna(".").astype(str).str.strip()
    work["GENE"] = work[gene_col].fillna(".").astype(str).str.strip()
    valid = ~work["SV_ID"].str.upper().isin(MISSING) & ~work["GENE"].str.upper().isin(MISSING)
    work = work.loc[valid].copy()
    panel_col = first_existing(work, ["PANEL_STATUS", "panel_gene"])
    raw = work[panel_col].fillna("").astype(str).str.strip().str.upper() if panel_col else pd.Series("", index=work.index)
    work["_panel_group"] = raw.map({
        "PANEL_GENE": "PANEL_GENE", "YES": "PANEL_GENE", "TRUE": "PANEL_GENE", "1": "PANEL_GENE",
        "NONPANEL_GENE": "NONPANEL_GENE", "NON_PANEL": "NONPANEL_GENE", "NO": "NONPANEL_GENE", "FALSE": "NONPANEL_GENE", "0": "NONPANEL_GENE",
    }).fillna("UNRESOLVED_PANEL")
    rank_col = first_existing(work, ["DIAGNOSTIC_REVIEW_RANK", "FINAL_EVENT_RANK_WITHIN_PANEL_STATUS", "EVENT_RANK_WITHIN_PANEL_STATUS", "FINAL_EVENT_RANK_GLOBAL", "EVENT_RANK_GLOBAL", "RANK_WITHIN_PANEL_STATUS"])
    score_col = first_existing(work, ["FINAL_GENE_RELEVANCE_DISPLAY_SCORE", "GENE_RELEVANCE_DISPLAY_SCORE", "GENE_RELEVANCE_SCORE", "EVENT_GENE_RELEVANCE_SCORE", "INTEGRATED_DISCOVERY_SCORE", "PHENOTYPE_SCORE"])
    work["_rank"] = numeric(work[rank_col]) if rank_col else float("nan")
    work["PLOT_ORDER_SCORE"] = numeric(work[score_col]) if score_col else float("nan")
    work["PLOT_ORDER_BASIS"] = rank_col or score_col or "SOURCE_ORDER_UNRANKED"
    work["_source_order"] = range(len(work))
    keys = ["_rank", "PLOT_ORDER_SCORE", "SV_ID", "GENE"] if rank_col or score_col else ["_source_order"]
    asc = [True, False, True, True] if rank_col or score_col else [True]
    return work.sort_values(keys, ascending=asc, na_position="last").drop_duplicates(["SV_ID", "GENE"]).reset_index(drop=True)


def select_events(df, top_n):
    if top_n < 1:
        raise ValueError("top-n must be at least 1")
    work = normalise_events(df)
    selected = [work[work["_panel_group"].eq(group)].head(top_n) for group in PANEL_GROUPS]
    result = pd.concat(selected, ignore_index=True)
    if result["PLOT_ORDER_BASIS"].eq("SOURCE_ORDER_UNRANKED").all():
        return result.sort_values("_source_order").reset_index(drop=True)
    return result.sort_values(["_rank", "PLOT_ORDER_SCORE", "SV_ID", "GENE"], ascending=[True, False, True, True], na_position="last").reset_index(drop=True)


def enrich_events(work, evidence_path):
    """Stream narrow evidence columns and fill missing selected-row fields only.

    Ambiguous SV-gene matches are reported, never resolved by choosing a row.
    Final ranks and patient-aware fields are not overwritten by upstream data.
    """
    work = work.copy()
    work["EVIDENCE_JOIN_STATUS"] = "NOT_REQUESTED"
    if not evidence_path or work.empty:
        return work
    header = pd.read_csv(evidence_path, sep="\t", nrows=0).columns
    id_col = first_existing(pd.DataFrame(columns=header), ["SV_ID", "ID"])
    gene_col = first_existing(pd.DataFrame(columns=header), ["GENE", "GENES", "ANNotsv_Gene", "Gene"])
    if id_col is None or gene_col is None:
        raise ValueError(f"Evidence table lacks an SV-gene key: {evidence_path}")
    columns = [c for c in header if c in EVIDENCE_FIELDS]
    wanted = set(zip(work["SV_ID"], work["GENE"]))
    matches = []
    for chunk in pd.read_csv(evidence_path, sep="\t", dtype=str, usecols=columns, chunksize=50000):
        keys = list(zip(chunk[id_col].fillna(".").str.strip(), chunk[gene_col].fillna(".").str.strip()))
        part = chunk.loc[[key in wanted for key in keys]].copy()
        if not part.empty:
            part["SV_ID"] = part[id_col].str.strip()
            part["GENE"] = part[gene_col].str.strip()
            matches.append(part)
    work["EVIDENCE_JOIN_STATUS"] = "NO_EXACT_SV_GENE_ROW"
    if not matches:
        return work
    detail = pd.concat(matches, ignore_index=True).drop_duplicates()
    duplicate = detail.duplicated(["SV_ID", "GENE"], keep=False)
    ambiguous = set(zip(detail.loc[duplicate, "SV_ID"], detail.loc[duplicate, "GENE"]))
    index = detail.loc[~duplicate].set_index(["SV_ID", "GENE"])
    for idx, row in work.iterrows():
        key = (row["SV_ID"], row["GENE"])
        if key in ambiguous:
            work.at[idx, "EVIDENCE_JOIN_STATUS"] = "AMBIGUOUS_SV_GENE_ROWS"
        elif key in index.index:
            work.at[idx, "EVIDENCE_JOIN_STATUS"] = "EXACT_SV_GENE_KEY"
            for col, value in index.loc[key].items():
                old = row.get(col)
                if old is None or pd.isna(old) or str(old).strip().upper() in MISSING:
                    work.at[idx, col] = value
    return work


def finding(value):
    value = str(value).upper()
    if value in {"YES", "TRUE", "1"}:
        return "REPORTED", "Match"
    if value in {"NO", "FALSE", "0", "NO_MATCH"}:
        return "NO_MATCH", "No match"
    if value == "NOT_APPLICABLE":
        return "NOT_APPLICABLE", "n/a"
    return "UNAVAILABLE", "Unavailable"


def af_values(value):
    values = []
    for token in re.split(r"[;,|]", str(value)):
        try:
            x = float(token.strip())
        except ValueError:
            continue
        if math.isfinite(x) and 0 <= x <= 1:
            values.append(x)
    return values


def population(row, rare_af=0.01):
    sources = []
    for col in ["NEEDLR_AF", "GNOMAD_SV_AF"]:
        values = af_values(row.get(col, "."))
        status = text(row, ["NEEDLR_STATUS", "POPULATION_STATUS"]).upper()
        if col == "NEEDLR_AF" and any(tag in status for tag in ["NO_MATCH", "NO_POPULATION_MATCH", "NOT_EVALUABLE", "NOT_AVAILABLE"]):
            continue
        # Explicitly rejected/absent gnomAD matches cannot supply an allele AF.
        if col == "GNOMAD_SV_AF" and text(row, ["GNOMAD_SV_EXACT_MATCH"]).upper() in {"NO", "NO_MATCH", "NOT_AVAILABLE"}:
            continue
        if values:
            sources.append((col, max(values)))
    if sources:
        maximum = max(x[1] for x in sources)
        state = "REVIEW" if maximum > rare_af else "REPORTED"
        label = f"max AF\n{maximum:.2g}"
        if text(row, ["EVENT_POPULATION_AF_CONFLICT"]).upper() in {"YES", "TRUE", "1"}:
            state = "REVIEW"
        return state, label
    status = text(row, ["NEEDLR_STATUS", "POPULATION_STATUS"])
    if "NO_MATCH" in status.upper() or "NO_POPULATION_MATCH" in status.upper():
        return "NO_MATCH", "No match\nAF unknown"
    if "NOT_EVALUABLE" in status.upper():
        return "UNAVAILABLE", "Not evaluable"
    return "UNAVAILABLE", "AF unknown"


def event_evidence(row, rare_af=0.01, platform="lrs"):
    """Return state and short observation for every domain; no total score."""
    result = {}
    tier = text(row, ["FINAL_GENE_RELEVANCE_TIER", "GENE_RELEVANCE_TIER"])
    result["Gene relevance"] = ("UNAVAILABLE", "Unknown") if tier == "." else ("REPORTED", tier.title())
    effect = text(row, ["SV_GENE_EFFECT", "SV_GENE_RELATIONSHIP"])
    result["Gene effect"] = ("UNRESOLVED", "Unresolved") if effect == "." or "UNRESOLVED" in effect else ("REPORTED", effect.replace("_", " ").lower())
    mechanism = text(row, ["FINAL_INHERITANCE_MECHANISM_CLASS", "INHERITANCE_MECHANISM_CLASS"])
    supported = {"STRONG_DOSAGE_OR_BIALLELIC_COMPATIBILITY", "SUPPORTED_DISEASE_MECHANISM", "AR_TRANS_SECOND_ALLELE_SUPPORTED"}
    if mechanism in supported:
        result["Mechanism"] = ("SUPPORTING", "Compatible")
    elif mechanism in MISSING or mechanism == "UNRESOLVED":
        result["Mechanism"] = ("UNRESOLVED", "Unresolved")
    else:
        result["Mechanism"] = ("REVIEW", mechanism.replace("_", " ").lower())
    count = number(row, ["CALLER_COUNT", "SUPP"])
    support = text(row, ["READ_SUPPORT", "CALLER_READ_SUPPORT", "EVENT_MAX_CALLER_READ_SUPPORT"])
    flags = text(row, ["CALLER_EVIDENCE_FLAGS"])
    match = text(row, ["CALLER_EVIDENCE_MATCH"])
    review = text(row, ["EVENT_TECHNICAL_REVIEW"])
    flagged = flags.upper() not in MISSING | {"PASS", "NONE", "NO_FLAGS"} or "AMBIGUOUS" in match or review.startswith("REVIEW_REQUIRED")
    if count is None and support == "." and flagged:
        result["Call support"] = ("REVIEW", "QC flag\nsupport unknown")
    elif count is None and support == ".":
        result["Call support"] = ("UNAVAILABLE", "Unknown")
    else:
        reads = number(row, ["EVENT_MAX_CALLER_READ_SUPPORT"])
        if reads is None:
            reads = number(row, ["READ_SUPPORT", "CALLER_READ_SUPPORT"])
        label = f"{count:g} callers" if count is not None else "Reads reported"
        if reads is not None:
            label += f"\nreads {reads:g}"
        if flagged:
            label += "\nQC review"
        result["Call support"] = ("REVIEW" if flagged else "REPORTED", label)
    result["Population"] = population(row, rare_af)
    svtype = text(row, ["SVTYPE"]).upper()
    depth_status = text(row, ["DEPTH_EVALUATION"])
    pattern = text(row, ["DEPTH_PATTERN"])
    ratio = number(row, ["DEPTH_RATIO"])
    if svtype in {"INV", "INS", "BND", "TRA"} or depth_status == "NOT_APPLICABLE":
        result["Local depth"] = ("NOT_APPLICABLE", "n/a")
    elif svtype not in {"DEL", "DUP", "CNV"} or depth_status in {"NOT_AVAILABLE", "NO_RESULT", "NOT_EVALUATED", "NO_DEPTH_SUMMARY"}:
        result["Local depth"] = ("UNAVAILABLE", "No depth result")
    elif pattern in {"CONSISTENT_WITH_LOSS", "CONSISTENT_WITH_GAIN"}:
        expected = "CONSISTENT_WITH_LOSS" if svtype == "DEL" else "CONSISTENT_WITH_GAIN" if svtype == "DUP" else pattern
        state = "SUPPORTING" if pattern == expected else "REVIEW"
        result["Local depth"] = (state, f"ratio {ratio:.2f}" if ratio is not None else pattern.replace("CONSISTENT_WITH_", "").title())
    elif ratio is not None or pattern == "NOT_CLEAR":
        result["Local depth"] = ("UNRESOLVED", f"ratio {ratio:.2f}" if ratio is not None else "Not clear")
    else:
        result["Local depth"] = ("UNAVAILABLE", "No depth result")
    phased = text(row, ["LONGPHASE_PHASED"])
    gt = text(row, ["LONGPHASE_GT"])
    ps = text(row, ["LONGPHASE_PS"])
    match = text(row, ["LONGPHASE", "LONGPHASE_MATCH"])
    if (match in {"NO", "NOT_AVAILABLE", "NO_MATCH"} or phased == "NO") and "|" in gt and ps != ".":
        result["SV phase"] = ("UNRESOLVED", "Phase/link conflict")
    elif (phased == "YES" or "|" in gt) and ps != ".":
        result["SV phase"] = ("REPORTED", gt if gt != "." else "Phased")
    elif match in {"YES", "NO"} or gt != "." or phased in {"YES", "NO"}:
        result["SV phase"] = ("UNRESOLVED", "Phase unresolved")
    else:
        result["SV phase"] = ("UNAVAILABLE", "Unavailable")
    pair = text(row, ["RECESSIVE_PAIR_STATUS"])
    inheritance = text(row, ["GENE_INHERITANCE_CLASS"])
    if pair in {"AR_TRANS_SNV_SV_CANDIDATE", "AR_TRANS_SECOND_ALLELE_SUPPORTED"}:
        result["SNV + SV"] = ("SUPPORTING", "Trans candidate")
    elif pair == "AR_CIS_NOT_BIALLELIC_BY_PHASE":
        result["SNV + SV"] = ("REVIEW", "Cis pair")
    elif pair == "NOT_AR_PAIRING_MODEL" or inheritance == "AUTOSOMAL_DOMINANT":
        result["SNV + SV"] = ("NOT_APPLICABLE", "n/a to AR model")
    elif "UNRESOLVED" in pair:
        result["SNV + SV"] = ("UNRESOLVED", "Phase unresolved")
    else:
        result["SNV + SV"] = ("UNAVAILABLE", "Not established")
    result["Straglr"] = finding(text(row, ["STRAGLR", "STRAGLR_MATCH"]))
    result["TLDR"] = finding(text(row, ["TLDR", "TLDR_MATCH"]))
    methyl = text(row, ["METHYLATION_CONTEXT"])
    if methyl == "EVALUATED":
        result["Methylation"] = ("REPORTED", "Context measured")
    elif methyl == "NO_CPG":
        result["Methylation"] = ("NO_MATCH", "No usable CpGs")
    else:
        result["Methylation"] = ("UNAVAILABLE", "Unavailable")
    if platform == "srs":
        # Use short-read evidence actually produced by this workflow. The
        # local small-variant phase is not an SV phase or an SV+SNV pairing.
        focus = text(row, ["DIAGNOSTIC_FOCUS"])
        if tier == "." and focus != ".":
            result["Gene relevance"] = ("REPORTED", focus.replace("_", " ").lower())
        depth_match = text(row, ["CNVPYTOR_READ_DEPTH_MATCH"]).upper()
        depth_qc = text(row, ["CNVPYTOR_STATUS"]).upper()
        level = number(row, ["CNVPYTOR_LEVEL"])
        if svtype in {"INV", "INS", "BND", "TRA"}:
            result["CNVpytor depth"] = ("NOT_APPLICABLE", "Not a depth CNV")
        elif depth_match in {"YES", "DEPTH_ONLY"}:
            direction_matches = level is not None and ((svtype == "DEL" and level < 1) or (svtype == "DUP" and level > 1))
            label = f"depth {level:.2f}" if level is not None else "Depth match"
            if depth_match == "DEPTH_ONLY":
                label += "\ndepth-only SV"
            result["CNVpytor depth"] = ("SUPPORTING" if depth_qc == "PASS" and direction_matches else "REVIEW", label)
        else:
            result["CNVpytor depth"] = finding(depth_match)
        gridss = text(row, ["GRIDSS_SUPPORT"]).upper()
        if gridss == "YES_BOTH_BREAKPOINTS":
            result["GRIDSS"] = ("SUPPORTING", "Both breakpoints")
        elif gridss == "YES_SINGLE_BREAKPOINT":
            result["GRIDSS"] = ("SUPPORTING", "Single-breakpoint event")
        elif gridss == "PARTIAL_ONE_BREAKPOINT":
            result["GRIDSS"] = ("UNRESOLVED", "One of two breakpoints")
        else:
            result["GRIDSS"] = finding(gridss)
        repeat = text(row, ["EXPANSIONHUNTER_LOCUS_OVERLAP"]).upper()
        result["Repeat locus"] = ("REPORTED", "Catalog locus overlap") if repeat == "YES" else finding(repeat)
        result["MELT"] = finding(text(row, ["MELT_SUPPORT"]))
        # A database search with no compatible match is different from a
        # resource that was not configured. Neither establishes novelty.
        scope = text(row, ["GNOMAD_SV_MATCH_SCOPE"]).upper()
        if "NOT_ATTEMPTED" in scope and result["Population"][0] == "UNAVAILABLE":
            result["Population"] = ("UNAVAILABLE", "Match not evaluated")
        elif text(row, ["GNOMAD_SV_EXACT_MATCH"]).upper() == "NO" and result["Population"][0] == "UNAVAILABLE":
            result["Population"] = ("NO_MATCH", "No compatible match\nAF unknown")
        return {domain: result[domain] for domain in SRS_DOMAINS}
    return result


def evidence_table(events, rare_af=0.01, platform="lrs"):
    rows = []
    for _, row in events.iterrows():
        for domain, (state, label) in event_evidence(row, rare_af, platform).items():
            rows.append({"SV_ID": row["SV_ID"], "GENE": row["GENE"], "PANEL_STATUS": row["_panel_group"], "DOMAIN": domain, "STATE": state, "OBSERVATION": label})
    return pd.DataFrame(rows, columns=["SV_ID", "GENE", "PANEL_STATUS", "DOMAIN", "STATE", "OBSERVATION"])
