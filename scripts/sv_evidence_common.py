"""Shared conservative parsers; no database absence or allele identity inferred."""
import math
import re

MISSING = {"", ".", "NA", "N/A", "NAN", "NONE", "NULL"}
AF_FIELDS = (
    "Allele_Freq_ALL_Control", "Pop_Freq_ALL", "Allele_Freq_ALL",
    "AlleleFreqAll", "AF_ALL", "AF", "MAX_AF", "AF_MAX", "SV_AF",
    "AF_1KGP", "1KGP_AF",
)


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def population_frequency(row):
    """Return (value, original field, status), with frequency units 0..1.

    A present but invalid higher-priority value is not silently replaced with a
    different population. Missing cells may fall through to a lower-priority field.
    """
    names = {k.lower(): k for k in row}
    for alias in AF_FIELDS:
        key = names.get(alias.lower())
        if key is None or str(row[key]).strip().upper() in MISSING:
            continue
        value = number(row[key])
        if value is None or not 0 <= value <= 1:
            return ".", key, "INVALID_FREQUENCY"
        return str(value), key, "REPORTED"
    return ".", ".", "UNKNOWN"


def breakend(row):
    """Read both endpoints and orientation. Unknown orientation stays unknown.

    Bracket orientation is encoded as local/remote retained-side signs. DELLY
    CT 3to5 corresponds to +-. Only matching representations are accepted.
    """
    chrom = str(row.get("CHROM", "."))
    pos = number(row.get("START", row.get("POS")))
    alt = str(row.get("ALT", "."))
    match = re.search(r"([\[\]])([^\[\]]+):(\d+)[\[\]]", alt)
    if match:
        chrom2, pos2 = match.group(2), int(match.group(3))
        orientation = ("-" if alt.startswith(("[", "]")) else "+") + (
            "+" if match.group(1) == "]" else "-"
        )
    else:
        chrom2 = str(row.get("CHR2", row.get("INFO_CHR2", ".")))
        pos2 = number(row.get("POS2"))
        if pos2 is None:
            pos2 = number(row.get("END"))
        orientation = str(row.get("BND_ORIENTATION", row.get("CALLER_STRANDS", row.get("INFO_CT", "."))))
        orientation = {"3to5": "+-", "5to3": "-+", "3to3": "++", "5to5": "--"}.get(orientation, orientation)
        if orientation not in {"++", "+-", "-+", "--"}:
            orientation = "."
    return chrom, pos, chrom2, pos2, orientation


def same_breakend(a, b, tolerance=500):
    a1, p1, a2, p2, ao = breakend(a)
    b1, q1, b2, q2, bo = breakend(b)
    if None in (p1, p2, q1, q2) or "." in (a2, b2, ao, bo):
        return None
    if a1 == b1 and a2 == b2 and ao == bo:
        distances = abs(p1-q1), abs(p2-q2)
    elif a1 == b2 and a2 == b1 and ao == bo[::-1]:
        distances = abs(p1-q2), abs(p2-q1)
    else:
        return None
    return sum(distances) if max(distances) <= tolerance else None
