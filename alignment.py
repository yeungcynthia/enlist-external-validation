"""Alignment and comparators extracted without algorithm changes from EnLiST.

See PROVENANCE.md for source revisions and extraction details.
"""


import re
from typing import Optional, List, Tuple, Dict
import numpy as np
import pandas as pd


GAP_PENALTY = -0.5


W_DRUG, W_START, W_END = 0.5, 0.3, 0.2


REGIMEN_LENIENT_JACCARD = 0.5


SYNONYMS = {
    "5fu": "5-fu", "5-fluorouracil": "5-fu", "fluorouracil": "5-fu", "5fu-": "5-fu",
    "5-fu/lv": "5-fu", "5fu/lv": "5-fu", "5-fu/leucovorin": "5-fu",
    "xeloda": "capecitabine", "cape": "capecitabine",
    "lv": "leucovorin", "folinic acid": "leucovorin",
    "avastin": "bevacizumab", "stivarga": "regorafenib", "erbitux": "cetuximab",
    "vectibix": "panitumumab", "keytruda": "pembrolizumab", "opdivo": "nivolumab",
    "yervoy": "ipilimumab", "cyramza": "ramucirumab", "zaltrap": "aflibercept",
    "tas-102": "trifluridine-tipiracil", "tas102": "trifluridine-tipiracil",
    "tas 102": "trifluridine-tipiracil", "lonsurf": "trifluridine-tipiracil",
    "trifluridine/tipiracil": "trifluridine-tipiracil",
    "trifluridine-tipiracil": "trifluridine-tipiracil", "trifluridine": "trifluridine-tipiracil",
    "medi4736": "durvalumab", "medi473": "durvalumab", "imfinzi": "durvalumab",
    "41bb": "pf-05082566", "pf-05082566": "pf-05082566",
    "zw25": "zanidatamab", "ds-8201a": "trastuzumab-deruxtecan", "ds8201a": "trastuzumab-deruxtecan",
    "ds-8201": "trastuzumab-deruxtecan", "enhertu": "trastuzumab-deruxtecan",
    "trastuzumab deruxtecan": "trastuzumab-deruxtecan", "trastuzumab-deruxtecan": "trastuzumab-deruxtecan",
    "agen2034": "balstilimab", "agen1181": "agen1181",
    "m7824": "bintrafusp-alfa", "bintrafusp alfa": "bintrafusp-alfa",
    "camptosar": "irinotecan", "irinotecan liposome": "irinotecan",
    "fruquintinib": "fruquintinib", "fruquintinib / placebo": "fruquintinib",
    "novo-ttf": "novo-ttf", "novo ttf": "novo-ttf", "ttfields": "novo-ttf",
    "ceralasertib": "ceralasertib", "azd6738": "ceralasertib",
    "ly3962673": "ly3962673", "kras g12d inhibitor ly3962673": "ly3962673",
    "kras g12d inhibitor": "ly3962673",
    "nab-paclitaxel": "nab-paclitaxel", "abraxane": "nab-paclitaxel",
    "gemcitabine": "gemcitabine", "gemzar": "gemcitabine",
    "mk-8353": "mk-8353", "mk8353": "mk-8353",
    "selumetinib": "selumetinib", "azd6244": "selumetinib",
    "cobimetinib": "cobimetinib", "cotellic": "cobimetinib",
    "lvgn6051": "lvgn6051",
}


REGIMEN_EXPANSIONS = {
    "folfox": {"5-fu", "leucovorin", "oxaliplatin"},
    "mfolfox6": {"5-fu", "leucovorin", "oxaliplatin"},
    "mfolfox-6": {"5-fu", "leucovorin", "oxaliplatin"},
    "modified folfox": {"5-fu", "leucovorin", "oxaliplatin"},
    "modified folfox 6": {"5-fu", "leucovorin", "oxaliplatin"},
    "folfox6": {"5-fu", "leucovorin", "oxaliplatin"},
    "folfiri": {"5-fu", "leucovorin", "irinotecan"},
    "mfolfiri": {"5-fu", "leucovorin", "irinotecan"},
    "modified folfiri": {"5-fu", "leucovorin", "irinotecan"},
    "folfoxiri": {"5-fu", "leucovorin", "oxaliplatin", "irinotecan"},
    "folfirinox": {"5-fu", "leucovorin", "oxaliplatin", "irinotecan"},
    "modified folfirinox": {"5-fu", "leucovorin", "oxaliplatin", "irinotecan"},
    "modified folfoxiri": {"5-fu", "leucovorin", "oxaliplatin", "irinotecan"},
    "capox": {"capecitabine", "oxaliplatin"}, "capeox": {"capecitabine", "oxaliplatin"},
    "xelox": {"capecitabine", "oxaliplatin"},
    "capiri": {"capecitabine", "irinotecan"}, "xeliri": {"capecitabine", "irinotecan"},
    "capeiri": {"capecitabine", "irinotecan"}, "irox": {"irinotecan", "oxaliplatin"},
}


def regimen_to_drug_set(regimen_str) -> Optional[frozenset]:
    """Decompose a regimen string into a canonical drug set (leucovorin stripped)."""
    if regimen_str is None or (isinstance(regimen_str, float) and pd.isna(regimen_str)):
        return None
    reg = str(regimen_str).strip().lower()
    if not reg:
        return None
    reg = re.sub(r"\([^)]*\)", "", reg)          # drop parentheticals
    for mod in ("adjuvant", "rechallenge", "monotherapy", "single agent", "modified", "maintenance", "neoadjuvant"):
        reg = re.sub(r"\b" + mod + r"\b", "", reg)
    # arrow notation -> take first stage only
    for arrow in ("\u2192", "->"):
        if arrow in reg:
            reg = reg.split(arrow)[0]
    reg = reg.replace(" + ", "+").replace(" / ", "+").replace("/", "+")
    reg = re.sub(r",\s*", "+", reg)
    reg = re.sub(r"\s+plus\s+", "+", reg)
    parts = [p.strip().strip("-").strip() for p in reg.split("+") if p.strip()]

    expanded = []
    for part in parts:
        if part in REGIMEN_EXPANSIONS or part in SYNONYMS:
            expanded.append(part)
        else:
            tokens = part.split()
            if len(tokens) > 1:
                joined = "".join(tokens)
                if joined in REGIMEN_EXPANSIONS or joined in SYNONYMS:
                    expanded.append(joined)
                elif len(tokens) == 2 and re.match(r"^\d+$", tokens[1]):
                    with_num = tokens[0] + tokens[1]
                    if with_num in REGIMEN_EXPANSIONS or with_num in SYNONYMS:
                        expanded.append(with_num)
                    elif tokens[0] in REGIMEN_EXPANSIONS or tokens[0] in SYNONYMS:
                        expanded.append(tokens[0])
                    else:
                        expanded.extend(tokens)
                else:
                    expanded.extend(tokens)
            else:
                expanded.append(part)

    drugs = set()
    for part in expanded:
        part = part.strip()
        if not part:
            continue
        if part in REGIMEN_EXPANSIONS:
            drugs.update(REGIMEN_EXPANSIONS[part])
        else:
            drugs.add(SYNONYMS.get(part, part))
    for junk in ("leucovorin", "placebo", "protocol", "trial", "inv-", ""):
        drugs.discard(junk)
    return frozenset(drugs) if drugs else None


def jaccard(a, b) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def date_proximity(a, b, max_days) -> float:
    if a is None or b is None:
        return 0.0
    diff = abs((a - b).days)
    return 0.0 if diff >= max_days else 1.0 - (diff / max_days)


def dates_match(a, b, tol) -> Optional[bool]:
    """None/None -> True; one None -> False; else within tolerance."""
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    return abs((a - b).days) <= tol


U_TOLERANCE = {"U0": 7, "U1": 14, "U2": 31, "U3": 93, "U4": 183}


U_TOLERANCE_DEFAULT = 30


def tol_for_u(u) -> int:
    """Map a gold U-code (U0..U4) to a date-match tolerance in days."""
    if u is None or (isinstance(u, float) and pd.isna(u)):
        return U_TOLERANCE_DEFAULT
    return U_TOLERANCE.get(str(u).strip().upper(), U_TOLERANCE_DEFAULT)


def norm_lot_label(label) -> Optional[str]:
    """'aLoT 1.0' / 'a1.0' / 'eLoT2.1' -> 'a1.0' / 'e2.1'."""
    if not label or pd.isna(label):
        return None
    s = str(label).strip().lower().replace(" ", "").replace("lot", "")
    return s or None


def _composite(a: Dict, b: Dict) -> float:
    drug_sim = jaccard(a["drugs"], b["drugs"]) if a["drugs"] and b["drugs"] else 0.0
    start_sim = date_proximity(a["start"], b["start"], 60)
    end_sim = date_proximity(a["stop"], b["stop"], 90)
    return W_DRUG * drug_sim + W_START * start_sim + W_END * end_sim


def _align_score(sim: float) -> float:
    if sim >= 0.5:
        return sim
    if sim >= 0.2:
        return sim * 0.5 - 0.1
    return -0.4


def needleman_wunsch(gs: List[Dict], llm: List[Dict],
                     gap: float = GAP_PENALTY) -> List[Tuple[Optional[int], Optional[int]]]:
    """Global alignment. Returns [(gs_idx|None, llm_idx|None), ...]."""
    n, m = len(gs), len(llm)
    if n == 0 and m == 0:
        return []
    if n == 0:
        return [(None, j) for j in range(m)]
    if m == 0:
        return [(i, None) for i in range(n)]

    S = np.zeros((n + 1, m + 1))
    for i in range(1, n + 1):
        S[i][0] = S[i - 1][0] + gap
    for j in range(1, m + 1):
        S[0][j] = S[0][j - 1] + gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            a = _align_score(_composite(gs[i - 1], llm[j - 1]))
            S[i][j] = max(S[i - 1][j - 1] + a, S[i - 1][j] + gap, S[i][j - 1] + gap)

    aln = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            a = _align_score(_composite(gs[i - 1], llm[j - 1]))
            cur = S[i][j]
            if abs(cur - (S[i - 1][j - 1] + a)) < 1e-9:
                aln.append((i - 1, j - 1)); i -= 1; j -= 1
            elif abs(cur - (S[i - 1][j] + gap)) < 1e-9:
                aln.append((i - 1, None)); i -= 1
            else:
                aln.append((None, j - 1)); j -= 1
        elif i > 0:
            aln.append((i - 1, None)); i -= 1
        else:
            aln.append((None, j - 1)); j -= 1
    aln.reverse()
    return aln
