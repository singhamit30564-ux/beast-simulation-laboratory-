"""CRISPR target search, guide design and specificity scoring.

Design notes
------------
* PAM matching is real IUPAC matching (see :data:`beastlab.sequtils.IUPAC`).
* Cut geometry follows the published behaviour of each effector: Cas9 family
  = blunt cut 3 bp upstream of a 3' PAM, Cas12 family = staggered cut 18/23 nt
  from a 5' PAM (5 nt 5' overhang).
* The *on-target* score is a transparent heuristic built from the sequence
  features that are known to matter (GC window, position-specific base
  preferences, homopolymers, self-complementarity, specificity).  It is **not**
  the Doench/Azimuth Rule Set 2 model, which is trained on measured efficiencies
  and cannot be reproduced from first principles.  The UI says so explicitly.
* The off-target score is an MIT/Hsu-style weighted-mismatch model using the
  published position-weight vector.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

from .sequtils import (IUPAC, clean, gc_content, longest_homopolymer, melting_temp,
                       mismatch_positions, pattern_regex, revcomp)

#: MIT/Hsu position weights for a 20 nt spacer, index 0 = PAM-proximal base.
MIT_WEIGHTS_20 = [0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.079,
                  0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

#: Modular scaffolds used to build the guide RNA shown in the UI.
SCAFFOLDS = {
    "sgRNA (SpCas9, 76 nt)":
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "tracrRNA (SpCas9, 63 nt)":
        "GGAACCAUUCAAAACAGCAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
    "crRNA repeat (S. pyogenes, 36 nt)":
        "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC",
    "crRNA repeat (S. thermophilus CRISPR1, 36 nt)":
        "GTTTTTGTACTCTCAAGATTTAAGTAACTGTACAAC",
    "crRNA repeat (S. thermophilus CRISPR3, 36 nt)":
        "GTTTTAGAGCTGTGTTGTTTCGAATGGTTCCAAAAC",
    "crRNA repeat (E. coli I-E, 29 nt)":
        "GGTTTATCCCCGCTGGCGCGGGGAACTC",
    "Cas12a direct repeat (LbCas12a, 20 nt)":
        "TAATTTCTACTAAGTGTAGAT",
    "Cas13a direct repeat (LwaCas13a, 30 nt, consensus)":
        "GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC",
    "pegRNA 3' extension motif (evopreQ1, RNA)":
        "GCGUCGCCAUGCGC",
}

#: Composition of the ubiquitous 76 nt sgRNA scaffold (segments, for drawing).
SGRNA_SEGMENTS = [
    ("Spacer (targeting)", 20, "#d4af37"),
    ("Tetraloop / nexus", 24, "#3fb6a8"),
    ("Stem-loop 1", 12, "#7f8cff"),
    ("Stem-loop 2", 10, "#ff7a59"),
    ("Stem-loop 3 / 3' tail", 10, "#9aa7bd"),
]


# --------------------------------------------------------------------------- #
# PAM scanning
# --------------------------------------------------------------------------- #
def pam_regex(pam_iupac: str) -> str:
    """Regex for an IUPAC PAM."""
    return pattern_regex(pam_iupac) if pam_iupac else ""


def scan_pams(seq: str, enzyme: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Find every PAM in `seq` for `enzyme` and return the targetable sites.

    Each site dict carries the protospacer, PAM, strand, 1-based genomic span and
    the cut coordinates implied by the effector's cut geometry.
    """
    s = clean(seq)
    sites: List[Dict[str, Any]] = []
    pam = enzyme.get("pam", "")
    if not pam or enzyme.get("pam_side") == "none":
        return sites
    regex = pam_regex(pam)
    side = enzyme["pam_side"]
    spacer_len = int(enzyme["spacer_len"])

    # plus strand
    import re
    for m in re.finditer(f"(?=({regex}))", s):
        p_start, p_end = m.start(1), m.end(1)
        if side == "3":
            ps_start, ps_end = p_start - spacer_len, p_start
            if ps_start < 0:
                continue
        else:  # 5'
            ps_start, ps_end = p_end, p_end + spacer_len
            if ps_end > len(s):
                continue
        sites.append(_site(seq, ps_start, ps_end, p_start, p_end, "+", enzyme))

    # minus strand: the PAM is read on the reverse complement
    rs = revcomp(s)
    for m in re.finditer(f"(?=({regex}))", rs):
        p_start, p_end = m.start(1), m.end(1)
        if side == "3":
            ps_start, ps_end = p_start - spacer_len, p_start
            if ps_start < 0:
                continue
        else:
            ps_start, ps_end = p_end, p_end + spacer_len
            if ps_end > len(s):
                continue
        # project back onto the plus strand
        f_ps_start = len(s) - ps_end
        f_ps_end = len(s) - ps_start
        f_p_start = len(s) - p_end
        f_p_end = len(s) - p_start
        sites.append(_site(seq, f_ps_start, f_ps_end, f_p_start, f_p_end, "-", enzyme))
    return sorted(sites, key=lambda d: d["start"])


def _site(seq: str, ps_start: int, ps_end: int, p_start: int, p_end: int,
          strand: str, enzyme: Dict[str, Any]) -> Dict[str, Any]:
    """Assemble one target-site record with cut coordinates."""
    s = clean(seq)
    protospacer = s[ps_start:ps_end] if strand == "+" else revcomp(s[ps_start:ps_end])
    pam_seq = s[p_start:p_end] if strand == "+" else revcomp(s[p_start:p_end])
    nts_cut, ts_cut = enzyme.get("nts_cut", 17), enzyme.get("ts_cut", 17)
    side = enzyme.get("pam_side", "3")
    if side == "3":
        # counted from the PAM-distal end of the protospacer
        cut_a = ps_start + nts_cut + 1
        cut_b = ps_start + ts_cut + 1
    else:
        # counted from the PAM-proximal end of the protospacer
        cut_a = ps_start + nts_cut + 1
        cut_b = ps_start + ts_cut + 1
    lo, hi = sorted((cut_a, cut_b))
    return {
        "start": ps_start, "end": ps_end,               # 0-based, half-open
        "start_1": ps_start + 1, "end_1": ps_end,       # 1-based inclusive
        "pam_start": p_start, "pam_end": p_end,
        "protospacer": protospacer,
        "pam": pam_seq,
        "strand": strand,
        "cut_left": lo, "cut_right": hi,                # 1-based cut boundaries
        "cut_site_1": (lo + hi) // 2,
        "blunt": enzyme.get("cut_geometry") == "blunt",
        "overhang": abs(cut_b - cut_a),
        "enzyme": enzyme["id"],
    }


def protospacer_for_guide(seq: str, guide: str, enzyme: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Locate a user-supplied spacer in `seq` and validate its PAM."""
    s = clean(seq)
    g = clean(guide)
    for strand in ("+", "-"):
        probe = g if strand == "+" else revcomp(g)
        start = s.find(probe)
        while start != -1:
            ps_start, ps_end = start, start + len(probe)
            if enzyme.get("pam_side") == "3":
                pam_zone = s[ps_end:ps_end + 8]
            else:
                pam_zone = s[max(0, ps_start - 8):ps_start]
            pam = enzyme.get("pam", "")
            if pam:
                import re
                hit = re.search(pam_regex(pam), pam_zone if enzyme.get("pam_side") == "3"
                                else revcomp(pam_zone))
                if hit:
                    if enzyme.get("pam_side") == "3":
                        p_start, p_end = ps_end + hit.start(), ps_end + hit.end()
                    else:
                        off = len(pam_zone) - hit.end()
                        p_start = ps_start - len(pam_zone) + off
                        p_end = ps_start - len(pam_zone) + off + (hit.end() - hit.start())
                    return _site(s, ps_start, ps_end, p_start, p_end, strand, enzyme)
            start = s.find(probe, start + 1)
    return None


# --------------------------------------------------------------------------- #
# on-target heuristics
# --------------------------------------------------------------------------- #
def _gc_component(guide: str) -> float:
    g = gc_content(guide)
    if 0.40 <= g <= 0.70:
        return 1.0
    if g < 0.40:
        return max(0.0, 1.0 - (0.40 - g) / 0.40)
    return max(0.0, 1.0 - (g - 0.70) / 0.30)


def _position_component(guide: str) -> float:
    """Preferred bases at the PAM-proximal and PAM-distal ends (empirical rules)."""
    g = clean(guide)
    if len(g) < 20:
        return 0.5
    score = 1.0
    if g[19] in "GC":                      # PAM-proximal position 20
        score += 0.05
    if g[18] in "GC":                      # position 19
        score += 0.05
    if g[0] == "G":                        # PAM-distal purine is favourable
        score += 0.05
    if g[1:4].count("G") >= 2:
        score -= 0.15
    if g[16] == "C":
        score -= 0.10
    return max(0.0, min(1.0, score / 1.15))


def _homopolymer_component(guide: str) -> float:
    base, run = longest_homopolymer(guide)
    if run >= 6:
        return 0.0
    if run == 5:
        return 0.35
    if run == 4:
        return 0.75
    return 1.0


def _polyt_component(guide: str) -> float:
    """U6-driven guides terminate at TTTT — a real design constraint."""
    g = clean(guide)
    return 0.0 if "TTTT" in g else 1.0


def _selfcomplementarity(guide: str) -> float:
    """Crude hairpin score: best reverse-complement match inside the spacer."""
    g = clean(guide)
    best = 0
    for i in range(len(g)):
        for j in range(i + 4, len(g) + 1):
            sub = g[i:j]
            if revcomp(sub) in g[:i] or revcomp(sub) in g[j:]:
                best = max(best, j - i)
    if best >= 8:
        return 0.0
    if best >= 6:
        return 0.4
    if best >= 5:
        return 0.7
    return 1.0


def specificity_score(off_targets: Sequence[Dict[str, Any]]) -> float:
    """MIT/Hsu-style specificity score in 0-100 (higher = more specific)."""
    total = 0.0
    for ot in off_targets:
        score = 100.0
        for pos in ot["mismatch_positions"]:
            idx = min(pos - 1, len(MIT_WEIGHTS_20) - 1)
            score *= (1.0 - MIT_WEIGHTS_20[idx])
        if ot.get("pam_valid") is False:
            score *= 0.2
        total += max(score, 0.0)
    return 100.0 * 100.0 / (100.0 + total)


def score_guide(guide: str, offtargets: Sequence[Dict[str, Any]] = ()) -> Dict[str, Any]:
    """Combine the sequence heuristics into one transparent on-target score."""
    g = clean(guide)
    spec = specificity_score(offtargets)
    comps = {
        "gc": _gc_component(g),
        "position": _position_component(g),
        "homopolymer": _homopolymer_component(g),
        "polyT": _polyt_component(g),
        "selfcomp": _selfcomplementarity(g),
        "specificity": min(1.0, spec / 70.0),
    }
    weights = {"gc": 0.30, "position": 0.20, "homopolymer": 0.15,
               "polyT": 0.10, "selfcomp": 0.10, "specificity": 0.15}
    total = sum(comps[k] * weights[k] for k in weights)
    seed_mm = sum(1 for ot in offtargets if ot["seed_mismatches"] == 0)
    return {
        "guide": g,
        "gc": gc_content(g),
        "tm": melting_temp(g),
        "on_target_score": round(100 * total, 1),
        "specificity_score": round(spec, 1),
        "components": {k: round(v, 3) for k, v in comps.items()},
        "n_offtargets": len(offtargets),
        "n_seed_perfect_offtargets": seed_mm,
        "homopolymer": longest_homopolymer(g)[1],
        "has_polyT": "TTTT" in g,
        "warnings": _warnings(g, seed_mm),
    }


def _warnings(guide: str, seed_perfect_offtargets: int) -> List[str]:
    out: List[str] = []
    if "TTTT" in guide:
        out.append("Contains TTTT — a U6-driven sgRNA may terminate early.")
    base, run = longest_homopolymer(guide)
    if run >= 5:
        out.append(f"Homopolymer run {base}×{run} — synthesis and Cas9 loading can suffer.")
    g = gc_content(guide)
    if g < 0.30:
        out.append("Very low GC — weak R-loop formation.")
    if g > 0.80:
        out.append("Very high GC — risk of secondary structure and slow turnover.")
    if seed_perfect_offtargets:
        out.append(f"{seed_perfect_offtargets} off-target site(s) with a perfect seed — high risk.")
    if guide.endswith("GG") or "GG" in guide[-4:]:
        out.append("G-rich PAM-proximal end can increase off-target tolerance.")
    return out


# --------------------------------------------------------------------------- #
# off-target search
# --------------------------------------------------------------------------- #
def find_offtargets(seq: str, guide: str, enzyme: Dict[str, Any],
                    max_mismatches: int = 4, include_pamless: bool = False) -> List[Dict[str, Any]]:
    """Scan `seq` for near-matches of `guide` (PAM-proximal seed first).

    A site is scored as an off-target when it carries the enzyme's PAM
    (``pam_valid=True``) or, when ``include_pamless`` is set, when it is a
    substitution-only near-match without a canonical PAM. Bulges are not modeled.
    """
    s, g = clean(seq), clean(guide)
    if not g or len(s) < len(g):
        return []
    seed_len = int(enzyme.get("seed_len", 12))
    side = enzyme.get("pam_side", "3")
    hits = []
    # Enumerate exact, adjacent PAM sites, then compare guide in its own orientation.
    sites = scan_pams(s, {**enzyme, "spacer_len": len(g)})
    if include_pamless:
        present = {(site["start"], site["strand"]) for site in sites}
        for strand in ("+", "-"):
            for start in range(len(s)-len(g)+1):
                if (start, strand) not in present:
                    window = s[start:start+len(g)]
                    sites.append({"start": start, "end": start+len(g), "strand": strand,
                                  "protospacer": window if strand == "+" else revcomp(window), "pam": ""})
    for site in sites:
        window = site["protospacer"]
        mm = mismatch_positions(g, window)
        # Existing ranking excludes all exact matches; external scanner reports them.
        if not mm or len(mm) > max_mismatches:
            continue
        seed_mm = sum(p > len(g)-seed_len if side == "3" else p <= seed_len for p in mm)
        hits.append({"start": site["start"], "end": site["end"], "strand": site["strand"],
                     "sequence": window, "mismatches": len(mm), "mismatch_positions": mm,
                     "seed_mismatches": seed_mm, "pam": site["pam"], "pam_valid": bool(site["pam"]),
                     "mismatch_string": "".join("." if a == b else "|" for a, b in zip(g, window))})
    return sorted(hits, key=lambda h: (h["mismatches"], h["seed_mismatches"], h["start"]))


# --------------------------------------------------------------------------- #
# guide design over a region
# --------------------------------------------------------------------------- #
def design_guides(seq: str, enzyme: Dict[str, Any], region: Optional[Tuple[int, int]] = None,
                  max_offtarget_mismatches: int = 4, limit: int = 200,
                  with_offtargets: bool = True, offtarget_top_n: int = 30,
                  reference_for_offtargets: Optional[str] = None) -> List[Dict[str, Any]]:
    """Design and score every guide in `region` (0-based half-open) of `seq`.

    Off-target scans are expensive, so a cheap heuristic pass ranks the candidates
    first and only the best ``offtarget_top_n`` are searched against the reference.
    """
    s = clean(seq)
    sites = scan_pams(s, enzyme)
    if region:
        lo, hi = region
        sites = [x for x in sites if lo <= x["start"] and x["end"] <= hi]
    ref = clean(reference_for_offtargets) if reference_for_offtargets else s
    prelim: List[Tuple[float, Dict[str, Any]]] = []
    for site in sites[:limit]:
        cheap = score_guide(site["protospacer"], [])
        prelim.append((cheap["on_target_score"], site))
    prelim.sort(key=lambda x: -x[0])
    if with_offtargets:
        prelim = prelim[:max(1, offtarget_top_n)]
    guides: List[Dict[str, Any]] = []
    for _, site in prelim:
        g = site["protospacer"]
        ot = find_offtargets(ref, g, enzyme, max_offtarget_mismatches) if with_offtargets else []
        scored = score_guide(g, ot)
        scored.update({
            "site": site,
            "cut_site_1": site["cut_site_1"],
            "strand": site["strand"],
            "pam": site["pam"],
            "start_1": site["start_1"],
            "end_1": site["end_1"],
            "n_candidate_sites": len(sites),
            "distance_to_region_start": site["start"] - (region[0] if region else 0),
        })
        guides.append(scored)
    guides.sort(key=lambda d: (-d["on_target_score"], d["n_offtargets"]))
    return guides


def pam_coverage(seq: str, pam_iupac: str, spacer_len: int = 20,
                 window: int = 50) -> Dict[str, Any]:
    """Fraction of positions in `seq` that have a targetable PAM within `window`."""
    s = clean(seq)
    regex = pam_regex(pam_iupac)
    import re
    pam_pos = [m.start() for m in re.finditer(regex, s)]
    pam_pos_rc = [len(s) - (m.start() + len(pam_iupac)) for m in re.finditer(regex, revcomp(s))]
    allp = sorted(pam_pos + pam_pos_rc)
    covered = 0
    dists: List[int] = []
    for i in range(len(s)):
        d = [abs(p - i) for p in allp]
        if d:
            dists.append(min(d))
            if min(d) <= window:
                covered += 1
    expected_gap = 4 ** len([c for c in pam_iupac if c == "N"]) * (3 / 4) if pam_iupac else 0
    return {
        "pam": pam_iupac,
        "n_pam_sites": len(allp),
        "pam_per_kb": 1000 * len(allp) / max(1, len(s)),
        "pam_every_bp": len(s) / max(1, len(allp)),
        "coverage_within_%d_bp" % window: (covered / max(1, len(s))),
        "median_distance_bp": sorted(dists)[len(dists) // 2] if dists else None,
        "random_expectation_every_bp": (4 ** pam_iupac.count("N")) / max(0.01, _pam_prob(pam_iupac)),
    }


def _pam_prob(pam: str) -> float:
    p = 1.0
    for c in pam:
        p *= len(IUPAC.get(c, "ACGT")) / 4.0
    return p


# --------------------------------------------------------------------------- #
# guide RNA construction
# --------------------------------------------------------------------------- #
def build_guide_rna(protospacer: str, enzyme: Dict[str, Any],
                    scaffold: str = "sgRNA (SpCas9, 76 nt)") -> Dict[str, Any]:
    """Assemble the RNA the effector actually loads."""
    ps = clean(protospacer)
    fam = enzyme.get("family", "")
    if fam.startswith("Cas12") or fam.startswith("Cas12f") or fam.startswith("Cas12j"):
        dr = "TAATTTCTACTAAGTGTAGAT"
        crrna = dr + ps
        return {"name": "crRNA (single RNA)", "rna": crrna, "spacer": ps,
                "repeat": dr, "len": len(crrna),
                "note": "Cas12a-type: one crRNA, 5' direct repeat + spacer, no tracrRNA."}
    if fam.startswith("Cas13"):
        dr = "GATTTAGACTACCCCAAAAACGAAGGGGACTAAAAC"
        crrna = dr + ps
        return {"name": "crRNA (Cas13)", "rna": crrna, "spacer": ps, "repeat": dr,
                "len": len(crrna),
                "note": "Type VI: targets RNA, no tracrRNA and no DNA PAM (PFS instead)."}
    scaf = SCAFFOLDS.get(scaffold, SCAFFOLDS["sgRNA (SpCas9, 76 nt)"])
    sgrna = ps + scaf
    return {"name": "sgRNA (spacer + scaffold)", "rna": sgrna, "spacer": ps,
            "scaffold": scaf, "len": len(sgrna),
            "note": "Type II: crRNA fused to the tracrRNA scaffold — the classic single guide."}


# --------------------------------------------------------------------------- #
# donor template / HDR design
# --------------------------------------------------------------------------- #
def design_hdr_donor(seq: str, cut_site_1: int, edit: str,
                     arm: int = 60, edit_offset: Optional[int] = None) -> Dict[str, Any]:
    """Build a symmetric HDR donor around `cut_site_1` (1-based) that installs `edit`.

    ``edit`` replaces ``[cut_site_1 + edit_offset, ... + len(edit))``.
    """
    s = clean(seq)
    cut0 = cut_site_1 - 1
    off = edit_offset if edit_offset is not None else 0
    e_start = cut0 + off
    e_end = e_start + len(edit)
    left = s[max(0, cut0 - arm):e_start]
    right = s[e_end:min(len(s), cut0 + arm)]
    donor = left + clean(edit) + right
    return {
        "donor": donor, "left_arm": left, "right_arm": right,
        "edit_string": edit, "length": len(donor),
        "gc": gc_content(donor),
        "left_tm": melting_temp(left[-25:]) if left else 0.0,
        "right_tm": melting_temp(right[:25]) if right else 0.0,
        "silent_pam_mutation": _pam_mutation_hint(s, cut0),
        "note": "Blocking the PAM or the seed in the donor prevents re-cutting of the edited allele.",
    }


def _pam_mutation_hint(seq: str, cut0: int) -> str:
    """Suggest a synonymous PAM disruption by mutating the NGG's second G."""
    return "e.g. change the PAM GG -> GA/GC (or the seed) in the donor to stop re-cutting."


# --------------------------------------------------------------------------- #
# prime editing design
# --------------------------------------------------------------------------- #
def design_pegRNA(seq: str, protospacer: str, enzyme: Dict[str, Any],
                  edit: str, edit_offset: int = 17,
                  pbs_len: int = 13, rtt_len: int = 16,
                  add_motif: bool = True) -> Dict[str, Any]:
    """Design a pegRNA: spacer + scaffold + 3' extension (PBS + RT template).

    ``edit`` is written into the RT template at ``edit_offset`` nt downstream of
    the nick (default = the Cas9 nick position, i.e. 3 bp 5' of the PAM).
    """
    s = clean(seq)
    ps = clean(protospacer)
    start = s.find(ps)
    if start == -1:
        start = s.find(revcomp(ps))
    if start == -1:
        raise ValueError("protospacer not found in the supplied sequence")
    nick = start + (enzyme.get("ts_cut", 17) + 1)  # 0-based nick position
    pbs = revcomp(s[max(0, nick - pbs_len):nick])
    rtt_region = list(s[nick:nick + max(rtt_len, edit_offset + len(edit))])
    while len(rtt_region) < edit_offset + len(edit):
        rtt_region.append("N")
    rtt_region[edit_offset:edit_offset + len(edit)] = list(clean(edit))
    rtt = "".join(rtt_region)
    motif = SCAFFOLDS["pegRNA 3' extension motif (evopreQ1, RNA)"] if add_motif else ""
    scaffold = SCAFFOLDS["sgRNA (SpCas9, 76 nt)"]
    return {
        "spacer": ps, "scaffold": scaffold, "pbs": pbs, "pbs_len": len(pbs),
        "pbs_tm": melting_temp(pbs),
        "pbs_gc": gc_content(pbs),
        "rtt": rtt, "rtt_len": len(rtt),
        "motif": motif,
        "pegRNA": ps + scaffold + rtt + motif,
        "nick_1": nick + 1,
        "note": "PBS 8-15 nt (GC 40-60% is comfortable); keep the 3' end of the PBS "
                "from ending in a run of Ts. PE3 adds a second nick ~50 bp from the edit.",
    }


def recommend_pbs(rtt_target: str, minimum: float = 40.0) -> List[Tuple[int, float, str]]:
    """Rank PBS lengths by melting temperature for the given RT-template 3' end."""
    out = []
    for L in range(8, 16):
        pbs = revcomp(clean(rtt_target)[:L])
        out.append((L, melting_temp(pbs), pbs))
    return sorted(out, key=lambda x: abs(x[1] - max(minimum, 45.0)))
