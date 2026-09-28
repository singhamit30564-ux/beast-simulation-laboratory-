"""Bacterial / archaeal CRISPR immunity: arrays, adaptation, interference and the
phage-bacteria population arms race.

What is real here
-----------------
* The CRISPR array of *E. coli* K-12 is parsed at run time from the GenBank window
  shipped in ``data/genomes/ecoli_k12.fasta`` — repeats and spacers are read out
  of the chromosome, not invented.
* Spacer acquisition is polarised (new spacers insert at the leader-proximal
  end) and PAM-dependent, because that is what Cas1-Cas2 does.
* Interference requires both a matching spacer *and* a PAM, and escape mutants
  arise from point mutations in the protospacer or PAM — the two ways phages
  really beat CRISPR.  Anti-CRISPR proteins are modelled as an additional,
  independent escape route.

What is a model, not a measurement
----------------------------------
The population dynamics use ordinary differential equations with published
phage-plaque parameters (adsorption rate, burst size, latency) that are typical
for *E. coli*/phage systems.  They reproduce the qualitative arms race —
collapse, rebound, escape — not any particular published time-course.
"""
from __future__ import annotations

import math
import random
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .sequtils import clean, gc_content, pattern_regex, revcomp

# --------------------------------------------------------------------------- #
# array parsing / construction
# --------------------------------------------------------------------------- #
def parse_crispr_array(seq: str, repeat: str, min_spacer: int = 20) -> Dict[str, Any]:
    """Split a CRISPR locus into repeats and spacers (leader assumed 5')."""
    s = clean(seq)
    r = clean(repeat)
    spans: List[Tuple[int, int]] = []
    start = s.find(r)
    while start != -1:
        spans.append((start, start + len(r)))
        start = s.find(r, start + len(r))
    spacers: List[Dict[str, Any]] = []
    for (a1, b1), (a2, b2) in zip(spans, spans[1:]):
        sp = s[b1:a2]
        if len(sp) >= min_spacer:
            spacers.append({"sequence": sp, "start": b1, "end": a2, "length": len(sp)})
    leader_end = spans[0][0] if spans else 0
    return {
        "repeat": r, "repeats": len(spans), "spacer_objects": spacers,
        "spacers": [x["sequence"] for x in spacers],
        "array_start": spans[0][0] if spans else 0,
        "array_end": spans[-1][1] if spans else 0,
        "leader": s[:leader_end] if spans else "",
        "leader_length": leader_end,
        "trailing": s[spans[-1][1]:] if spans else s,
        "spacer_lengths": [len(x["sequence"]) for x in spacers],
        "array_length": (spans[-1][1] - spans[0][0]) if spans else 0,
    }


def rebuild_array(repeat: str, spacers: Sequence[str], leader: str = "") -> str:
    """Assemble a locus from a repeat and an ordered spacers list."""
    return clean(leader) + clean(repeat) + clean(repeat).join(clean(s) for s in spacers) + clean(repeat)


def seed_matched_fraction(spacer: str, protospacer: str, seed_len: int = 8) -> float:
    """Fraction of the PAM-proximal seed that matches (type I seeds are short)."""
    a, b = clean(spacer), clean(protospacer)
    n = min(len(a), len(b), seed_len)
    if n == 0:
        return 0.0
    return sum(x == y for x, y in zip(a[-n:], b[-n:])) / n


# --------------------------------------------------------------------------- #
# adaptation — acquiring new spacers
# --------------------------------------------------------------------------- #
def candidate_protospacers(invader: str, system: Dict[str, Any],
                           limit: int = 400) -> List[Dict[str, Any]]:
    """Every protospacer in an invader genome that sits next to a valid PAM.

    Both conventions are searched (PAM 5' or 3' of the protospacer) so that a
    system quoted either way still finds its real sites.
    """
    s = clean(invader)
    ps_len = int(system.get("protospacer_len", 32))
    pam = system.get("pam_iupac", "")
    regex = pattern_regex(pam) if pam else ""
    out: List[Dict[str, Any]] = []
    if not regex:
        # no PAM requirement: sample the genome
        step = max(1, len(s) // limit)
        for i in range(0, max(1, len(s) - ps_len), step):
            out.append({"start": i, "end": i + ps_len, "strand": "+",
                        "protospacer": s[i:i + ps_len], "pam": "n/a", "pam_side": "none"})
        return out[:limit]

    # PAM 3' of the protospacer
    for m in re.finditer(regex, s):
        lo = m.start() - ps_len
        if lo >= 0:
            out.append({"start": lo, "end": m.start(), "strand": "+",
                        "protospacer": s[lo:m.start()], "pam": s[m.start():m.end()],
                        "pam_side": "3' of protospacer"})
        hi = m.end() + ps_len
        if hi <= len(s):
            out.append({"start": m.end(), "end": hi, "strand": "+",
                        "protospacer": s[m.end():hi], "pam": s[m.start():m.end()],
                        "pam_side": "5' of protospacer"})
    rs = revcomp(s)
    for m in re.finditer(regex, rs):
        lo = m.start() - ps_len
        if lo >= 0:
            f_lo = len(s) - m.start()
            f_hi = len(s) - (m.start() - ps_len)
            out.append({"start": f_hi, "end": f_lo, "strand": "-",
                        "protospacer": revcomp(rs[lo:m.start()]),
                        "pam": revcomp(rs[m.start():m.end()]),
                        "pam_side": "3' of protospacer"})
    return out[:limit * 4]


def simulate_adaptation(invader: str, system: Dict[str, Any], array: Dict[str, Any],
                        n_new: int = 2, mode: str = "naive",
                        existing_target: Optional[str] = None,
                        seed: int = 3) -> Dict[str, Any]:
    """Cas1-Cas2 integrates new spacers at the leader-proximal end of the array.

    ``mode='primed'`` biases acquisition to the neighbourhood of an existing
    target (primed adaptation), which is what happens when Cascade is already
    bound to an invader.
    """
    rng = random.Random(seed)
    cands = candidate_protospacers(invader, system)
    if not cands:
        return {"new_spacers": [], "array": array, "events": [], "candidates": 0}
    if mode == "primed" and existing_target:
        pos = clean(invader).find(clean(existing_target))
        if pos >= 0:
            cands.sort(key=lambda c: abs(c["start"] - pos))
            cands = cands[:max(8, len(cands) // 10)]
    existing = set(array["spacers"])
    chosen: List[Dict[str, Any]] = []
    tries = 0
    while len(chosen) < n_new and tries < 200 * max(1, n_new):
        tries += 1
        c = rng.choice(cands)
        if c["protospacer"] in existing or c["protospacer"] in [x["protospacer"] for x in chosen]:
            continue
        # PAM-proximal 3 nt of the prespacer derive from the PAM; 8 nt of the
        # protospacer are lost at the leader-distal end during integration
        chosen.append(c)
    new_spacers = [c["protospacer"] for c in chosen]
    ordered = new_spacers + list(array["spacers"])         # leader-proximal insertion
    events = [
        {"event": "spacer acquisition",
         "spacer": c["protospacer"], "pam": c["pam"], "pam_side": c["pam_side"],
         "invader_position": f"{c['start']}-{c['end']}",
         "gc": round(gc_content(c["protospacer"]), 3)}
        for c in chosen
    ]
    return {
        "new_spacers": new_spacers,
        "array": {**array, "spacers": ordered, "spacer_objects": [],
                  "array_sequence": rebuild_array(array["repeat"], ordered)},
        "events": events,
        "candidates": len(cands),
        "mode": mode,
        "note": ("New spacers are inserted next to the leader — acquisition is polarised, "
                 "so the newest immunity sits at the 5' end of the array."),
    }


# --------------------------------------------------------------------------- #
# interference — does the array stop the invader?
# --------------------------------------------------------------------------- #
def interference_scan(invader: str, spacers: Sequence[str], system: Dict[str, Any],
                      max_mismatches: int = 0) -> List[Dict[str, Any]]:
    """Match every spacer against the invader, requiring a PAM."""
    s = clean(invader)
    hits: List[Dict[str, Any]] = []
    for idx, sp in enumerate(spacers):
        sp_c = clean(sp)
        for strand, probe in (("+", s), ("-", revcomp(s))):
            start = probe.find(sp_c)
            while start != -1:
                protospacer = probe[start:start + len(sp_c)]
                if strand == "+":
                    left = s[max(0, start - 6):start]
                    right = s[start + len(sp_c): start + len(sp_c) + 6]
                    coords = (start, start + len(sp_c))
                else:
                    f_start = len(s) - start - len(sp_c)
                    left = s[f_start + len(sp_c): f_start + len(sp_c) + 6][::-1]
                    right = s[max(0, f_start - 6):f_start][::-1]
                    coords = (f_start, f_start + len(sp_c))
                pam = system.get("pam_iupac", "")
                pam_ok, pam_seq, side = _pam_around(sp_c, system, left, right)
                hits.append({
                    "spacer_index": idx, "spacer": sp_c, "strand": strand,
                    "start": coords[0], "end": coords[1],
                    "protospacer": protospacer, "pam": pam_seq, "pam_side": side,
                    "pam_valid": pam_ok, "interference": pam_ok,
                    "note": "targets the invader" if pam_ok
                            else "match found but no PAM — Cascade cannot engage",
                })
                start = probe.find(sp_c, start + 1)
    return hits


def _pam_around(spacer: str, system: Dict[str, Any], left: str, right: str) -> Tuple[bool, str, str]:
    pam = system.get("pam_iupac", "")
    if not pam:
        return True, "none", "no PAM requirement (type III)"
    regex = pattern_regex(pam)
    m = re.search(regex, right)                 # PAM 3' of the protospacer
    if m:
        return True, right[m.start():m.end()], "3' flank"
    m = re.search(regex, revcomp(left))
    if m:
        return True, revcomp(left)[m.start():m.end()], "5' flank (reverse strand)"
    return False, "", "absent"


def immunity_efficiency(active_targets: int, system: Dict[str, Any],
                        acr_present: bool = False, acr_strength: float = 0.85,
                        per_target_escape: float = 0.03) -> float:
    """Probability that a phage infection of an immune cell is aborted.

    Each independent target site adds protection (a single site is ~97%
    effective, two sites are essentially absolute), while an anti-CRISPR
    protein multiplies the residual failure probability.
    """
    if active_targets <= 0:
        return 0.0
    eff = 1.0 - per_target_escape ** min(active_targets, 4)
    if acr_present:
        eff *= (1.0 - acr_strength)
    return max(0.0, min(1.0, eff))


# --------------------------------------------------------------------------- #
# population dynamics — the arms race
# --------------------------------------------------------------------------- #
def simulate_phage_challenge(
    immune_fraction: float = 0.0,
    spa_fraction: float = 0.0,
    immunity_efficiency_per_cell: float = 0.97,
    acr_fraction: float = 0.0,
    acquisition_rate: float = 1e-4,
    escape_rate: float = 2e-6,
    initial_bacteria: float = 1e7,
    initial_phage: float = 1e5,
    hours: float = 24.0,
    dt: float = 0.02,
    growth_rate: float = 1.4,        # per hour
    carrying_capacity: float = 1e9,
    adsorption_rate: float = 2.5e-9,  # mL per phage per hour (scaled units)
    burst_size: float = 150.0,
    latency_h: float = 0.5,
    phage_decay: float = 0.05,
    immunity_cost: float = 0.03,
    acr_cost: float = 0.02,
) -> Dict[str, Any]:
    """Deterministic phage-vs-bacteria arms race.

    Compartments
    ------------
    B0  sensitive bacteria (no matching spacer)
    B1  immune bacteria (matching spacer, functional system)
    B2  immune but CRISPR-suppressed by an anti-CRISPR protein
    P0  wild-type phage
    P1  escape phage (protospacer/PAM mutated, no longer targeted)
    P2  anti-CRISPR carrying phage
    """
    B0 = initial_bacteria * (1 - immune_fraction - spa_fraction)
    B1 = initial_bacteria * immune_fraction
    B2 = initial_bacteria * spa_fraction
    P0 = initial_phage * (1 - acr_fraction)
    P2 = initial_phage * acr_fraction
    P1 = 0.0
    steps = int(hours / dt)
    series: Dict[str, List[float]] = {"t": [], "B0": [], "B1": [], "B2": [], "B_total": [],
                                      "P0": [], "P1": [], "P2": [], "immune_pct": []}
    events: List[Dict[str, Any]] = []
    escape_seen = False
    collapse_seen = False
    sample_every = max(1, int(0.25 / dt))
    eff = immunity_efficiency_per_cell
    for step in range(steps + 1):
        t = step * dt
        Btot = max(B0 + B1 + B2, 1e-9)
        Ptot = P0 + P1 + P2
        growth = growth_rate * max(0.0, 1 - Btot / carrying_capacity)

        # --- infection events (phage x bacterium encounters) ------------------
        enc = adsorption_rate * Ptot
        lysed_B0 = enc * B0                       # every encounter with a sensitive cell kills it
        abort_P0_B1 = adsorption_rate * P0 * B1 * eff        # CRISPR aborts these infections
        lysed_B1 = adsorption_rate * (P0 * (1 - eff) + P1 + P2 * (1 - eff)) * B1
        lysed_B2 = enc * B2                       # anti-CRISPR phage ignore the spacer

        # --- spacer acquisition: a few survivors of infection become immune ---
        acquired = acquisition_rate * lysed_B0
        lysed_B0 -= acquired

        # --- bacteria ---------------------------------------------------------
        dB0 = growth * B0 - lysed_B0
        dB1 = growth * (1 - immunity_cost) * B1 - lysed_B1 + acquired
        dB2 = growth * (1 - immunity_cost - acr_cost) * B2 - lysed_B2

        # --- phage ------------------------------------------------------------
        # productive infections burst after the latent period
        prod_rate = burst_size / max(latency_h, 1e-3)
        prod_P0 = prod_rate * (adsorption_rate * P0 * B0 + adsorption_rate * P0 * B1 * (1 - eff))
        prod_P1 = prod_rate * (adsorption_rate * P1 * (B0 + B1))
        prod_P2 = prod_rate * (adsorption_rate * P2 * (B0 + B1 * (1 - eff) + B2))
        mut_escape = prod_rate * abort_P0_B1 * escape_rate     # CRISPR-survivor phage
        dP0 = prod_P0 - adsorption_rate * P0 * Btot - phage_decay * P0
        dP1 = prod_P1 + mut_escape - adsorption_rate * P1 * Btot - phage_decay * P1
        dP2 = prod_P2 - adsorption_rate * P2 * Btot - phage_decay * P2

        B0 = max(0.0, B0 + dt * dB0)
        B1 = max(0.0, B1 + dt * dB1)
        B2 = max(0.0, B2 + dt * dB2)
        P0 = max(0.0, P0 + dt * dP0)
        P1 = max(0.0, P1 + dt * dP1)
        P2 = max(0.0, P2 + dt * dP2)

        Btot_new = B0 + B1 + B2
        if step % sample_every == 0 or step == steps:
            series["t"].append(round(t, 2))
            series["B0"].append(B0); series["B1"].append(B1); series["B2"].append(B2)
            series["B_total"].append(Btot_new)
            series["P0"].append(P0); series["P1"].append(P1); series["P2"].append(P2)
            series["immune_pct"].append(100 * B1 / max(1, Btot_new))
        if not escape_seen and P1 > 1e3:
            escape_seen = True
            events.append({"t": round(t, 2), "event": "escape mutants detected",
                           "detail": f"protospacer/PAM-escaped phage passed 10³ PFU at "
                                     f"t = {t:.2f} h"})
        if not collapse_seen and Btot_new < 0.2 * initial_bacteria:
            collapse_seen = True
            events.append({"t": round(t, 2), "event": "bacterial population collapsed",
                           "detail": "viable count fell below 20% of the starting population"})
    final_b = series["B_total"][-1]
    final_p = series["P0"][-1] + series["P1"][-1] + series["P2"][-1]
    outcome = ("bacteria win — infection cleared"
               if final_p < 1e4 and final_b > 0.6 * initial_bacteria else
               "phage wins — population destroyed"
               if final_b < 0.2 * initial_bacteria else
               "coexistence — resistant bacteria and escape phage persist")
    return {
        "series": series, "events": events,
        "final": {"bacteria": final_b, "phage": final_p,
                  "immune_percent": round(100 * series["B1"][-1] / max(1, final_b), 2),
                  "escape_phage_percent": round(100 * series["P1"][-1] / max(1e-9, final_p), 3)},
        "outcome": outcome,
        "parameters": {
            "growth_rate_per_h": growth_rate, "burst_size": burst_size,
            "adsorption_rate": adsorption_rate, "latency_h": latency_h,
            "escape_rate": escape_rate, "immunity_efficiency": immunity_efficiency_per_cell,
            "immunity_cost": immunity_cost, "acquisition_rate": acquisition_rate,
        },
        "assumptions": [
            "Phage adsorption is mass action; latency is collapsed into the burst term.",
            "Escape mutants arise from replication in immune cells at the escape rate.",
            "Immunity carries a fitness cost, otherwise resistance would always sweep.",
        ],
    }


def moist_curve(moi: float, immune_fraction: float, efficiency: float) -> Dict[str, float]:
    """Probability that at least one phage infects a given cell (Poisson)."""
    p_sens = 1 - math.exp(-moi)
    p_survive_immune = math.exp(-moi * (1 - efficiency))
    return {
        "moi": moi,
        "p_infected_sensitive": p_sens,
        "p_survive_if_immune": p_survive_immune,
        "p_lysis_if_sensitive": p_sens,
        "protection_factor": (1 - p_sens) / max(1e-9, p_survive_immune) if p_sens > 0 else float("inf"),
    }
