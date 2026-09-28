"""DNA repair, editing outcome and safety modelling.

The models here are *mechanistically structured* rather than fully
first-principles: pathway competition follows the cell-cycle logic of
NHEJ/MMEJ/HDR, indel spectra follow the measured shape of Cas9 editing
(1-2 bp indels dominate, deletions outnumber insertions), microhomology-driven
deletions are computed from the actual sequence, and base/prime editor
outcomes are parameterised from published reporter loci.  Every assumption is
labelled in the returned dictionaries so the UI can show its work.
"""
from __future__ import annotations

import math
import random
from typing import Any, Dict, List

from .sequtils import clean, gc_content, microhomologies, revcomp, translate


# --------------------------------------------------------------------------- #
# biological context models
# --------------------------------------------------------------------------- #
CELL_TYPES: List[Dict[str, Any]] = [
    {
        "id": "hek293t", "name": "HEK293T (embryonic kidney)",
        "hdr_factor": 1.00, "nhej_factor": 1.00, "mmej_factor": 1.00,
        "s_g2_fraction": 0.50, "doubling_time_h": 24, "ploidy": "2n (hypertriploid in practice)",
        "p53_competent": False, "indel_spectrum": "1 bp >> 2-10 bp; small insertions common",
        "notes": "The reference line for reported editing efficiencies: transformed, "
                 "defective DNA-damage checkpoints, therefore unusually permissive.",
    },
    {
        "id": "k562", "name": "K562 (chronic myeloid leukaemia)",
        "hdr_factor": 0.85, "nhej_factor": 1.05, "mmej_factor": 1.10,
        "s_g2_fraction": 0.55, "doubling_time_h": 20, "ploidy": "3n (near-triploid)",
        "p53_competent": False, "indel_spectrum": "multi-allelic; large deletions more visible",
        "notes": "Suspension line with high MMEJ activity; three alleles complicate "
                 "genotyping.",
    },
    {
        "id": "ipsc", "name": "Human iPSC",
        "hdr_factor": 0.55, "nhej_factor": 0.95, "mmej_factor": 1.05,
        "s_g2_fraction": 0.45, "doubling_time_h": 20, "ploidy": "2n",
        "p53_competent": True, "indel_spectrum": "clonal selection biases toward smaller indels",
        "notes": "Low HDR unless synchronised; single-cell cloning after editing adds "
                 "selection bias and takes weeks.",
    },
    {
        "id": "primary_t", "name": "Primary human T cell (CAR-T workflow)",
        "hdr_factor": 0.45, "nhej_factor": 1.20, "mmej_factor": 1.20,
        "s_g2_fraction": 0.35, "doubling_time_h": 30, "ploidy": "2n",
        "p53_competent": True, "indel_spectrum": "efficient indels (TRAC/PDCD1); large "
                                                 "deletions measurable",
        "notes": "Activated T cells edit well as RNPs, but knock-in by HDR stays low; "
                 "AAV6 donor helps a little.",
    },
    {
        "id": "cd34_hspc", "name": "CD34+ haematopoietic stem cell",
        "hdr_factor": 0.40, "nhej_factor": 1.15, "mmej_factor": 1.15,
        "s_g2_fraction": 0.40, "doubling_time_h": 36, "ploidy": "2n",
        "p53_competent": True, "indel_spectrum": "indels plus variable large deletions",
        "notes": "The cell type used by Casgevy: ex vivo editing of CD34+ cells, "
                 "myeloablative conditioning, autologous reinfusion.",
    },
    {
        "id": "hepatocyte", "name": "Human primary hepatocyte / liver",
        "hdr_factor": 0.20, "nhej_factor": 1.25, "mmej_factor": 1.10,
        "s_g2_fraction": 0.10, "doubling_time_h": 500, "ploidy": "2n (polyploid subsets)",
        "p53_competent": True, "indel_spectrum": "indels dominate; NHEJ is the only route "
                                                 "to a durable knockout",
        "notes": "Largely post-mitotic — the reason in-vivo programmes (NTLA-2001, "
                 "VERVE-101) knock genes out instead of repairing them.",
    },
    {
        "id": "neuron", "name": "Post-mitotic neuron",
        "hdr_factor": 0.03, "nhej_factor": 0.70, "mmej_factor": 0.60,
        "s_g2_fraction": 0.0, "doubling_time_h": 1e6, "ploidy": "2n (often 4n)",
        "p53_competent": True, "indel_spectrum": "few DSBs repaired by c-NHEJ; "
                                                 "microhomology and alt-EJ dominate",
        "notes": "No S/G2 phase means no HDR at all; editing strategies in neurons are "
                 "restricted to knockout, base editing or prime editing.",
    },
    {
        "id": "bacterium", "name": "Bacterium (E. coli replicon)",
        "hdr_factor": 0.90, "nhej_factor": 0.0, "mmej_factor": 0.0,
        "s_g2_fraction": 1.0, "doubling_time_h": 0.33, "ploidy": "1n (monoploid)",
        "p53_competent": False, "indel_spectrum": "no NHEJ at all — a DSB is lethal unless a "
                                                  "homology template is supplied, in which case "
                                                  "recombination is very efficient",
        "notes": "Bacteria have no NHEJ pathway at all: a single cut is usually a "
                 "deAD(B) event unless an HR template is supplied. This is exactly why "
                 "CRISPR works so well for curing bacteria of plasmids.",
    },
]

DELIVERY: List[Dict[str, Any]] = [
    {"id": "rnp", "name": "RNP electroporation",
     "onset_h": 0.5, "duration_h": 48, "efficiency": 1.00, "toxicity": 0.25,
     "immunogenicity": 0.20, "payload_kb": 0.2,
     "notes": "Pre-formed Cas9-guide complex; the fastest, cleanest route, but no "
              "sustained expression for HDR donors."},
    {"id": "mrna_lipid", "name": "mRNA + sgRNA lipid nanoparticle (LNP)",
     "onset_h": 3, "duration_h": 72, "efficiency": 0.95, "toxicity": 0.35,
     "immunogenicity": 0.45, "payload_kb": 5.0,
     "notes": "The in-vivo delivery of choice (liver tropism); transient, so "
              "specificity is time-limited."},
    {"id": "plasmid", "name": "Plasmid transfection",
     "onset_h": 12, "duration_h": 168, "efficiency": 0.80, "toxicity": 0.40,
     "immunogenicity": 0.55, "payload_kb": 12.0,
     "notes": "Cheap and long-lasting — and exactly that persistence drives "
              "off-target editing and innate immune sensing."},
    {"id": "aav6", "name": "AAV6 donor + RNP",
     "onset_h": 6, "duration_h": 336, "efficiency": 0.90, "toxicity": 0.30,
     "immunogenicity": 0.70, "payload_kb": 4.7,
     "notes": "The classic HDR knock-in route in T cells and HSPCs; capsid "
              "immunity and the 4.7 kb packaging limit are the constraints."},
    {"id": "lenti", "name": "Lentiviral delivery",
     "onset_h": 24, "duration_h": 1e6, "efficiency": 0.85, "toxicity": 0.60,
     "immunogenicity": 0.75, "payload_kb": 8.0,
     "notes": "Stable integration of the editor — powerful for screens, unacceptable "
              "for most therapies."},
]

#: Asynchronous culture distribution (G1, S, G2, M).
PHASE_FRACTION = {"G1": 0.45, "S": 0.30, "G2": 0.20, "M": 0.05}


# --------------------------------------------------------------------------- #
# pathway choice
# --------------------------------------------------------------------------- #
def repair_pathways(cell: Dict[str, Any], phase: str = "G1",
                    template: bool = False, template_type: str = "none",
                    mmej_score: float = 0.0, cas_kinetics: float = 1.0,
                    syncing: float = 0.0) -> Dict[str, Any]:
    """Competition between c-NHEJ, MMEJ, HDR and SSA at one DSB.

    Returns fractional probabilities that sum to 1. ``mmej_score`` (0-1) is the
    microhomology availability computed from the real sequence around the cut;
    ``syncing`` (0-1) shifts the population toward S/G2 (e.g. nocodazole).
    """
    s_g2 = PHASE_FRACTION.get(phase, 0.45) if phase else 0.45
    if phase in ("S", "G2"):
        s_g2 = 1.0
    elif phase in ("G1", "M"):
        s_g2 = 0.0
    else:  # asynchronous
        cycle = cell.get("s_g2_fraction", 0.5)
        s_g2 = min(1.0, cycle + syncing * (1 - cycle))

    c_nhej = 0.82 * cell.get("nhej_factor", 1.0)
    mmej = (0.10 + 0.25 * mmej_score) * cell.get("mmej_factor", 1.0) * (0.6 + 0.6 * s_g2)
    hdr = 0.0
    if template:
        donor_factor = {"none": 0.0, "ssODN": 0.55, "dsDNA": 0.80,
                        "plasmid": 0.70, "aav": 0.95, "lentivirus": 0.90}.get(template_type, 0.4)
        hdr = 0.95 * donor_factor * cell.get("hdr_factor", 1.0) * s_g2
        c_nhej *= (1 - 0.35 * s_g2)      # HR factors compete with c-NHEJ in S/G2
    ssa = 0.02 * cell.get("mmej_factor", 1.0) if mmej_score > 0.5 else 0.0

    total = c_nhej + mmej + hdr + ssa
    probs = {
        "c-NHEJ": c_nhej / total,
        "MMEJ": mmej / total,
        "HDR": hdr / total,
        "SSA": ssa / total,
    }
    return {
        "probabilities": {k: round(v, 4) for k, v in probs.items()},
        "percent": {k: round(100 * v, 2) for k, v in probs.items()},
        "s_g2_fraction": round(s_g2, 3),
        "reasoning": _pathway_reasoning(cell, phase, template, template_type, mmej_score, s_g2),
    }


def _pathway_reasoning(cell: Dict[str, Any], phase: str, template: bool,
                       template_type: str, mmej_score: float, s_g2: float) -> List[str]:
    out = [
        f"{cell['name']}: S/G2 fraction used by the model = {s_g2:.0%}.",
        "c-NHEJ operates throughout the cycle and is the default fate of a DSB in human cells.",
    ]
    if phase in ("G1", "M"):
        out.append("The simulated cell is outside S/G2, so HDR is set to zero by definition.")
    if mmej_score > 0.2:
        out.append(f"Flanking microhomology gives MMEJ an available substrate (score {mmej_score:.2f}).")
    if template:
        out.append(f"Donor present ({template_type}) — HDR gains a template term.")
    else:
        out.append("No donor template: precise repair is impossible, only indels are produced.")
    if cell.get("nhej_factor", 1) == 0:
        out.append("This organism lacks NHEJ entirely — an unrepaired DSB is lethal.")
    return out


def microhomology_score(seq: str, cut: int, min_len: int = 3) -> Dict[str, Any]:
    """Quantify MMEJ substrate at a cut from the real sequence."""
    mh = microhomologies(seq, cut, min_len=min_len)
    longest = max((h[2] for h in mh), default=0)
    score = min(1.0, longest / 10.0) if longest else 0.0
    return {"microhomologies": mh, "longest": longest, "score": score,
            "sequence_context": seq[max(0, cut - 25):cut + 25]}


# --------------------------------------------------------------------------- #
# indel spectra
# --------------------------------------------------------------------------- #
def _indel_size_distribution(cell: Dict[str, Any], pathway: str, rng: random.Random,
                             n: int) -> List[int]:
    """Draw indel size deltas: negative = deletion, positive = insertion."""
    mmej_bias = 0.55 if pathway == "MMEJ" else 0.0
    out: List[int] = []
    for _ in range(n):
        u = rng.random()
        if u < 0.32 - 0.10 * mmej_bias:                       # -1
            out.append(-1)
        elif u < 0.50 - 0.08 * mmej_bias:                     # -2
            out.append(-2)
        elif u < 0.60:                                        # -3..-5
            out.append(-rng.randint(3, 5))
        elif u < 0.68:                                        # -6..-20
            out.append(-rng.randint(6, 20))
        elif u < 0.74:                                        # microhomology-sized deletions
            out.append(-rng.choice([7, 8, 10, 11, 13, 17, 21, 25]))
        elif u < 0.80:
            out.append(rng.randint(-120, -30))                # large deletion
        elif u < 0.93:                                        # +1
            out.append(1)
        elif u < 0.97:
            out.append(rng.randint(2, 5))
        else:
            out.append(rng.randint(6, 40))
    return out


def simulate_indels(seq: str, cut: int, n_alleles: int, cell: Dict[str, Any],
                    pathway: str = "c-NHEJ", seed: int = 0) -> Dict[str, Any]:
    """Generate an allele population of indel outcomes at a cut site."""
    rng = random.Random(seed)
    sizes = _indel_size_distribution(cell, pathway, rng, n_alleles)
    alleles: List[Dict[str, Any]] = []
    for size in sizes:
        if size < 0:                                     # deletion
            lo = max(0, cut + size)
            deleted = seq[lo:cut]
            allele = seq[:lo] + seq[cut:]
            kind = "deletion"
            detail = deleted
        else:                                            # insertion
            ins = "".join(rng.choice("ACGT") for _ in range(size))
            allele = seq[:cut] + ins + seq[cut:]
            kind = "insertion"
            detail = ins
        alleles.append({"size": size, "kind": kind, "detail": detail,
                        "sequence": allele, "in_frame": size % 3 == 0})
    counts: Dict[int, int] = {}
    for s in sizes:
        counts[s] = counts.get(s, 0) + 1
    frameshift = sum(1 for s in sizes if s % 3 != 0) / max(1, len(sizes))
    return {
        "alleles": alleles,
        "size_histogram": dict(sorted(counts.items())),
        "frameshift_fraction": frameshift,
        "mean_size": sum(sizes) / max(1, len(sizes)),
        "pct_1bp": 100 * sum(1 for s in sizes if abs(s) == 1) / max(1, len(sizes)),
        "pct_large_del": 100 * sum(1 for s in sizes if s <= -30) / max(1, len(sizes)),
    }


def mmej_products(seq: str, cut: int, min_len: int = 3) -> List[Dict[str, Any]]:
    """Deletions that microhomology actually predicts at this cut (sequence-specific)."""
    out: List[Dict[str, Any]] = []
    for left, right, length, repeat in microhomologies(seq, cut, min_len=min_len):
        if right <= left:
            continue
        deleted = seq[left:right]
        product = seq[:left] + seq[right:]
        out.append({
            "repeat": repeat, "repeat_len": length,
            "deletion_len": right - left,
            "deleted_sequence": deleted,
            "product": product,
            "left_flank": seq[max(0, left - 12):left],
            "right_flank": seq[right:right + 12],
        })
    return sorted(out, key=lambda d: -d["repeat_len"])[:6]


# --------------------------------------------------------------------------- #
# full editing simulation
# --------------------------------------------------------------------------- #
def simulate_nuclease_edit(seq: str, cut: int, cell: Dict[str, Any],
                           dose: float = 1.0, phase: str = "asynchronous",
                           template: bool = False, template_type: str = "none",
                           desired_edit: str = "", mmej_min: int = 3,
                           n_cells: int = 1000, seed: int = 42) -> Dict[str, Any]:
    """Simulate a full nuclease experiment and return allele + outcome fractions."""
    mh = microhomology_score(seq, cut, mmej_min)
    pathways = repair_pathways(cell, phase if phase != "asynchronous" else None,
                              template, template_type, mh["score"])
    rng = random.Random(seed)
    p = pathways["probabilities"]
    # delivery/dose scaling: fraction of cells that get a cut at all
    cut_fraction = min(0.98, 0.85 * dose) if cell["id"] != "bacterium" else min(0.99, 0.95 * dose)
    n_cut = int(n_cells * cut_fraction)
    outcomes = {"unedited": n_cells - n_cut, "c-NHEJ": 0, "MMEJ": 0, "HDR": 0, "SSA": 0}
    hdr_precise = 0
    hdr_imprecise = 0
    alleles: List[Dict[str, Any]] = []
    for _ in range(n_cut):
        u = rng.random()
        if u < p["c-NHEJ"]:
            outcomes["c-NHEJ"] += 1
        elif u < p["c-NHEJ"] + p["MMEJ"]:
            outcomes["MMEJ"] += 1
        elif u < p["c-NHEJ"] + p["MMEJ"] + p["HDR"]:
            outcomes["HDR"] += 1
            if rng.random() < 0.92:
                hdr_precise += 1
            else:
                hdr_imprecise += 1
        else:
            outcomes["SSA"] += 1
    indels = simulate_indels(seq, cut, outcomes["c-NHEJ"] + outcomes["MMEJ"], cell,
                             "MMEJ" if outcomes["MMEJ"] > outcomes["c-NHEJ"] else "c-NHEJ",
                             seed=seed)
    mmej_pred = mmej_products(seq, cut)
    indel_fraction = (outcomes["c-NHEJ"] + outcomes["MMEJ"] + outcomes["SSA"]) / n_cells
    ko_any = 1 - (1 - indel_fraction) ** 2 if cell["ploidy"].startswith("2") else indel_fraction
    # a knockout needs a frameshift (or a large deletion) in all alleles
    fs = indels["frameshift_fraction"]
    p_null_allele = indel_fraction * fs + 0.15 * (1 - fs) * indel_fraction
    if cell["ploidy"].startswith("2"):
        biallelic_ko = p_null_allele ** 2 + 2 * p_null_allele * (1 - p_null_allele) * 0.0
    else:
        biallelic_ko = p_null_allele
    return {
        "cell": cell["name"],
        "pathways": pathways,
        "counts": outcomes,
        "outcome_percent": {k: round(100 * v / n_cells, 2) for k, v in outcomes.items()},
        "indels": indels,
        "mmej_predictions": mmej_pred,
        "microhomology": {"longest": mh["longest"], "score": mh["score"]},
        "indel_fraction": indel_fraction,
        "hdr_fraction": outcomes["HDR"] / n_cells,
        "hdr_precise_fraction": hdr_precise / n_cells,
        "hdr_imprecise_fraction": hdr_imprecise / n_cells,
        "knockout_allele_fraction": p_null_allele,
        "biallelic_knockout_fraction": biallelic_ko,
        "alleles": alleles,
        "assumptions": [
            "Pathway competition is cell-cycle gated: HDR is only available in S/G2.",
            "Indel sizes are drawn from an empirical spectrum, not from a repair biophysics model.",
            "Microhomology-driven deletions are computed from the actual flanking sequence.",
        ],
    }


# --------------------------------------------------------------------------- #
# base editing
# --------------------------------------------------------------------------- #
#: context preference multipliers (rAPOBEC1 prefers TC, then CC, AC, GC)
CBE_CONTEXT = {"TC": 1.00, "CC": 0.55, "AC": 0.35, "GC": 0.18}
ABE_CONTEXT_5 = {"A": 0.55, "C": 1.00, "G": 0.95, "T": 1.05}   # preceding base
ABE_CONTEXT_3 = {"A": 0.85, "C": 1.10, "G": 0.95, "T": 0.95}   # following base


def base_editing_window(site: Dict[str, Any], editor: Dict[str, Any]) -> List[Dict[str, Any]]:
    """List editable bases inside the editor window for one target site."""
    ps = site["protospacer"]
    lo, hi = editor["window"]
    want = "C" if editor["family"].startswith("CBE") else ("A" if editor["family"].startswith("ABE")
                                                           else "A")
    out: List[Dict[str, Any]] = []
    for i, base in enumerate(ps, start=1):
        if base != want:
            continue
        if not (lo <= i <= hi):
            continue
        left = ps[i - 2] if i >= 2 else "N"
        right = ps[i] if i < len(ps) else "N"
        if editor["family"].startswith("CBE"):
            w = CBE_CONTEXT.get(left + base, 0.2)
        else:
            w = ABE_CONTEXT_5.get(left, 1.0) * ABE_CONTEXT_3.get(right, 1.0)
        # window edge taper
        centre = (lo + hi) / 2
        taper = max(0.15, 1.0 - 0.35 * abs(i - centre))
        out.append({
            "position": i, "base": base, "context": f"{left}{base}{right}",
            "relative_efficiency": round(w * taper, 3),
            "is_bystander": False,
        })
    if out:
        best = max(o["relative_efficiency"] for o in out)
        for o in out:
            o["relative_efficiency"] = round(o["relative_efficiency"] / best, 3)
            o["is_bystander"] = o["relative_efficiency"] < 0.99
    return out


def simulate_base_editing(site: Dict[str, Any], editor: Dict[str, Any], cell: Dict[str, Any],
                          dose: float = 1.0, n_cells: int = 1000, seed: int = 7) -> Dict[str, Any]:
    """Simulate CBE/ABE population outcomes, including bystander editing."""
    rng = random.Random(seed)
    editable = base_editing_window(site, editor)
    if not editable:
        return {"editable": [], "products": {}, "editing_efficiency": 0.0,
                "note": "No editable base of the right identity inside this editor's window.",
                "allele_composition": {}}
    peak = min(0.95, editor["max_efficiency"] * dose * cell.get("hdr_factor", 1.0) ** 0
               * (0.75 + 0.25 * dose))
    products: Dict[str, int] = {}
    conversions = {e["position"]: 0 for e in editable}
    for _ in range(n_cells):
        pattern = []
        for e in editable:
            p = max(0.02, peak * e["relative_efficiency"])
            hit = rng.random() < p
            if hit:
                conversions[e["position"]] += 1
            pattern.append("1" if hit else "0")
        key = "".join(pattern)
        products[key] = products.get(key, 0) + 1
    # indels that accompany base editing
    n_indel = sum(1 for _ in range(n_cells) if rng.random() < editor["indel_rate"] * dose)
    edited_fraction = sum(v for k, v in products.items() if k != "0" * len(editable)) / n_cells
    purity = (products.get("0" * len(editable), 0)
              + sum(v for k, v in products.items() if k.count("1") == 1)) / n_cells
    allele_map: Dict[str, int] = {}
    for key, count in products.items():
        if key == "0" * len(editable):
            allele_map["unedited"] = count
        elif key.count("1") == 1:
            allele_map["precise single conversion"] = allele_map.get("precise single conversion", 0) + count
        else:
            allele_map["multiply converted (bystander)"] = allele_map.get("multiply converted (bystander)", 0) + count
    if n_indel:
        allele_map["indel"] = n_indel
    return {
        "editable": editable,
        "products": dict(sorted(products.items(), key=lambda kv: -kv[1])[:12]),
        "editing_efficiency": round(100 * edited_fraction, 2),
        "product_purity": round(100 * purity, 2),
        "per_position_conversion": {p: round(100 * c / n_cells, 2) for p, c in conversions.items()},
        "indel_rate_percent": round(100 * n_indel / n_cells, 2),
        "allele_composition": allele_map,
        "conversion": editor["conversion"],
        "transversion_risk": "CBE can also produce C>G / C>A transversions (~1-3% of alleles)."
        if editor["family"].startswith("CBE") else "ABE byproducts are mostly A>G; C>T side-products are rare.",
    }


# --------------------------------------------------------------------------- #
# prime editing
# --------------------------------------------------------------------------- #
def prime_editing_efficiency(pbs_tm: float, pbs_gc: float, rtt_len: int,
                             editor: Dict[str, Any], cell: Dict[str, Any],
                             pe3: bool = False, motif: bool = True) -> Dict[str, Any]:
    """Predict PE efficiency from pegRNA design parameters (heuristic, documented)."""
    tm_term = math.exp(-((pbs_tm - 45.0) ** 2) / (2 * 12.0 ** 2))         # best around 45 °C
    gc_term = math.exp(-((pbs_gc - 0.50) ** 2) / (2 * 0.22 ** 2))
    rtt_term = 1.0 if 12 <= rtt_len <= 34 else max(0.25, 1 - abs(rtt_len - 20) / 40)
    motif_term = 1.12 if motif else 1.0
    base = editor.get("max_efficiency", 0.2) * cell.get("hdr_factor", 1.0)
    eff = base * (0.45 + 0.55 * tm_term) * (0.6 + 0.4 * gc_term) * rtt_term * motif_term
    if pe3:
        eff *= 1.7
    eff = min(eff, 0.85)
    indels = editor.get("indel_rate", 0.02) * (3.0 if pe3 else 1.0)
    return {
        "efficiency": eff,
        "efficiency_percent": round(100 * eff, 2),
        "indel_rate": round(100 * indels, 2),
        "precise_edit_percent": round(100 * eff * (1 - indels), 2),
        "terms": {"pbs_tm": round(tm_term, 3), "pbs_gc": round(gc_term, 3),
                  "rtt_length": round(rtt_term, 3), "motif": motif_term},
        "notes": ["Flap resolution is stochastic: the same pegRNA gives precise edits, "
                  "indels and unedited alleles in one population.",
                  "PE3 roughly doubles efficiency and triples indels — that is the trade.",
                  "The 3' extension motif (evopreQ1) protects the pegRNA from exonuclease."],
    }


def simulate_prime_editing(seq: str, cut: int, edit: str, peg: Dict[str, Any],
                           editor: Dict[str, Any], cell: Dict[str, Any],
                           pe3: bool = False, n_cells: int = 1000, seed: int = 11) -> Dict[str, Any]:
    """Allele-level prime editing simulation for a point edit or small indel."""
    eff = prime_editing_efficiency(peg["pbs_tm"], peg["pbs_gc"], peg["rtt_len"],
                                   editor, cell, pe3, bool(peg.get("motif")))
    rng = random.Random(seed)
    precise = int(n_cells * eff["efficiency"] * (1 - eff["indel_rate"]))
    indel = int(n_cells * eff["efficiency"] * eff["indel_rate"])
    unedited = n_cells - precise - indel
    edited_seq = seq[:cut] + clean(edit) + seq[cut:]
    return {
        "summary": {"precise_edit": precise, "indel": indel, "unedited": unedited},
        "percent": {"precise_edit": round(100 * precise / n_cells, 2),
                    "indel": round(100 * indel / n_cells, 2),
                    "unedited": round(100 * unedited / n_cells, 2)},
        "edited_sequence": edited_seq,
        "original_sequence": seq,
        "model": eff,
    }


# --------------------------------------------------------------------------- #
# safety / burden
# --------------------------------------------------------------------------- #
def safety_panel(cell: Dict[str, Any], n_dsbs: float, off_target_hits: int,
                 on_target_indel_fraction: float, dose: float = 1.0,
                 hours_expressed: float = 48.0) -> Dict[str, Any]:
    """Estimate the collateral damage profile of an editing experiment.

    * off-target indel burden scales with the number of expressed editor-hours
    * large deletions and translocations scale with the number of DSBs
    * the p53 response only exists in checkpoint-competent cells
    """
    exposure = (hours_expressed / 48.0) * dose
    off_target_indels = min(0.95, 0.012 * off_target_hits * exposure)
    large_deletions = min(0.60, 0.008 * n_dsbs * exposure)
    translocations = min(0.30, 0.0016 * n_dsbs ** 1.6 * exposure) if n_dsbs > 1 else 0.0
    chromothripsis = min(0.05, 0.0004 * n_dsbs ** 1.8 * exposure)
    p53 = 0.0
    if cell.get("p53_competent"):
        p53 = min(0.95, 0.20 * n_dsbs * exposure)
    senescence = min(0.70, p53 * 0.5)
    return {
        "off_target_indel_fraction": round(off_target_indels, 4),
        "large_deletion_fraction": round(large_deletions, 4),
        "translocation_fraction": round(translocations, 4),
        "chromothripsis_fraction": round(chromothripsis, 5),
        "p53_response": round(p53, 4),
        "senescence_or_apoptosis": round(senescence, 4),
        "on_target_indel_fraction": round(on_target_indel_fraction, 4),
        "therapeutic_index": round(on_target_indel_fraction / max(1e-6, off_target_indels + 1e-6), 1),
        "notes": [
            "Transient editors (RNP, mRNA/LNP) reduce every number here — exposure is the lever.",
            "Two simultaneous cuts within one cell are when rearrangements become real.",
            "p53 activation is a genuine safety consideration in clinical editing.",
        ],
    }


def allele_frequency_curve(initial: float, generations: int = 20,
                           selection: float = 0.0, drift: float = 0.0,
                           seed: int = 5) -> List[Dict[str, float]]:
    """Simple allele fixation / loss trajectory (used for population panels)."""
    rng = random.Random(seed)
    f = initial
    out: List[Dict[str, float]] = []
    for gen in range(generations):
        f = f + selection * f * (1 - f)
        if drift:
            f += rng.gauss(0, drift) * math.sqrt(max(f * (1 - f), 1e-6))
        f = min(1.0, max(0.0, f))
        out.append({"generation": gen, "frequency": f})
    return out
