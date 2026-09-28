"""Clinical programmes, disease response models and reference data.

Numbers come from the cited publications, the FDA labels and the trial press
releases listed in ``data/references.md``.  Where a value is a modelled
projection rather than a measured one, the response functions below say so and
the UI marks the output as a projection.
"""
from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Tuple

TRIALS: List[Dict[str, Any]] = [
    {
        "id": "casgevy", "product": "Casgevy (exagamglogene autotemcel, exa-cel)",
        "sponsor": "CRISPR Therapeutics / Vertex", "modality": "CRISPR-Cas9 nuclease",
        "target": "BCL11A erythroid enhancer", "edit": "disruption of a GATA1 site -> HbF de-repression",
        "delivery": "ex vivo autologous CD34+ HSPCs, myeloablative busulfan conditioning",
        "indication": "Sickle-cell disease (VOC) and transfusion-dependent beta-thalassaemia",
        "phase": "Approved (FDA Dec 2023 for SCD, Jan 2024 for TDT; first CRISPR medicine)",
        "n_patients": "44 (SCD, CLIMB-121) / 52 (TDT, CLIMB-111)",
        "dose": "single infusion, minimum 3 x 10^6 CD34+ cells/kg",
        "outcomes": "SCD: 97% free of severe vaso-occlusive crises for >= 12 months. "
                    "TDT: 91% transfusion-independent for >= 12 months with mean Hb >= 9 g/dL.",
        "notes": "Editing is confined to the patient's own stem cells; no donor needed. "
                 "The HbF rise is the pharmacodynamics read-out.",
        "ref": "Frangoul 2021 NEJM 384:252; FDA STN 125787 (Dec 2023); FDA STN 125788 (Jan 2024)",
    },
    {
        "id": "ntla2001", "product": "NTLA-2001 (nexiguran ziclumeran / nex-z)",
        "sponsor": "Intellia / Regeneron", "modality": "CRISPR-Cas9 nuclease, in vivo",
        "target": "TTR (transthyretin)", "edit": "gene knockout in hepatocytes",
        "delivery": "single intravenous lipid nanoparticle carrying mRNA + guide",
        "indication": "ATTR amyloidosis with cardiomyopathy and polyneuropathy",
        "phase": "Phase 3",
        "n_patients": "> 60 patients dosed in the phase 1 study by 2023",
        "dose": "single ascending dose (0.1, 0.3, 0.7, 1.0 mg/kg in the first cohorts)",
        "outcomes": "Deep and durable serum TTR reduction after a single infusion: "
                    "up to ~90% median reduction at day 28 at the higher dose levels; "
                    "sustained beyond 12 months in the first cohort.",
        "notes": "The first systemically delivered CRISPR therapy in humans, and a proof "
                 "that knockout (not HDR) is the tractable in-vivo edit.",
        "ref": "Gillmore 2021 NEJM 385:493; Intellia press releases (2023-2024)",
    },
    {
        "id": "verve101", "product": "VERVE-101 (and successor VERVE-102)",
        "sponsor": "Verve Therapeutics", "modality": "Base editing (CBE), in vivo",
        "target": "PCSK9", "edit": "splice-site disruption -> permanent PCSK9 inactivation",
        "delivery": "single intravenous lipid nanoparticle",
        "indication": "Heterozygous familial hypercholesterolaemia with ASCVD",
        "phase": "Phase 1b (heart-1); superseded by VERVE-102",
        "n_patients": "10 patients across 4 dose cohorts (0.1-0.6 mg/kg)",
        "dose": "0.1, 0.3, 0.45, 0.6 mg/kg single infusion",
        "outcomes": "PCSK9 reduced by 47-84%; time-averaged LDL-C reduced by 39-55% "
                    "(up to 6 months at the highest dose). Transient ALT rises at higher "
                    "doses; two cardiovascular serious adverse events were reported.",
        "notes": "First in-vivo base-editing trial in humans. Non-human primate data showed "
                 "46-70% liver PCSK9 editing and 69% LDL-C reduction durable to 476 days.",
        "ref": "heart-1 NCT05398029; Kathiresan 2022 Circulation 147:1002; "
               "Han 2024 review; Verve AHA presentation 2023",
    },
    {
        "id": "edit101", "product": "EDIT-101",
        "sponsor": "Editas Medicine", "modality": "CRISPR-Cas9 nuclease (AAV5)",
        "target": "CEP290 IVS26", "edit": "deletion of the cryptic splice site",
        "delivery": "subretinal AAV5 injection",
        "indication": "Leber congenital amaurosis 10 (LCA10)",
        "phase": "Phase 1/2 (BRILLIANCE); programme deprioritised after readout",
        "n_patients": "14 (12 adults, 2 children) at low, intermediate and high dose",
        "dose": "single subretinal injection",
        "outcomes": "No dose-limiting toxicity; a subset of patients improved on at least "
                    "one of several vision tests (full-field light sensitivity, mobility, "
                    "visual acuity).",
        "notes": "The first in-vivo CRISPR trial for a genetic blindness. The tough part is "
                 "efficiency: the target is deep intronic and the tissue cannot be re-dosed freely.",
        "ref": "Maeder 2019 Nat Med 25:229; BRILLIANCE NCT03872479",
    },
    {
        "id": "ctxtcells", "product": "CTX110 / CTX130 / CRISPR-edited CAR-T programmes",
        "sponsor": "CRISPR Therapeutics, Caribou, Allogene and academic groups",
        "modality": "multiplex CRISPR-Cas9 editing",
        "target": "TRAC (to remove the endogenous TCR), PDCD1, B2M, CD70, CD19",
        "edit": "knockouts that make allogeneic T cells usable and exhaustion-resistant",
        "delivery": "ex vivo electroporation of RNPs into donor or patient T cells",
        "indication": "B-cell malignancies, renal cell carcinoma, solid tumours",
        "phase": "Phase 1/2",
        "n_patients": "dozens across programmes (CTX110: objective responses in relapsed "
                      "large B-cell lymphoma)",
        "dose": "single infusion after lymphodepletion",
        "outcomes": "Feasibility established: multiplex knockouts of TRAC and PDCD1 are "
                    "tolerated; allogeneic CAR-T avoids autologous manufacturing delay.",
        "notes": "Where editing is used as a manufacturing tool rather than as a cure: "
                 "three simultaneous edits per cell, and chromosome rearrangements between "
                 "cut sites are the key safety question.",
        "ref": "Lu 2020 Nature 579:262; Stadtmauer 2020 Science 367:eaba7365; CTX110 NCT04035434",
    },
    {
        "id": "ccr5trial", "product": "CCR5-edited CD4+ T cells (DdPc)",
        "sponsor": "Affiliated Hospital of Sun Yat-sen University",
        "modality": "CRISPR-Cas9 nuclease",
        "target": "CCR5", "edit": "knockout of the HIV-1 co-receptor",
        "delivery": "ex vivo electroporated autologous CD4+ T cells, then re-infusion",
        "indication": "HIV-1 infection",
        "phase": "Phase 1 (NCT03164135)",
        "n_patients": "1 patient",
        "dose": "single infusion",
        "outcomes": "Engraftment and long-term persistence of edited cells; HIV-1 viral load "
                    "rebounded after antiretroviral treatment interruption.",
        "notes": "Historic and ethically contested (He Jiankui's germline experiment used the "
                 "same target). Somatic CCR5 editing in consenting adults is legitimate; "
                 "germline editing is prohibited.",
        "ref": "Xu 2019 NEJM 381:1240; Tebas 2014 NEJM 370:901",
    },
    {
        "id": "beam101", "product": "BEAM-101 / BEAM-302",
        "sponsor": "Beam Therapeutics", "modality": "Base editing (ABE)",
        "target": "HBG1/HBG2 promoter HPFH sites; SERPINA1 (AATD)",
        "edit": "A>G at the -198/-175 HPFH sites; SERPINA1 Z allele correction",
        "delivery": "ex vivo CD34+ cells (BEAM-101); LNP in vivo (BEAM-302)",
        "indication": "Sickle-cell disease; alpha-1 antitrypsin deficiency",
        "phase": "Phase 1/2",
        "n_patients": "first patients dosed 2024-2025",
        "dose": "single dose",
        "outcomes": "Reported increases in HbF with no double-strand breaks, and correction "
                    "of the SERPINA1 Z allele in vivo in the first AATD patient.",
        "notes": "Base editing avoids translocations and large deletions because there is no "
                 "DSB — the trade-off is being restricted to transitions (A>G / C>T).",
        "ref": "Gaudelli 2017 Nature 551:464; Beam press releases 2024-2025",
    },
]

# --------------------------------------------------------------------------- #
# disease response models
# --------------------------------------------------------------------------- #
DISEASES: Dict[str, Dict[str, Any]] = {
    "scd": {
        "name": "Sickle-cell disease (HbSS)",
        "market": "BCL11A enhancer disruption (Casgevy) or HBG promoter base editing",
        "baseline": {"HbF_percent": 6.0, "HbS_percent": 80.0, "Hb_g_dL": 8.0,
                     "VOC_per_year": 3.5},
        "response": "Each percent of HbF (pan-cellular, F-cells) reduces polymerisation, "
                    "haemolysis and VOC frequency; clinically meaningful benefit starts "
                    "around 20-30% HbF.",
        "endpoint": "VOC-free for 12 months",
        "ref": "Frangoul 2021 NEJM 384:252; Platt 1994 NEJM 330:1639 (HbF threshold)",
    },
    "tbt": {
        "name": "Transfusion-dependent beta-thalassaemia",
        "market": "BCL11A enhancer disruption (Casgevy)",
        "baseline": {"HbF_percent": 10.0, "HbS_percent": 0.0, "Hb_g_dL": 7.5,
                     "transfusions_per_year": 20.0},
        "response": "Effective erythropoiesis returns once HbF replaces the missing beta-globin; "
                    "transfusion independence is the regulatory endpoint.",
        "endpoint": "Transfusion independence for 12 months",
        "ref": "Frangoul 2021 NEJM 384:252; Locatelli 2022 NEJM 386:1703",
    },
    "attr": {
        "name": "ATTR amyloidosis",
        "market": "TTR knockout (NTLA-2001)",
        "baseline": {"TTR_mg_dL": 25.0, "walk_test_m": 350.0},
        "response": "Serum TTR reduction tracks with clinical benefit in polyneuropathy; "
                    "a ~90% knockdown is the target of the in-vivo programmes.",
        "endpoint": "Reduction in serum TTR and slowed neurological progression",
        "ref": "Gillmore 2021 NEJM 385:493; Adams 2018 NEJM 379:11 (patisiran benchmark)",
    },
    "fh": {
        "name": "Heterozygous familial hypercholesterolaemia",
        "market": "PCSK9 inactivation by base editing (VERVE-101)",
        "baseline": {"LDL_mg_dL": 193.0, "PCSK9_percent": 100.0},
        "response": "A one-off edit permanently lowers PCSK9 and therefore LDL-C; a 50% "
                    "LDL-C reduction is similar to what antibody inhibitors achieve, but durable.",
        "endpoint": "Time-averaged LDL-C reduction",
        "ref": "Kathiresan 2022 Circulation 147:1002; heart-1 NCT05398029",
    },
    "lca10": {
        "name": "LCA10 (CEP290)",
        "market": "CEP290 IVS26 deletion (EDIT-101)",
        "baseline": {"vision_score": 0.15, "retinal_thickness_um": 120.0},
        "response": "Partial splice correction in a fraction of photoreceptors can rescue "
                    "measurable light sensitivity — the ceiling is set by how many cells are edited.",
        "endpoint": "Improvement in at least one of: FST, mobility, visual acuity",
        "ref": "Maeder 2019 Nat Med 25:229; BRILLIANCE NCT03872479",
    },
}


def editing_to_hbf(editing_percent: float, baseline_hbf: float = 6.0,
                   ceiling: float = 60.0) -> float:
    """Map BCL11A/HBG editing to pan-cellular HbF (saturating, from clinical data)."""
    e = max(0.0, editing_percent) / 100.0
    return baseline_hbf + (ceiling - baseline_hbf) * (1 - math.exp(-3.0 * e))


def voc_rate(hbf_percent: float, baseline_voc: float = 3.5) -> float:
    """Vaso-occlusive crisis rate as a function of HbF (steep above ~25% HbF)."""
    if hbf_percent >= 30:
        return baseline_voc * 0.05
    if hbf_percent <= 6:
        return baseline_voc
    frac = (hbf_percent - 6.0) / 24.0
    return baseline_voc * (1 - frac ** 1.4) * 0.9 + baseline_voc * 0.05 * frac


def ttr_knockdown(editing_percent: float, baseline: float = 25.0) -> float:
    return baseline * (1 - min(0.97, editing_percent / 100.0) * 1.0)


def ldl_response(editing_percent: float, baseline_ldl: float = 193.0) -> Tuple[float, float]:
    """Return (LDL-C, PCSK9) after editing; ~1:1 relationship at partial editing."""
    pcs = max(5.0, 100.0 * (1 - min(0.95, editing_percent / 100.0)))
    ldl = baseline_ldl * (0.35 + 0.65 * pcs / 100.0)
    return ldl, pcs


def simulate_cohort(disease: str, editing_percent: float, n: int = 40,
                    variability: float = 0.28, seed: int = 13) -> Dict[str, Any]:
    """Project a patient cohort: editing efficiency varies between patients."""
    rng = random.Random(seed)
    d = DISEASES[disease]
    rows: List[Dict[str, float]] = []
    for i in range(n):
        e = max(0.0, min(100.0, rng.gauss(editing_percent, editing_percent * variability)))
        if disease in ("scd", "tbt"):
            hbf = editing_to_hbf(e, d["baseline"]["HbF_percent"])
            voc = voc_rate(hbf, d["baseline"]["VOC_per_year"]) if disease == "scd" else 0.0
            rows.append({"editing_percent": e, "hbF_percent": hbf, "voc_per_year": voc,
                         "hemoglobin_g_dL": min(14.5, d["baseline"]["Hb_g_dL"] + 0.55 * (hbf - d["baseline"]["HbF_percent"]) / 3.0),
                         "responder": 1.0 if (voc <= 0.2 if disease == "scd" else hbf >= 25) else 0.0})
        elif disease == "attr":
            ttr = ttr_knockdown(e, d["baseline"]["TTR_mg_dL"])
            rows.append({"editing_percent": e, "ttr_mg_dL": ttr,
                         "walk_test_m": d["baseline"]["walk_test_m"] + 90 * (1 - ttr / d["baseline"]["TTR_mg_dL"]),
                         "responder": 1.0 if ttr < 0.25 * d["baseline"]["TTR_mg_dL"] else 0.0})
        elif disease == "fh":
            ldl, pcs = ldl_response(e, d["baseline"]["LDL_mg_dL"])
            rows.append({"editing_percent": e, "LDL_mg_dL": ldl, "PCSK9_percent": pcs,
                         "responder": 1.0 if ldl < 0.7 * d["baseline"]["LDL_mg_dL"] else 0.0})
        else:  # lca10
            rescue = min(0.75, e / 100.0 * 1.4)
            rows.append({"editing_percent": e, "vision_score": d["baseline"]["vision_score"] * (1 + 3.5 * rescue),
                         "responder": 1.0 if rescue > 0.15 else 0.0})
    responders = sum(r["responder"] for r in rows)
    return {
        "disease": d["name"], "rows": rows, "n": n,
        "response_rate_percent": round(100 * responders / n, 1),
        "mean_editing_percent": round(sum(r["editing_percent"] for r in rows) / n, 1),
        "model": "Between-patient variability modelled as a 28% coefficient of variation on "
                 "editing efficiency; response thresholds use the trial endpoints.",
        "disclaimer": "Projection, not a prediction: real cohorts are smaller and far more "
                      "heterogeneous than this model.",
    }
