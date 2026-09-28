"""Repair lab — the decision tree behind every editing outcome."""
from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd
import streamlit as st

from .. import theme, viz
from ..data.targets import all_loci, coding_strand, get_locus
from ..repair import CELL_TYPES, microhomology_score, repair_pathways, simulate_indels

CELL_BY_ID = {c["id"]: c for c in CELL_TYPES}

DECISION_TREE = f"""
<div style="display:flex;gap:12px;flex-wrap:wrap">
  <div class="card red" style="flex:1 1 210px">
    <h4>① The break</h4>
    <p>A blunt or staggered DSB. The <b>geometry</b> already biases the outcome: blunt ends
    (Cas9) are Ku70/80's favourite substrate, staggered ends (Cas12a) leave a short
    overhang that must be processed first.</p>
  </div>
  <div class="card" style="flex:1 1 210px">
    <h4>② End recognition</h4>
    <p><b>Ku70/80</b> binds first and blocks resection → c-NHEJ.<br>
    <b>MRN/CtIP</b> competes: it resects the 5' ends and commits the break to
    homology-directed routes.</p>
  </div>
  <div class="card teal" style="flex:1 1 210px">
    <h4>③ Resection + cycle</h4>
    <p>Short resection → <b>MMEJ</b> (anneals microhomologies, always deleting).<br>
    Long resection in <b>S/G2</b> → <b>HDR</b> with a sister chromatid or a supplied donor.<br>
    Long resection with no template → <b>SSA</b>, which deletes everything between the repeats.</p>
  </div>
  <div class="card violet" style="flex:1 1 210px">
    <h4>④ Outcome</h4>
    <p>c-NHEJ → small indels (mostly 1-2 bp).<br>
    MMEJ → predictable microhomology-bounded deletions.<br>
    HDR → precise, templated change (what everyone wants and rarely gets).<br>
    SSA → large deletions between repeats.</p>
  </div>
</div>
"""

BIOLOGY_CARDS = [
    ("c-NHEJ — non-homologous end joining",
     "Ku70/80 → DNA-PKcs → XRCC4/LIG4. Active in every phase of the cycle, error-prone at the "
     "junction, and the reason a nuclease knockout works at all. Yields mostly 1-2 bp indels; "
     "the reading frame is destroyed ~2/3 of the time.",
     "plasma"),
    ("MMEJ / alt-EJ — microhomology-mediated end joining",
     "PARP1, POLθ and the MRN complex anneal 3-12 bp of flanking homology. Always deletes, and "
     "the deletion size is decided by where the microhomology sits — which is why it can be "
     "predicted from the sequence.",
     "violet"),
    ("HDR — homology-directed repair",
     "RAD51 coats the resected 3' overhang and invades a template: the sister chromatid in S/G2, "
     "or a supplied ssODN/dsDNA/AAV donor. The only route to a precise change, and the reason "
     "cell-cycle stage dominates every HDR experiment.",
     "green"),
    ("SSA — single-strand annealing",
     "Long resected repeats (>30 bp) anneal to each other and the intervening DNA is deleted "
     "with the flanking repeats. The mechanism behind many of the large deletions seen in "
     "editing experiments.",
     "teal"),
]

HDR_LEVERS: List[Dict[str, Any]] = [
    {"id": "sync", "name": "Synchronise the population into S/G2",
     "effect": "raises the fraction of cells where HDR is even possible",
     "nhej": 1.0, "hdr": 1.0, "mechanism": "Cell-cycle gating: HDR simply does not exist in G1."},
    {"id": "53bp1", "name": "Inhibit 53BP1 (i53 / shRNA)",
     "effect": "releases the block on end resection",
     "nhej": 0.45, "hdr": 1.65, "mechanism": "53BP1 shields the break from resection; removing it "
                                             "shifts the balance toward BRCA1/RAD51."},
    {"id": "tether", "name": "Tether the donor to Cas9",
     "effect": "raises the local donor concentration at the break",
     "nhej": 1.0, "hdr": 1.45, "mechanism": "Donor-Cas9 fusions keep the template where the "
                                            "resection happens."},
    {"id": "cold", "name": "Cold shock / low-temperature incubation",
     "effect": "slows NHEJ more than HDR",
     "nhej": 0.9, "hdr": 1.25, "mechanism": "Empirical: 32-33 °C delays the Ku-dependent "
                                            "pathway in several cell types."},
    {"id": "rs1", "name": "RAD51 stimulation (RS-1)",
     "effect": "stabilises the RAD51 filament",
     "nhej": 1.0, "hdr": 1.30, "mechanism": "Small molecules that lock the recombinase filament "
                                            "increase templated repair."},
    {"id": "aav", "name": "Switch to an AAV6 donor",
     "effect": "better template than a naked oligo for long edits",
     "nhej": 1.0, "hdr": 1.35, "mechanism": "Viral delivery protects the donor and improves "
                                            "nuclear delivery."},
]


def render() -> None:
    st.markdown(theme.hero(
        "🩹 Repair lab",
        "The DSB is not the edit — the repair is. This bench takes the cell's decision apart and "
        "shows which levers actually move the outcome, with the sequence-level consequences "
        "computed from real DNA.",
        ["pathway biology", "microhomology finder", "HDR levers"],
    ), unsafe_allow_html=True)

    tab_tree, tab_mh, tab_spectrum, tab_hdr = st.tabs(
        ["Decision tree", "Microhomology finder", "Indel spectra by cell type", "HDR levers"])

    with tab_tree:
        st.markdown(DECISION_TREE, unsafe_allow_html=True)
        st.markdown("#### Where each pathway is allowed to act")
        phases = ["G1", "S", "G2", "M"]
        st.dataframe(pd.DataFrame([{
            "pathway": "c-NHEJ", **{p: ("yes" if p != "M" else "limited") for p in phases}},
            {"pathway": "MMEJ", **{p: ("yes" if p in ("S", "G2") else "limited") for p in phases}},
            {"pathway": "HDR", **{p: ("yes" if p in ("S", "G2") else "no") for p in phases}},
            {"pathway": "SSA", **{p: ("yes" if p in ("S", "G2") else "no") for p in phases}},
        ]), width="stretch", hide_index=True)
        for title, body, kind in BIOLOGY_CARDS:
            st.markdown(theme.card(title, body, kind), unsafe_allow_html=True)
        st.caption("Selected references: Ceccaldi 2016 Trends Cell Biol; Sfeir & Symington 2015; "
                   "Yeh 2019 Nat Cell Biol (MMEJ); Richardson 2018 Nat Genet (large deletions).")

    with tab_mh:
        st.markdown("##### Microhomology finder")
        st.caption("Pick a cut position and see exactly which deletions MMEJ can generate. This "
                   "is computed from the real sequence, not sampled from a distribution.")
        c1, c2, c3 = st.columns([1.3, 1, 1])
        with c1:
            locus_id = st.selectbox("Locus", [l["id"] for l in all_loci()],
                                    format_func=lambda i: get_locus(i)["name"])
        locus = get_locus(locus_id)
        seq = coding_strand(locus)
        with c2:
            cut = st.slider("Cut position", 30, max(31, len(seq) - 30), len(seq) // 2)
        with c3:
            min_len = st.slider("Minimum microhomology length", 2, 8, 3)
        mh = microhomology_score(seq, cut, min_len)
        st.markdown(theme.metrics([
            ("Longest microhomology", f"{mh['longest']} bp", "flanking the cut"),
            ("MMEJ substrate score", f"{mh['score']:.2f}", "0 = none, 1 = very favourable"),
        ]), unsafe_allow_html=True)

        window = seq[max(0, cut - 45):cut + 45]
        offset = max(0, cut - 45)
        hl = []
        for (l, r, length, repeat) in mh["microhomologies"][:6]:
            hl.append((l - offset + 1, l - offset + length, theme.C["gold"], f"MH {length} bp"))
            hl.append((r - offset + 1, r - offset + length, theme.C["plasma"], ""))
        st.markdown(theme.render_sequence(window, start=offset + 1, width=100, highlights=hl),
                    unsafe_allow_html=True)
        if mh["microhomologies"]:
            st.dataframe(pd.DataFrame([{
                "repeat": repeat, "length": length,
                "left position": l + 1, "right position": r + 1,
                "deleted bp": r - l,
                "predicted deletion": seq[l:r][:60],
                "flanking": f"{seq[max(0,l-8):l]}|{seq[l:r]}|{seq[r:r+8]}",
            } for l, r, length, repeat in mh["microhomologies"]]),
                width="stretch", hide_index=True)
            st.markdown(theme.callout(
                "Each row is a deletion the cell can actually produce: the two copies of the "
                "repeat anneal and everything between them is lost. Longer repeats are used "
                "preferentially, which is why MMEJ outcomes are sequence-predictable.", "info"),
                unsafe_allow_html=True)
        else:
            st.info("No microhomology of this length flanks the cut: MMEJ has nothing to anneal, "
                    "so c-NHEJ small indels dominate.")

    with tab_spectrum:
        st.markdown("##### How the cell type changes the indel spectrum")
        picks = st.multiselect("Cell types", list(CELL_BY_ID),
                               default=["hek293t", "primary_t", "ipsc", "neuron"],
                               format_func=lambda i: CELL_BY_ID[i]["name"])
        rows = []
        charts = {}
        for cid in picks:
            cell = CELL_BY_ID[cid]
            sim = simulate_indels("ACGT" * 40, 80, 1000, cell, "c-NHEJ", seed=7)
            rows.append({
                "cell type": cell["name"], "S/G2 fraction": cell["s_g2_fraction"],
                "mean indel (bp)": round(sim["mean_size"], 2),
                "1 bp indels": f"{sim['pct_1bp']:.1f}%",
                "frameshift share": f"{100*sim['frameshift_fraction']:.1f}%",
                "large deletions (≥30 bp)": f"{sim['pct_large_del']:.2f}%",
                "indel spectrum note": cell["indel_spectrum"],
            })
            charts[cell["short"] if "short" in cell else cell["name"].split(" (")[0]] = \
                [sim["size_histogram"].get(k, 0) for k in range(-10, 6)]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        if charts:
            st.plotly_chart(viz.fig_dose_response(
                [str(k) for k in range(-10, 6)], charts,
                title="Indel size profile (counts per 1000 alleles, sizes -10 … +5)",
                xlabel="indel size (bp)", ylabel="alleles"), width="stretch")
        st.markdown(theme.callout(
            "Cell type shifts the <i>composition</i> more than the total: which pathway dominates, "
            "and therefore whether you get a tidy -1 bp or a messy -30 bp deletion.", "info"),
            unsafe_allow_html=True)

    with tab_hdr:
        st.markdown("##### HDR levers — what actually moves the number")
        c1, c2 = st.columns([1.2, 1.4])
        with c1:
            cell_id = st.selectbox("Cell type", list(CELL_BY_ID),
                                   format_func=lambda i: CELL_BY_ID[i]["name"], index=0)
            cell = CELL_BY_ID[cell_id]
            template_type = st.selectbox("Donor format", ["ssODN", "dsDNA", "plasmid", "aav",
                                                          "lentivirus"], index=0)
            syncing = st.slider("Population sync into S/G2", 0.0, 1.0, 0.0, step=0.05)
            active = st.multiselect("Interventions", [l["id"] for l in HDR_LEVERS],
                                    default=["53bp1"],
                                    format_func=lambda i: next(x["name"] for x in HDR_LEVERS
                                                               if x["id"] == i))
            phase = st.selectbox("Cell-cycle phase", ["asynchronous", "G1", "S", "G2"], index=0)

        def pathway_with(levers: List[str], sync: float) -> Dict[str, float]:
            eff_cell = dict(cell)
            nhej_mult, hdr_mult = 1.0, 1.0
            for lv in HDR_LEVERS:
                if lv["id"] in levers:
                    nhej_mult *= lv["nhej"]
                    hdr_mult *= lv["hdr"]
            eff_cell["nhej_factor"] = cell.get("nhej_factor", 1.0) * nhej_mult
            eff_cell["hdr_factor"] = cell.get("hdr_factor", 1.0) * hdr_mult
            res = repair_pathways(eff_cell, phase if phase != "asynchronous" else None,
                                  template=True, template_type=template_type,
                                  mmej_score=0.2, syncing=sync)
            return res["percent"]

        base_pct = pathway_with([], 0.0)
        with c2:
            st.markdown(theme.metrics([
                ("Baseline HDR", f"{base_pct['HDR']:.1f}%", "no interventions"),
                ("Optimised HDR", f"{pathway_with(active, syncing)['HDR']:.1f}%",
                 f"{len(active)} intervention(s) active"),
                ("Fold change", f"{pathway_with(active, syncing)['HDR'] / max(0.01, base_pct['HDR']):.1f}×",
                 "compared with baseline"),
            ]), unsafe_allow_html=True)

        st.plotly_chart(viz.fig_pathway(pathway_with(active, syncing)), width="stretch")

        ladder = []
        for lv in HDR_LEVERS:
            ladder.append((lv["name"], pathway_with([lv["id"]], syncing)["HDR"]))
        ladder.append(("ALL combined", pathway_with([l["id"] for l in HDR_LEVERS], syncing)["HDR"]))
        st.plotly_chart(viz.fig_bar([n for n, _ in ladder], [v for _, v in ladder],
                                    "HDR share with each lever used alone", theme.C["green"],
                                    height=340, ytitle="HDR (% of repair events)"),
                        width="stretch")
        for lv in HDR_LEVERS:
            st.markdown(theme.card(lv["name"], f"<b>Effect:</b> {lv['effect']}<br>"
                                               f"<b>Mechanism:</b> {lv['mechanism']}", "green"),
                        unsafe_allow_html=True)
        st.markdown(theme.callout(
            "Read the numbers as directions, not promises: real HDR rates vary by locus, by "
            "delivery, and by how healthy the cells are. The one rule that survives every "
            "experiment is that HDR needs S/G2 — no donor, no synchronisation, no template "
            "choice can change that.", "info"), unsafe_allow_html=True)
