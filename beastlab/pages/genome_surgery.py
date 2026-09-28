"""Genome surgery — the experiment bench: cut, edit, repair, measure."""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import streamlit as st

from .. import theme, viz
from ..crispr import design_guides, design_pegRNA, find_offtargets
from ..data.cas_enzymes import BASE_EDITORS, EFFECTOR_BY_ID, NUCLEASES, PRIME_EDITORS
from ..data.targets import all_loci, coding_strand, get_locus
from ..repair import (CELL_TYPES, DELIVERY, base_editing_window, prime_editing_efficiency,
                      safety_panel, simulate_base_editing, simulate_nuclease_edit,
                      simulate_prime_editing)
from ..sequtils import clean

CELL_BY_ID = {c["id"]: c for c in CELL_TYPES}
DELIVERY_BY_ID = {d["id"]: d for d in DELIVERY}
EDITOR_BY_ID = {b["id"]: b for b in BASE_EDITORS}


def _scan_enzyme(editor: Dict[str, Any]) -> Dict[str, Any]:
    """Base and prime editors reuse a Cas nickase, but only the nucleases carry the
    PAM/site metadata the scanner needs. Fall back to a nuclease with the same PAM."""
    if "pam_side" in editor and "spacer_len" in editor:
        return editor
    pam = editor.get("pam", "NGG")
    for nuc in NUCLEASES:
        if nuc.get("pam") == pam and "spacer_len" in nuc:
            return nuc
    return EFFECTOR_BY_ID["SpCas9"]


def _best_guide(sequence: str, enzyme: Dict[str, Any], region: Tuple[int, int]) -> Optional[Dict[str, Any]]:
    guides = design_guides(sequence, _scan_enzyme(enzyme), region=region,
                           with_offtargets=False, limit=120)
    return guides[0] if guides else None


def _representative_allele(seq: str, cut: int, size: int, seed: int = 1) -> str:
    if size < 0:
        lo = max(0, cut + size)
        return seq[:lo] + seq[cut:]
    import random
    rng = random.Random(seed)
    ins = "".join(rng.choice("ACGT") for _ in range(size))
    return seq[:cut] + ins + seq[cut:]


def _gel_lanes(seq_len: int, cut: int, edited_allele_diff: int,
               hdr_done: bool) -> List[Dict[str, Any]]:
    left, right = cut, seq_len - cut
    lanes = [
        {"name": "uncut", "bands": [(seq_len, 0.95)]},
        {"name": "cut", "bands": [(left, 0.9), (right, 0.9)]},
    ]
    if edited_allele_diff:
        lanes.append({"name": "indel allele",
                      "bands": [(seq_len + edited_allele_diff, 0.85)]})
    if hdr_done:
        lanes.append({"name": "HDR allele", "bands": [(seq_len, 0.7), (left, 0.55)]})
    lanes.append({"name": "ladder", "bands": [(b, 0.35) for b in (100, 250, 500, 1000,
                                                                  2000, 5000)]})
    return lanes


def render() -> None:
    st.markdown(theme.hero(
        "✂️ Genome surgery",
        "Set up the full experiment — target, effector, cell type, delivery, dose — then look at "
        "what the cell does with the damage. Allele tables, a simulated genotyping gel, and a "
        "collateral-damage radar.",
        ["4 modalities", "8 cell types", "5 delivery routes"],
    ), unsafe_allow_html=True)

    # ------------------------------------------------------------- setup ----
    with st.form("surgery_setup"):
        st.markdown("### Experiment setup")
        r1 = st.columns([1.3, 1.2, 1])
        with r1[0]:
            locus_id = st.selectbox("Target locus", [l["id"] for l in all_loci()],
                                    format_func=lambda i: get_locus(i)["name"])
            locus = get_locus(locus_id)
            sequence = coding_strand(locus)
        with r1[1]:
            modality = st.selectbox("Modality", [
                "Nuclease knockout (NHEJ indels)",
                "HDR knock-in / precise repair",
                "Base editing (no DSB)",
                "Prime editing (no DSB)",
            ])
            cell_id = st.selectbox("Cell type", list(CELL_BY_ID),
                                   format_func=lambda i: CELL_BY_ID[i]["name"],
                                   index=4)
        with r1[2]:
            dose = st.slider("Delivered dose (fraction of optimised)", 0.05, 2.0, 1.0, step=0.05)
            n_cells = st.select_slider("Alleles simulated", [200, 500, 1000, 2000], value=1000)

        r2 = st.columns([1.2, 1.1, 1.2])
        with r2[0]:
            delivery = st.selectbox("Delivery route", list(DELIVERY_BY_ID),
                                    format_func=lambda i: DELIVERY_BY_ID[i]["name"],
                                    index=0)
        with r2[1]:
            phase = st.selectbox("Cell-cycle context",
                                 ["asynchronous", "G1", "S", "G2", "M"], index=0)
            syncing = st.slider("Population sync into S/G2", 0.0, 1.0, 0.0, step=0.05)
        with r2[2]:
            if modality.startswith("Nuclease") or modality.startswith("HDR"):
                effector_id = st.selectbox("Nuclease", [e["id"] for e in NUCLEASES],
                                           index=0,
                                           format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
            elif modality.startswith("Base"):
                effector_id = st.selectbox("Base editor", [b["id"] for b in BASE_EDITORS],
                                           format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
            else:
                effector_id = st.selectbox("Prime editor", [p["id"] for p in PRIME_EDITORS],
                                           index=1,
                                           format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
            enzyme = EFFECTOR_BY_ID[effector_id]

        r3 = st.columns([1, 1, 1])
        with r3[0]:
            template_type = st.selectbox("Donor template (if any)",
                                         ["none", "ssODN", "dsDNA", "plasmid", "aav", "lentivirus"],
                                         index=1 if modality.startswith("HDR") else 0)
        with r3[1]:
            desired_edit = st.text_input("Sequence to write (optional)",
                                         placeholder="e.g. T for the sickle correction")
        with r3[2]:
            windows_bp = st.number_input("Amplicon length for the gel (bp)", 300, 5000, 1200, 100)

        submitted = st.form_submit_button("▶ Run the experiment", type="primary")

    if not submitted and "surgery_result" not in st.session_state:
        st.info("Configure the experiment and press **Run**. Everything below is generated "
                "on demand from the shipped sequences.")
        return

    if submitted:
        region = (0, len(sequence))
        guide = _best_guide(sequence, enzyme, region)
        if guide is None:
            st.error(f"No {enzyme['pam']} PAM in this locus for {enzyme['name']}. "
                     f"Pick a PAM-flexible variant (SpG, SpRY) or a different locus.")
            return
        cut = guide["site"]["cut_left"]
        st.session_state["surgery_result"] = {
            "guide": guide, "cut": cut, "sequence": sequence, "locus": locus,
            "modality": modality, "cell": CELL_BY_ID[cell_id], "dose": dose,
            "phase": phase, "syncing": syncing, "template": template_type,
            "desired_edit": clean(desired_edit) if desired_edit else "",
            "n_cells": n_cells, "enzyme": enzyme, "amplicon": int(windows_bp),
            "delivery": DELIVERY_BY_ID[delivery],
        }

    R = st.session_state["surgery_result"]
    locus, sequence, enzyme = R["locus"], R["sequence"], R["enzyme"]
    cell, cut = R["cell"], R["cut"]
    modality, dose, phase = R["modality"], R["dose"], R["phase"]
    n_cells = R["n_cells"]

    # ------------------------------------------------------------- context --
    st.markdown("### The cut")
    start = max(0, cut - 30)
    end = min(len(sequence), cut + 30)
    gsite = R["guide"]["site"]
    hl = [(gsite["start"] - start + 1, gsite["end"] - start, theme.C["violet"], "protospacer"),
          (gsite["pam_start"] - start + 1, gsite["pam_end"] - start, theme.C["gold"], "PAM")]
    st.markdown(theme.render_sequence(sequence[start:end], start=start + 1, width=110,
                                      highlights=hl), unsafe_allow_html=True)
    st.markdown(theme.kv({
        "locus": locus["name"], "cut position": f"{cut} (1-based)",
        "guide": R["guide"]["guide"], "PAM": R["guide"]["pam"],
        "effector": enzyme["name"], "modality": modality,
        "cell": cell["name"], "delivery": R["delivery"]["name"],
    }), unsafe_allow_html=True)

    off_hits = find_offtargets(sequence, R["guide"]["guide"], enzyme, 3)
    n_dsbs = 1.0 if modality.startswith(("Nuclease", "HDR")) else 0.0

    # ------------------------------------------------------- nuclease path --
    if modality.startswith("Nuclease") or modality.startswith("HDR"):
        with st.spinner("simulating repair…"):
            result = simulate_nuclease_edit(
                sequence, cut, cell, dose=dose, phase=phase,
                template=(modality.startswith("HDR") or R["template"] != "none"),
                template_type=R["template"] if R["template"] != "none" else "ssODN",
                desired_edit=R["desired_edit"], n_cells=n_cells)
        pct = result["outcome_percent"]

        st.markdown("### Outcome")
        st.markdown(theme.metrics([
            ("Indels", f"{100*result['indel_fraction']:.1f}%", "of all alleles carry an indel"),
            ("Precise HDR", f"{100*result['hdr_precise_fraction']:.1f}%",
             "templated, error-free repair"),
            ("Biallelic knockout", f"{100*result['biallelic_knockout_fraction']:.1f}%",
             "both alleles disrupted (frameshift)"),
            ("Mean indel size", f"{result['indels']['mean_size']:.1f} bp", "population average"),
            ("Frameshift share", f"{100*result['indels']['frameshift_fraction']:.0f}%",
             "indels that shift the reading frame"),
        ]), unsafe_allow_html=True)

        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(viz.fig_pathway(pct), width="stretch")
        with c2:
            st.plotly_chart(viz.fig_indel_spectrum(result["indels"]["size_histogram"]),
                            width="stretch")
        st.plotly_chart(viz.fig_allele_pie({k: v for k, v in result["counts"].items() if v},
                                           "Allele fate of the simulated population"),
                        width="stretch")

        st.markdown("#### Why the pathway split came out this way")
        for line in result["pathways"]["reasoning"]:
            st.markdown(f"- {line}")

        # sequence-level outcome
        modes = sorted(result["indels"]["size_histogram"].items(),
                       key=lambda kv: -kv[1])
        if modes:
            size = int(modes[0][0])
            edited = _representative_allele(sequence, cut, size)
            st.markdown("#### Most common allele, at the sequence level")
            st.markdown(viz.sequence_diff_html(sequence, edited, context=22),
                        unsafe_allow_html=True)
            st.caption(f"The dominant product here is a {abs(size)} bp "
                       f"{'deletion' if size < 0 else 'insertion'} — the modal outcome of the "
                       f"drawn spectrum, shown against the reference.")

        if R["desired_edit"]:
            hdr_seq = sequence[:cut] + R["desired_edit"] + sequence[cut:]
            st.markdown("#### The intended HDR product")
            st.markdown(viz.sequence_diff_html(sequence, hdr_seq, context=22),
                        unsafe_allow_html=True)

        with st.expander("Sequence-specific microhomology deletions (computed from this locus)"):
            if result["mmej_predictions"]:
                st.dataframe(pd.DataFrame([{
                    "microhomology": p["repeat"], "length": p["repeat_len"],
                    "deleted bp": p["deletion_len"],
                    "deleted sequence": p["deleted_sequence"][:48],
                } for p in result["mmej_predictions"]]), width="stretch", hide_index=True)
                st.caption("These deletions are not random: they are exactly the products that "
                           "annealing the flanking microhomologies would generate.")
            else:
                st.write("No microhomology ≥3 bp flanks this cut — MMEJ has little to work with, "
                         "so small NHEJ indels dominate.")

        # gel
        st.markdown("#### Simulated genotyping gel")
        diff = int(modes[0][0]) if modes else 0
        st.plotly_chart(viz.fig_gel(_gel_lanes(R["amplicon"], max(1, cut % R["amplicon"]),
                                               diff, result["hdr_precise_fraction"] > 0.02)),
                        width="stretch")

        # safety
        st.markdown("#### Collateral damage")
        safety = safety_panel(cell, n_dsbs=max(n_dsbs, 1.0),
                              off_target_hits=max(1, len(off_hits)),
                              on_target_indel_fraction=result["indel_fraction"],
                              dose=dose,
                              hours_expressed=R["delivery"]["duration_h"] * dose)
        c3, c4 = st.columns([1, 1])
        with c3:
            st.plotly_chart(viz.fig_radar({
                "off-target indels": min(1.0, safety["off_target_indel_fraction"] * 8),
                "large deletions": min(1.0, safety["large_deletion_fraction"] * 4),
                "translocations": min(1.0, safety["translocation_fraction"] * 6),
                "p53 response": safety["p53_response"],
                "senescence": safety["senescence_or_apoptosis"],
                "on-target indels": min(1.0, safety["on_target_indel_fraction"]),
            }), width="stretch")
        with c4:
            st.markdown(theme.kv({
                "off-target indels": f"{100*safety['off_target_indel_fraction']:.2f}%",
                "large deletions": f"{100*safety['large_deletion_fraction']:.2f}%",
                "translocations": f"{100*safety['translocation_fraction']:.3f}%",
                "p53 response": f"{100*safety['p53_response']:.0f}%" if cell["p53_competent"] else "not applicable",
                "therapeutic index": safety["therapeutic_index"],
            }), unsafe_allow_html=True)
            for note in safety["notes"]:
                st.markdown(f"- {note}")

    # ----------------------------------------------------- base editing -----
    elif modality.startswith("Base"):
        editor = enzyme
        st.markdown("### Base editing window")
        st.caption(f"{editor['name']} converts {editor['conversion']} in protospacer positions "
                   f"{editor['window'][0]}–{editor['window'][1]}.")
        gsite = R["guide"]["site"]
        site = gsite
        editable = base_editing_window(site, editor)
        if not editable:
            st.warning("No editable base of the right identity sits in this editor's window at "
                       "this site. Move the target window or choose a different editor — this is "
                       "the single most common failure mode in base editing design.")
            return
        st.dataframe(pd.DataFrame([{
            "position": e["position"], "base": e["base"], "context": e["context"],
            "relative efficiency": e["relative_efficiency"],
            "bystander": "yes" if e["is_bystander"] else "target (best)",
        } for e in editable]), width="stretch", hide_index=True)

        with st.spinner("simulating conversion…"):
            be = simulate_base_editing(site, editor, cell, dose=dose, n_cells=n_cells)
        st.markdown(theme.metrics([
            ("Editing efficiency", f"{be['editing_efficiency']:.1f}%", "alleles with ≥1 conversion"),
            ("Product purity", f"{be['product_purity']:.1f}%", "alleles with a single conversion"),
            ("Indels", f"{be['indel_rate_percent']:.2f}%", "DSB-free, so far lower than nuclease"),
        ]), unsafe_allow_html=True)
        c1, c2 = st.columns(2)
        with c1:
            st.plotly_chart(viz.fig_bar([str(k) for k in be["per_position_conversion"]],
                                        list(be["per_position_conversion"].values()),
                                        "Per-position conversion (%)", theme.C["teal"], 320),
                            width="stretch")
        with c2:
            st.plotly_chart(viz.fig_allele_pie(be["allele_composition"]), width="stretch")
        st.markdown(theme.callout(be["transversion_risk"], "info"), unsafe_allow_html=True)

        best_edit = max(editable, key=lambda e: e["relative_efficiency"])
        target_base = best_edit["base"]
        new_base = "T" if target_base == "C" else "G"
        pos_in_locus = site["start"] + best_edit["position"] - 1
        edited_seq = sequence[:pos_in_locus] + new_base + sequence[pos_in_locus + 1:]
        st.markdown("#### The edit, on the chromosome")
        st.markdown(viz.sequence_diff_html(sequence, edited_seq, context=22),
                    unsafe_allow_html=True)

        safety = safety_panel(cell, n_dsbs=0.0, off_target_hits=max(1, len(off_hits)),
                              on_target_indel_fraction=be["indel_rate_percent"] / 100.0,
                              dose=dose, hours_expressed=R["delivery"]["duration_h"] * dose)
        st.markdown(theme.callout(
            f"No double-strand break, so the translocation term is exactly zero and large "
            f"deletions are essentially absent. The residual risks are RNA off-target editing "
            f"and, for CBE, cytosine transversions. Off-target DNA burden from the nickase "
            f"component stays at ~{100*safety['off_target_indel_fraction']:.2f}%.", "good"),
            unsafe_allow_html=True)

    # ---------------------------------------------------- prime editing -----
    else:
        editor = enzyme
        st.markdown("### Prime editing")
        gsite = R["guide"]["site"]
        edit = R["desired_edit"] or "AG"
        p1, p2, p3 = st.columns(3)
        with p1:
            pbs_len = st.slider("PBS length", 6, 18, 13)
        with p2:
            rtt_len = st.slider("RT template length", 8, 40, 16)
        with p3:
            pe3 = st.toggle("Add a second, PE3-style nick", value=False)
        peg = design_pegRNA(sequence, R["guide"]["guide"], enzyme, edit,
                            edit_offset=len(R["guide"]["guide"]) // 2,
                            pbs_len=pbs_len, rtt_len=rtt_len, add_motif=True)
        eff = prime_editing_efficiency(peg["pbs_tm"], peg["pbs_gc"], peg["rtt_len"], editor,
                                       cell, pe3)
        sim = simulate_prime_editing(sequence, peg["nick_1"] - 1, edit, peg, editor, cell,
                                     pe3=pe3, n_cells=n_cells)
        st.markdown(theme.metrics([
            ("Precise edit", f"{sim['percent']['precise_edit']:.1f}%", "exactly the intended change"),
            ("Indels", f"{sim['percent']['indel']:.2f}%", "flap-resolution byproducts"),
            ("Unedited", f"{sim['percent']['unedited']:.1f}%", "no change"),
            ("PBS Tm", f"{peg['pbs_tm']:.1f} °C", f"PBS {peg['pbs']}"),
        ]), unsafe_allow_html=True)
        st.markdown(theme.kv({
            "pegRNA": f"{len(peg['pegRNA'])} nt",
            "RT template": peg["rtt"],
            "nick": peg["nick_1"],
            "3' motif": "present" if peg["motif"] else "absent",
        }), unsafe_allow_html=True)
        st.markdown(viz.sequence_diff_html(sequence, sim["edited_sequence"], context=22),
                    unsafe_allow_html=True)
        for note in eff["notes"]:
            st.markdown(f"- {note}")

    # ------------------------------------------------------------ narrative -
    st.markdown("### What just happened, in order")
    steps = [
        f"The editor ({enzyme['name']}) is delivered by {R['delivery']['name']} — onset about "
        f"{R['delivery']['onset_h']} h, activity lasting about "
        f"{'indefinitely' if R['delivery']['duration_h'] > 1e5 else str(int(R['delivery']['duration_h'])) + ' h'}.",
        f"It locates {R['guide']['guide']} next to {R['guide']['pam']} and invades the duplex.",
    ]
    if modality.lower().startswith("base") or modality.lower().startswith("prime"):
        steps.append("No double-strand break is made: the chemistry (deamination) or the "
                     "reverse transcriptase writes the change at the nick.")
        steps.append("Consequence: translocation, large-deletion and p53 terms stay at zero — "
                     "but the edit is limited to what the editor's chemistry can write.")
    else:
        steps.append(f"The DSB at position {cut} is detected by the cell and sent down one of "
                     f"four repair routes.")
        steps.append("In this cell type the split is dominated by c-NHEJ, which is why indels, "
                     "not precision, are the default outcome.")
        steps.append("Every allele that survives with a frameshift contributes to the knockout "
                     "fraction.")
    for i, step in enumerate(steps, 1):
        st.markdown(f"**{i}.** {step}")
    st.caption("Model assumptions: " + " · ".join(result["assumptions"])
               if modality.startswith(("Nuclease", "HDR")) else
               "Base/prime editing outcomes come from the editor parameters in the Cas explorer "
               "plus a sequence context model for the window.")
