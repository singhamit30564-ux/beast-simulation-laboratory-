"""Flagship sickle-cell teaching workflow; never a patient-level predictor."""
import streamlit as st
from .data.targets import get_locus, coding_strand
from .data.clinical import editing_to_hbf
from .io import parse_dna, export_center
from .offtarget import scan

# Synthetic enhancer fragment: annotated classroom motif, not a clinical guide.
DEMO = "AACCTT" + "ACGTACGTACGTACAGATAA" + "TGG" + "AACCTTAGCTAGCTAGCTAA"
DEMO_GUIDE = DEMO[6:26]


def experiment(sequence, guide, deletion, edited_fraction, motif_start):
    hits = scan(sequence, guide)
    exact = [h for h in hits if h["mismatches"] == 0]
    # Projection from reported 1-based forward protospacer interval.
    hit = exact[0] if exact else None
    cut = None
    if hit:
        start = int(hit["position"].split(":")[1].split("-")[0])-1
        cut = start+17 if hit["strand"] == "+" else start+3
    lo = max(0, cut-deletion//2) if cut is not None else 0
    hi = min(len(sequence), lo+deletion) if cut is not None else 0
    disrupted = bool(hit and deletion > 0 and lo < motif_start+6 and hi > motif_start)
    hbf = editing_to_hbf(edited_fraction if disrupted else 0)
    return {"binding_site": hit, "cut_boundary_0based": cut, "deletion_interval_0based": [lo, hi],
            "edited_allele": sequence[:lo]+sequence[hi:] if hit else sequence,
            "motif_disrupted_assumption": disrupted, "assumed_edited_cells_percent": edited_fraction,
            "predicted_HbF_percent": round(hbf, 2),
            "verdict": "✅ Functional cure predicted" if disrupted and hbf >= 30 else "⚠️ needs work",
            "verdict_scope": "Toy classroom threshold HbF >=30%, NOT a clinical outcome or proof of cure",
            "assumptions": "Exact 20nt match + adjacent NGG; first exact site; chosen deletion; motif overlap implies disruption; HbF=6+54*(1-exp(-3*edited_fraction/100)). Ignores delivery, mosaicism, HbF distribution, safety and durability."}


def render():
    with st.expander("🩸 Sickle cell case study · BCL11A → HbF", expanded=True):
        st.warning("Educational simulation only. Even the ‘Functional cure predicted’ card is a classroom verdict, not a medical prediction.")
        st.markdown("**1 · Load HBB disease context**")
        hbb = get_locus("hbb_sickle")
        st.caption(f"Bundled reference: {hbb['locus']}. Reference sequence does not establish a patient's HbS genotype.")
        st.code(coding_strand(hbb)[:100] + "…", language=None)
        st.write("HbS disease involves HBB. This experiment does not correct HBB: it disrupts the erythroid BCL11A enhancer to model release of fetal haemoglobin (HbF) repression.")
        st.markdown("**2 · Design a gRNA for the BCL11A enhancer**")
        mode = st.radio("Case sequence", ["Labeled synthetic demo", "User enhancer fragment"], key="case_mode")
        if mode == "Labeled synthetic demo":
            sequence, guide = DEMO, DEMO_GUIDE
            st.caption("Synthetic teaching sequence, not a patient sequence or validated BCL11A guide. AGATAA is the annotated teaching motif.")
            st.code(sequence, language=None)
            guide = st.text_input("Case gRNA (20 bases, no PAM)", value=guide)
        else:
            sequence = st.text_area("Enhancer FASTA/raw DNA (max 10,000 bases)", key="case_sequence")
            guide = st.text_input("Case gRNA (20 bases, no PAM)", value="")
        if mode == "Labeled synthetic demo":
            st.caption("Demo design: the 20 bases immediately before TGG are selected; the predicted cut overlaps the teaching motif after the chosen deletion.")
        try:
            sequence = parse_dna(sequence)
            guide = parse_dna(guide, guide=True)
        except ValueError as exc:
            st.info(str(exc)); return
        motif = sequence.find("AGATAA")
        if motif < 0:
            st.warning("No annotated AGATAA teaching motif. This case cannot infer enhancer disruption.")
            return
        st.caption(f"Teaching motif: bases {motif+1}–{motif+6}, first AGATAA. Motif presence alone does not establish functional enhancer identity.")
        deletion = st.slider("Assumed NHEJ deletion size (bp)", 0, 12, 6, key="case_deletion")
        fraction = st.slider("Assumed edited erythroid cells (%)", 0, 100, 75, key="case_fraction")
        fingerprint = (sequence, guide, deletion, fraction)
        if st.button("Run case: Cas9 cut → NHEJ → HbF", key="case_run"):
            st.session_state["case_result"] = (fingerprint, experiment(sequence, guide, deletion, fraction, motif))
        saved = st.session_state.get("case_result")
        if saved and saved[0] == fingerprint:
            result = saved[1]
            st.markdown("**3 · Cas9 cut → 4 · NHEJ indel**")
            if result["binding_site"]:
                st.write(f"Exact match + NGG; cut boundary {result['cut_boundary_0based']} (0-based). Chosen deletion interval {result['deletion_interval_0based']}.")
                st.code(result["edited_allele"], language=None)
            else:
                st.warning("No exact guide + adjacent NGG site: no simulated cut.")
            st.markdown("**5 · HbF reactivation → 6 · Verdict**")
            st.metric("Toy HbF prediction", f"{result['predicted_HbF_percent']}%")
            st.info(result["verdict"] + " — educational simulation only")
            st.caption("Assumptions: " + result["assumptions"])
            st.caption("What this does not claim: clinical efficacy, safety, true enhancer activity or a patient's cure. The 30% threshold is a teaching rule; no single HbF threshold proves a cure.")
            export_center({"sequence": sequence, "guide": guide}, result, mode, "case")
