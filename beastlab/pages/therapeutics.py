"""Therapeutics — from an edit to a clinical endpoint, with the trial data attached."""
from __future__ import annotations

from typing import Dict, List

import pandas as pd
import streamlit as st
from ..presentation import table, chart

from .. import theme, viz
from ..data.clinical import (DISEASES, TRIALS, editing_to_hbf, ldl_response,
                             simulate_cohort, ttr_knockdown, voc_rate)
from ..data.targets import get_locus
from ..sequtils import gc_content

DISEASE_LOCUS = {"scd": "bcl11a_enhancer", "tbt": "hbg1_promoter", "attr": "ttr_exon3",
                 "fh": "pcsk9_start", "lca10": "cep290_ivs26"}
DISEASE_EDIT = {
    "scd": "Disrupt the erythroid enhancer of BCL11A (or edit the HBG1/2 promoter) so that "
           "foetal haemoglobin is made again after birth.",
    "tbt": "The same HbF de-repression: replace the missing beta-globin with a globin chain "
           "that still pairs with alpha.",
    "attr": "Knock the TTR gene out in hepatocytes — the protein is made almost only in the "
            "liver, so a knockout is a whole-body knockdown.",
    "fh": "Install a stop codon in PCSK9 (base editing) so the receptor that clears LDL is not "
          "degraded.",
    "lca10": "Cut inside the CEP290 intron 26 mutation to delete the cryptic splice site and "
             "restore the correct transcript in photoreceptors.",
}
OUTCOME_KEY = {"scd": "voc_per_year", "tbt": "hbF_percent", "attr": "ttr_mg_dL",
               "fh": "LDL_mg_dL", "lca10": "vision_score"}
OUTCOME_LABEL = {"scd": "vaso-occlusive crises per year", "tbt": "HbF (%)",
                 "attr": "serum TTR (mg/dL)", "fh": "LDL-C (mg/dL)",
                 "lca10": "visual function score (model units)"}


def _response_curves(disease: str) -> Dict[str, List[float]]:
    """Sweep editing 0-100% and report the modelled pharmacodynamic read-outs."""
    xs = list(range(0, 101, 5))
    if disease in ("scd", "tbt"):
        hbf = [editing_to_hbf(x, DISEASES[disease]["baseline"]["HbF_percent"]) for x in xs]
        series = {"HbF (%)": hbf}
        if disease == "scd":
            series["VOC / year"] = [voc_rate(h, DISEASES["scd"]["baseline"]["VOC_per_year"])
                                    * 10 for h in hbf]
        return series
    if disease == "attr":
        return {"TTR (mg/dL)": [ttr_knockdown(x, DISEASES["attr"]["baseline"]["TTR_mg_dL"])
                                for x in xs]}
    if disease == "fh":
        ldl = [ldl_response(x, DISEASES["fh"]["baseline"]["LDL_mg_dL"])[0] for x in xs]
        pcs = [ldl_response(x, DISEASES["fh"]["baseline"]["LDL_mg_dL"])[1] for x in xs]
        return {"LDL-C (mg/dL)": ldl, "PCSK9 (% of baseline)": pcs}
    return {"visual function score": [DISEASES["lca10"]["baseline"]["vision_score"]
                                      * (1 + 3.5 * min(0.75, x / 100.0 * 1.4))
                                      for x in xs]}


def render() -> None:
    from ..case_study import render as sickle_case
    sickle_case()
    st.markdown(theme.hero(
        "💊 Therapeutics",
        "Every programme on this page is a real trial or an approved medicine. Move the editing "
        "efficiency slider and watch the clinical read-out that the regulators actually measured "
        "move with it.",
        ["5 indications", "7 programmes", "cohort projection"],
    ), unsafe_allow_html=True)

    tab_curve, tab_cohort, tab_trials = st.tabs(
        ["Efficiency → endpoint", "Patient cohort", "The programmes"])

    # ------------------------------------------------------- dose response --
    with tab_curve:
        c1, c2 = st.columns([1, 2])
        with c1:
            disease = st.selectbox("Indication", list(DISEASES),
                                   format_func=lambda k: DISEASES[k]["name"])
            d = DISEASES[disease]
            editing = st.slider("Editing efficiency (% of alleles)", 0, 100, 60)
            st.markdown(theme.kv({
                "programme": d["market"],
                "cell target": get_locus(DISEASE_LOCUS[disease])["name"],
                "regulatory endpoint": d["endpoint"],
                "reference": d["ref"],
            }), unsafe_allow_html=True)
        with c2:
            curves = _response_curves(disease)
            chart(viz.fig_dose_response(
                list(range(0, 101, 5)), curves,
                title=f"{d['name']} — modelled endpoint response",
                xlabel="editing efficiency (%)", ylabel="value"), width="stretch")

        st.markdown("#### At the chosen efficiency")
        if disease in ("scd", "tbt"):
            hbf = editing_to_hbf(editing, d["baseline"]["HbF_percent"])
            items = [("HbF", f"{hbf:.1f}%", f"baseline {d['baseline']['HbF_percent']:.0f}%"),
                     ("Pan-cellular?", "yes" if hbf > 20 else "partial",
                      "benefit rises steeply above ~20-30%")]
            if disease == "scd":
                voc = voc_rate(hbf, d["baseline"]["VOC_per_year"])
                items.append(("VOC / year", f"{voc:.2f}",
                              f"baseline {d['baseline']['VOC_per_year']:.1f}"))
            st.markdown(theme.metrics(items), unsafe_allow_html=True)
        elif disease == "attr":
            ttr = ttr_knockdown(editing, d["baseline"]["TTR_mg_dL"])
            st.markdown(theme.metrics([
                ("Serum TTR", f"{ttr:.1f} mg/dL", f"baseline {d['baseline']['TTR_mg_dL']:.0f}"),
                ("Knockdown", f"{100*(1 - ttr/d['baseline']['TTR_mg_dL']):.0f}%",
                 "target of the in-vivo programmes: ~90%"),
            ]), unsafe_allow_html=True)
        elif disease == "fh":
            ldl, pcs = ldl_response(editing, d["baseline"]["LDL_mg_dL"])
            st.markdown(theme.metrics([
                ("LDL-C", f"{ldl:.0f} mg/dL", f"baseline {d['baseline']['LDL_mg_dL']:.0f}"),
                ("PCSK9", f"{pcs:.0f}%", "of baseline protein"),
                ("Reduction", f"{100*(1 - ldl/d['baseline']['LDL_mg_dL']):.0f}%",
                 "durable after a single dose"),
            ]), unsafe_allow_html=True)
        else:
            rescue = min(0.75, editing / 100.0 * 1.4)
            st.markdown(theme.metrics([
                ("Photoreceptors rescued", f"{100*rescue:.0f}%", "modelled fraction with a "
                                                                 "corrected transcript"),
                ("Light sensitivity", "improved" if rescue > 0.15 else "unlikely to be measurable",
                 "the BRILLIANCE endpoint was exactly this threshold"),
            ]), unsafe_allow_html=True)

        st.markdown(theme.callout(d["response"], "info"), unsafe_allow_html=True)
        st.markdown(theme.card("What the edit does in this disease",
                               f"{DISEASE_EDIT[disease]}", "teal"), unsafe_allow_html=True)
        st.caption("Caveat: this is an illustrative mapping built from the published "
                   "population-level relationship between the biomarker and the endpoint. It "
                   "does not model trial eligibility, conditioning, or comedication.")

    # ------------------------------------------------------------ cohort ----
    with tab_cohort:
        c1, c2 = st.columns([1, 1.6])
        with c1:
            disease = st.selectbox("Indication", list(DISEASES), key="cohort_disease",
                                   format_func=lambda k: DISEASES[k]["name"])
            mean_editing = st.slider("Mean editing across the cohort (%)", 5, 95, 70)
            variability = st.slider("Between-patient variability (CV %)", 5, 60, 28,
                                    help="How uneven the editing is between patients — the "
                                         "single biggest reason a good average still leaves "
                                         "non-responders.")
            n_patients = st.slider("Simulated patients", 10, 200, 40, step=10)
            seed = st.number_input("Random seed", 0, 9999, 13)
            run = st.button("👥 Simulate the cohort", type="primary")
        result = None
        if run or "cohort_result" in st.session_state:
            if run:
                st.session_state["cohort_result"] = simulate_cohort(
                    disease, mean_editing, n_patients, variability / 100.0, int(seed))
                st.session_state["cohort_disease_key"] = disease
            result = st.session_state["cohort_result"]
        if result is None:
            st.info("Press **Simulate the cohort** to project patient-level outcomes.")
        else:
            df = pd.DataFrame(result["rows"])
            key = OUTCOME_KEY[disease]
            with c2:
                st.markdown(theme.metrics([
                    ("Response rate", f"{result['response_rate_percent']:.0f}%",
                     "of simulated patients"),
                    ("Mean editing", f"{result['mean_editing_percent']:.1f}%", "realised"),
                    ("Non-responders", f"{result['n'] - round(result['response_rate_percent']*result['n']/100)}",
                     "patients below the endpoint threshold"),
                ]), unsafe_allow_html=True)
                chart(viz.fig_histogram(df["editing_percent"].tolist(),
                                                  "Distribution of realised editing",
                                                  "editing efficiency (%)"), width="stretch")
            chart(viz.fig_scatter(
                df["editing_percent"].tolist(), df[key].tolist(),
                f"{OUTCOME_LABEL[disease]} versus realised editing", "editing efficiency (%)",
                OUTCOME_LABEL[disease], color=theme.C["gold"]), width="stretch")
            table(df.round(2), width="stretch", hide_index=True, height=260)
            st.markdown(theme.callout(result["model"], "info"), unsafe_allow_html=True)
            st.markdown(theme.callout(result["disclaimer"], "warn"), unsafe_allow_html=True)

    # ------------------------------------------------------------ trials ----
    with tab_trials:
        st.markdown("#### Programmes in the clinic")
        modalities = sorted({t["modality"] for t in TRIALS})
        pick = st.multiselect("Filter by modality", modalities, default=[])
        rows = [t for t in TRIALS if not pick or t["modality"] in pick]
        table(pd.DataFrame([{
            "product": t["product"], "sponsor": t["sponsor"], "modality": t["modality"],
            "target": t["target"], "indication": t["indication"], "phase": t["phase"],
            "delivery": t["delivery"],
        } for t in rows]), width="stretch", hide_index=True, height=300)
        st.caption(f"{len(rows)} of {len(TRIALS)} programmes shown · data as compiled in "
                   f"`data/references.md`.")

        st.markdown("#### Programme detail")
        trial_id = st.selectbox("Programme", [t["id"] for t in TRIALS],
                                format_func=lambda i: next(t["product"] for t in TRIALS
                                                           if t["id"] == i))
        t = next(x for x in TRIALS if x["id"] == trial_id)
        st.markdown(theme.card(t["product"], f"<b>{t['indication']}</b><br>{t['outcomes']}",
                               "gold"), unsafe_allow_html=True)
        st.markdown(theme.kv({
            "sponsor": t["sponsor"], "modality": t["modality"], "target": t["target"],
            "edit": t["edit"], "delivery": t["delivery"], "phase": t["phase"],
            "patients": t["n_patients"], "dose": t["dose"], "reference": t["ref"],
        }), unsafe_allow_html=True)
        st.markdown(theme.callout(t["notes"], "info"), unsafe_allow_html=True)

        # link the trial back to the editing mechanics
        locus_id = DISEASE_LOCUS.get(
            next((k for k, v in DISEASES.items() if v["name"].split(" (")[0].lower() in
                  t["indication"].lower()), ""), "")
        if locus_id:
            locus = get_locus(locus_id)
            st.markdown(f"**In this app:** the target sits in the *{locus['name']}* record — "
                        f"{locus['length']:,} bp of real sequence "
                        f"({100*gc_content(locus['sequence']):.1f}% GC). Open "
                        f"**Guide design bench** → *Genome library* → `{locus_id}` to design the "
                        f"actual guides for it.")
        st.markdown(theme.callout(
            "Approval status is a moving target — check the current label or trial registry "
            "before quoting any of this in a document.", "warn"), unsafe_allow_html=True)

    st.markdown("---")
    st.caption("The clinical layer is a model layered over published averages. The sequence "
               "layer below it (guide design, repair simulation, immunity) is where you can "
               "change the inputs and see the consequences yourself.")
