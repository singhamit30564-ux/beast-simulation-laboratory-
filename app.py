"""BEAST Simulation Laboratory — interactive CRISPR laboratory.

Run with::

    streamlit run app.py --server.address 0.0.0.0 --server.port 8501

Eight workbenches (mission control, immunity arena, Cas explorer, guide design,
genome surgery, repair lab, therapeutics, database) share one simulation core in
``beastlab``.
"""
from __future__ import annotations

import streamlit as st

from beastlab import theme
from beastlab.pages import (cas_explorer, database, genome_surgery, guide_design,
                            home, immunity_arena, repair_lab, therapeutics)

st.set_page_config(
    page_title="BEAST Simulation Laboratory",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto",
)

st.markdown(theme.CSS, unsafe_allow_html=True)

nav = st.navigation(
    {
        "Laboratory": [
            st.Page(home.render, title="Mission control", icon="🛰️", default=True,
                    url_path="mission-control"),
            st.Page(immunity_arena.render, title="Immunity arena", icon="🦠",
                    url_path="immunity-arena"),
            st.Page(genome_surgery.render, title="Genome surgery", icon="✂️",
                    url_path="genome-surgery"),
            st.Page(repair_lab.render, title="Repair lab", icon="🩹", url_path="repair-lab"),
        ],
        "Design bench": [
            st.Page(guide_design.render, title="Guide design", icon="🎯", url_path="guide-design"),
            st.Page(cas_explorer.render, title="Cas explorer", icon="🔬", url_path="cas-explorer"),
        ],
        "Translation": [
            st.Page(therapeutics.render, title="Therapeutics", icon="⚕️",
                    url_path="therapeutics"),
            st.Page(database.render, title="Database & sources", icon="📚", url_path="database"),
        ],
    }
)

with st.sidebar:
    st.markdown(
        theme.badge("BEAST v2.0", theme.C["gold"]) +
        theme.badge("in-silico only", theme.C["teal"]),
        unsafe_allow_html=True)
    st.caption(
        "A simulation laboratory. Every sequence in the genome library is real "
        "(NCBI / Ensembl); every model states its assumptions on screen."
    )

st.sidebar.toggle("Lite mode · low-memory phone", key="lite_mode", value=False)
st.sidebar.caption("Static cut frames, no supplementary Plotly charts, max 50 displayed rows. Setting persists in this session.")
def clear_session():
    lite = st.session_state.get("lite_mode", False)
    st.session_state.clear()
    st.session_state["lite_mode"] = lite


st.sidebar.button("Clear session inputs & results", on_click=clear_session)
from beastlab.field_notes import render as field_notes
field_notes(nav.title)
st.caption("Assumptions: educational, sequence-only or toy models; bench-specific details below. What this does not claim: experimental validation, clinical safety, efficacy or proof of cure.")
nav.run()
st.caption("Zero Data Retention — computed in RAM. User inputs stay in the live session until cleared or expired; no app-level disk storage. Public references alone are cached. Hosting/network logs are outside this app’s control.")

from pathlib import Path  # noqa: E402

_reference_file = Path(__file__).resolve().parent / "data" / "references.md"
_n_lines = len(_reference_file.read_text().splitlines()) if _reference_file.exists() else 0

st.markdown(
    f"""<hr style="border-color:{theme.C['line']};margin-top:38px">
    <div style="color:{theme.C['muted']};font-size:0.76rem;line-height:1.5">
    BEAST Simulation Laboratory · digital genome surgery · no wet-lab data uploaded,
    nothing here is medical advice. Traceable citations:
    <code>data/references.md</code> ({_n_lines} lines).
    </div>""",
    unsafe_allow_html=True,
)
