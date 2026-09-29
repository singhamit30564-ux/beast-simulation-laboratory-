"""Low-bandwidth presentation helpers; raw results remain available for export."""
import streamlit as st


def table(data, *args, **kwargs):
    if st.session_state.get("lite_mode", False):
        n = len(data)
        data = data.iloc[:50] if hasattr(data, "iloc") else data[:50]
        if n > 50:
            st.caption(f"Lite mode: showing 50 of {n:,} rows. Calculations are unchanged.")
    return st.dataframe(data, *args, **kwargs)


def chart(fig, *args, **kwargs):
    if st.session_state.get("lite_mode", False):
        st.caption("Lite mode: supplementary chart skipped; numeric outputs remain below.")
        return None
    return st.plotly_chart(fig, *args, **kwargs)
