"""Six SVG-rendered Plotly frames: no WebGL, no auto-play."""
import textwrap
import plotly.graph_objects as go
import streamlit as st

STEPS = ["Cas9 approaches", "PAM scan", "DNA unzip", "Guide RNA locks", "CUT — scissor snap", "Repair begins"]
LINES = ["Find the doorway before turning the key.", "PAM is the doorman — no PAM, no entry.",
         "Open a small window, not the whole genome.", "Pairing checks the address.",
         "The scissors act; the cell decides what follows.", "Repair can leave many different endings."]


def frame_data(step):
    xs = list(range(12))
    upper = [1 + (0.5 if step >= 2 and 4 <= x <= 8 else 0) for x in xs]
    lower = [-v for v in upper]
    # Fixed trace topology ensures reliable animation, including the broken backbone.
    data = []
    for ys, color in ((upper, "#52c7ca"), (lower, "#b2a0ff")):
        xx, yy = [], []
        for x, y in zip(xs, ys):
            if x == 7 and step >= 4:
                xx.append(None); yy.append(None)
            xx.append(x); yy.append(y)
        data.append(go.Scatter(x=xx, y=yy, mode="lines+markers", marker=dict(size=13), line=dict(color=color, width=3)))
    bx, by = [], []
    for x in xs:
        if step < 2 or x < 4 or x > 8:
            bx += [x, x, None]; by += [upper[x], lower[x], None]
    data.append(go.Scatter(x=bx, y=by, mode="lines", line=dict(color="#8995a5", width=1)))
    data.append(go.Scatter(x=[2 if step == 0 else 7], y=[3], mode="markers+text", text=["Cas9"], textposition="middle center", marker=dict(size=65, color="#e8be62")))
    data.append(go.Scatter(x=[9, 10, 11], y=[upper[9], upper[10], upper[11]], mode="markers", marker=dict(size=20, symbol="circle-open", color="#e8be62"), opacity=1 if step >= 1 else 0))
    data.append(go.Scatter(x=[4, 5, 6, 7, 8], y=[0.8]*5, mode="lines+markers", line=dict(color="#ef8a61", width=4), opacity=1 if step >= 3 else 0))
    data.append(go.Scatter(x=[6.5], y=[0], mode="text", text=["✂" if step == 4 else "repair → indel" if step == 5 else ""], textfont=dict(size=25)))
    return data


def figure(step=0, animated=True):
    fig = go.Figure(data=frame_data(step))
    fig.update_layout(height=330, margin=dict(l=5, r=5, t=90, b=60), showlegend=False,
                      title=dict(text=f"{step+1}. {STEPS[step]}<br><sup>Dr. Titan: {'<br>'.join(textwrap.wrap(LINES[step], 42))}</sup>", font_size=15),
                      xaxis=dict(range=[-1, 12], visible=False, fixedrange=True),
                      yaxis=dict(range=[-2.5, 4], visible=False, fixedrange=True))
    if animated:
        fig.frames = [go.Frame(name=str(i), data=frame_data(i), layout=dict(title=dict(text=f"{i+1}. {STEPS[i]}<br><sup>Dr. Titan: {'<br>'.join(textwrap.wrap(LINES[i], 42))}</sup>"))) for i in range(6)]
        opts = {"frame": {"duration": 1400, "redraw": True}, "transition": {"duration": 0}, "fromcurrent": True}
        fig.update_layout(updatemenus=[dict(type="buttons", direction="left", x=0, y=-0.05, buttons=[
            dict(label="▶ Play", method="animate", args=[None, opts]),
            dict(label="Ⅱ Pause", method="animate", args=[[None], {"mode": "immediate", "frame": {"duration": 0, "redraw": False}}]),
            dict(label="↻ Replay", method="animate", args=[[str(i) for i in range(6)], {**opts, "fromcurrent": False, "mode": "immediate"}]),
        ])])
    return fig


def render():
    st.subheader("🎬 How Cas9 cuts · six moments")
    st.caption("Assumptions: simplified SpCas9/NGG cartoon; circles and lines stand for molecular parts. What this does not claim: atomic structure, timing, or guaranteed repair.")
    if st.session_state.get("lite_mode", False):
        for i in range(6):
            with st.expander(f"{i+1}. {STEPS[i]}", expanded=i == 0):
                st.plotly_chart(figure(i, False), config={"staticPlot": True}, key=f"cut_static_{i}")
    else:
        st.plotly_chart(figure(), config={"displayModeBar": False}, key="cut_animation")
