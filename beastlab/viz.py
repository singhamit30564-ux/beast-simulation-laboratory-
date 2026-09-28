"""Figures: Plotly charts plus hand-built SVG biology."""
from __future__ import annotations

import html
from typing import Any, Dict, List, Optional, Sequence

import plotly.graph_objects as go

from .theme import C

PLOT_LAYOUT = dict(
    paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
    font=dict(color=C["text"], family="Inter, sans-serif", size=12),
    margin=dict(l=50, r=20, t=48, b=44),
    legend=dict(bgcolor="rgba(0,0,0,0)", font=dict(size=11)),
    hoverlabel=dict(bgcolor=C["panel"], bordercolor=C["line"], font=dict(color=C["text"])),
)
GRID = dict(gridcolor=C["line"], zerolinecolor=C["line"])


def _fig(title: str = "", height: int = 380) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(**PLOT_LAYOUT, title=dict(text=title, font=dict(size=14, color=C["muted"])),
                      height=height)
    return fig


def _axes(fig: go.Figure, xtitle: str = "", ytitle: str = "", xlog: bool = False,
          ylog: bool = False) -> go.Figure:
    fig.update_xaxes(title=xtitle, **GRID, type="log" if xlog else "-")
    fig.update_yaxes(title=ytitle, **GRID, type="log" if ylog else "-")
    return fig


# --------------------------------------------------------------------------- #
# pathway / outcome charts
# --------------------------------------------------------------------------- #
def fig_pathway(percent: Dict[str, float], height: int = 330) -> go.Figure:
    labels = list(percent.keys())
    colors = {"c-NHEJ": C["plasma"], "MMEJ": C["violet"], "HDR": C["green"], "SSA": C["teal"]}
    fig = go.Figure(go.Pie(labels=labels, values=[percent[k] for k in labels], hole=0.62,
                           textinfo="label+percent", textfont=dict(size=11),
                           marker=dict(colors=[colors.get(k, C["muted"]) for k in labels],
                                       line=dict(color=C["bg"], width=2))))
    fig.update_layout(**PLOT_LAYOUT, height=height,
                      title=dict(text="Repair pathway competition at one DSB",
                                 font=dict(size=14, color=C["muted"])))
    return fig


def fig_indel_spectrum(hist: Dict[Any, int], height: int = 360) -> go.Figure:
    keys = sorted(int(k) for k in hist)
    vals = [hist[k] for k in keys] if all(k in hist for k in keys) else \
        [hist[str(k)] if str(k) in hist else hist[k] for k in keys]
    fig = _fig("Indel size spectrum (alleles per size)", height)
    colors = [C["green"] if k > 0 else C["plasma"] for k in keys]
    fig.add_trace(go.Bar(x=[str(k) for k in keys], y=vals, marker_color=colors,
                         name="alleles",
                         hovertemplate="Δ%{x} bp: %{y} alleles<extra></extra>"))
    fig.add_vline(x=len([k for k in keys if k < 0]) - 0.5, line=dict(color=C["line"], dash="dot"))
    fig.add_annotation(x=0.02, y=1.06, xref="paper", yref="paper", showarrow=False,
                       text="◀ deletions", font=dict(color=C["plasma"], size=11))
    fig.add_annotation(x=0.85, y=1.06, xref="paper", yref="paper", showarrow=False,
                       text="insertions ▶", font=dict(color=C["green"], size=11))
    return _axes(fig, "indel size (bp)", "alleles")


def fig_allele_pie(counts: Dict[str, int], title: str = "Allele composition",
                   height: int = 320) -> go.Figure:
    palette = [C["muted"], C["plasma"], C["violet"], C["green"], C["gold"], C["teal"]]
    fig = go.Figure(go.Pie(labels=list(counts.keys()), values=list(counts.values()), hole=0.55,
                           marker=dict(colors=palette[:len(counts)], line=dict(color=C["bg"], width=2)),
                           textinfo="label+percent"))
    fig.update_layout(**PLOT_LAYOUT, height=height,
                      title=dict(text=title, font=dict(size=14, color=C["muted"])))
    return fig


def fig_guide_scatter(guides: Sequence[Dict[str, Any]], height: int = 430) -> go.Figure:
    fig = _fig("Guide ranking", height)
    xs = [g["on_target_score"] for g in guides]
    ys = [g["n_offtargets"] for g in guides]
    texts = [f"{g['guide']} | {g['pam']}" for g in guides]
    colors = [C["green"] if g["n_offtargets"] == 0 else
              (C["gold"] if g["n_offtargets"] < 3 else C["plasma"]) for g in guides]
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", text=texts,
                             marker=dict(size=13, color=colors, line=dict(color=C["bg"], width=1)),
                             hovertemplate="%{text}<br>on-target score %{x}<br>"
                                           "off-target sites %{y}<extra></extra>"))
    return _axes(fig, "on-target score (heuristic)", "off-target sites (≤4 mismatches)")


def fig_bar(labels: Sequence[str], values: Sequence[float], title: str = "",
            color: str = C["teal"], height: int = 380, horizontal: bool = False,
            ytitle: str = "", xtitle: str = "") -> go.Figure:
    fig = _fig(title, height)
    if horizontal:
        fig.add_trace(go.Bar(y=list(labels), x=list(values), orientation="h",
                             marker_color=color))
        fig = _axes(fig, xtitle or "", "", )
    else:
        fig.add_trace(go.Bar(x=list(labels), y=list(values), marker_color=color))
        fig = _axes(fig, xtitle, ytitle)
    return fig


# --------------------------------------------------------------------------- #
# immunity charts
# --------------------------------------------------------------------------- #
def fig_arms_race(series: Dict[str, List[float]], events: Sequence[Dict[str, Any]] = (),
                  height: int = 460) -> go.Figure:
    fig = _fig("Phage × bacteria — population arms race", height)
    fig.add_trace(go.Scatter(x=series["t"], y=series["B_total"], name="bacteria (total)",
                             line=dict(color=C["teal"], width=3)))
    fig.add_trace(go.Scatter(x=series["t"], y=series["B1"], name="immune bacteria",
                             line=dict(color=C["green"], width=2, dash="dot")))
    fig.add_trace(go.Scatter(x=series["t"], y=series["P0"], name="phage (wild type)",
                             line=dict(color=C["plasma"], width=2)))
    if any(v > 0 for v in series.get("P1", [])):
        fig.add_trace(go.Scatter(x=series["t"], y=series["P1"], name="phage (escaped)",
                                 line=dict(color=C["gold"], width=2, dash="dash")))
    if any(v > 0 for v in series.get("P2", [])):
        fig.add_trace(go.Scatter(x=series["t"], y=series["P2"], name="phage (anti-CRISPR)",
                                 line=dict(color=C["violet"], width=2, dash="dot")))
    for ev in events:
        fig.add_vline(x=ev["t"], line=dict(color=C["red"], dash="dashdot", width=1))
        fig.add_annotation(x=ev["t"], y=1.02, yref="paper", text=ev["event"],
                           showarrow=False, font=dict(size=10, color=C["red"]),
                           textangle=0)
    fig.update_yaxes(type="log")
    return _axes(fig, "time (hours)", "count (log₁₀)", ylog=True)


def fig_moi(curve: Dict[str, float], height: int = 300) -> go.Figure:
    mo = curve["moi"]
    fig = _fig("Poisson infection probability", height)
    fig.add_trace(go.Bar(x=["sensitive cell<br>is infected", "immune cell<br>survives"],
                         y=[100 * curve["p_infected_sensitive"], 100 * curve["p_survive_if_immune"]],
                         marker_color=[C["plasma"], C["green"]]))
    return _axes(fig, "", "percent")


def fig_array_map(parse: Dict[str, Any], height: int = 210) -> go.Figure:
    fig = _fig("CRISPR array — repeat / spacer architecture", height)
    x, w, colors, texts = [], [], [], []
    for s in parse["spacer_objects"]:
        x.append(s["start"]); w.append(s["length"]); colors.append(C["violet"])
        texts.append(f"spacer<br>{s['sequence'][:40]}…<br>{s['length']} nt")
    repeat = parse["repeat"]
    for i in range(parse["repeats"]):
        start = parse["array_start"] + i * (len(repeat) + (len(parse["spacers"][0]) if parse["spacers"] else 0))
        x.append(start); w.append(len(repeat)); colors.append(C["gold"])
        texts.append(f"repeat<br>{repeat}")
    fig.add_trace(go.Bar(x=x, y=[1] * len(x), width=w, marker=dict(color=colors),
                         text=texts, hovertemplate="%{text}<br>starts at %{x}<extra></extra>"))
    fig.update_yaxes(visible=False)
    fig.update_xaxes(title="position in the locus (bp)")
    return fig


# --------------------------------------------------------------------------- #
# genome editing charts
# --------------------------------------------------------------------------- #
def fig_pam_coverage(rows: Sequence[Dict[str, Any]], height: int = 380) -> go.Figure:
    fig = _fig("PAM availability across the target region", height)
    fig.add_trace(go.Bar(x=[r["pam"] for r in rows], y=[r["pam_per_kb"] for r in rows],
                         marker_color=C["gold"], name="PAM density",
                         hovertemplate="%{x}: %{y:.1f} sites per kb<extra></extra>"))
    return _axes(fig, "effector PAM", "PAM sites per kb")


def fig_dose_response(x: Sequence[float], series: Dict[str, Sequence[float]],
                      title: str = "Dose–response", xlabel: str = "dose",
                      ylabel: str = "effect", height: int = 400) -> go.Figure:
    fig = _fig(title, height)
    palette = [C["gold"], C["teal"], C["plasma"], C["violet"], C["green"]]
    for i, (name, ys) in enumerate(series.items()):
        fig.add_trace(go.Scatter(x=list(x), y=list(ys), name=name, mode="lines+markers",
                                 line=dict(color=palette[i % len(palette)], width=2)))
    return _axes(fig, xlabel, ylabel)


def fig_gel(lanes: Sequence[Dict[str, Any]], height: int = 420) -> go.Figure:
    """Simulated agarose gel.

    Each lane is ``{"name": ..., "bands": [(size_bp, intensity), ...]}``.
    Migration is log-linear in fragment size.
    """
    fig = _fig("Genotyping gel (simulated)",
               height)
    fig.update_layout(xaxis=dict(showgrid=False), yaxis=dict(showgrid=False))
    max_bp = max([b[0] for lane in lanes for b in lane["bands"]] + [1000])
    min_bp = min([b[0] for lane in lanes for b in lane["bands"]] + [100])

    def y_of(size: int) -> float:
        import math
        return (math.log10(max_bp + 50) - math.log10(size + 50)) / \
               (math.log10(max_bp + 50) - math.log10(min_bp * 0.4 + 50))

    shapes, annotations = [], []
    for i, lane in enumerate(lanes):
        shapes.append(dict(type="rect", x0=i - 0.42, x1=i + 0.42, y0=0, y1=1,
                           fillcolor="rgba(20,28,45,0.85)", line=dict(color=C["line"])))
        for size, intensity in lane["bands"]:
            y = y_of(size)
            shapes.append(dict(type="rect", x0=i - 0.34, x1=i + 0.34, y0=y - 0.022, y1=y + 0.022,
                               fillcolor=f"rgba(212,175,55,{min(1.0, max(0.15, intensity))})",
                               line=dict(width=0)))
            annotations.append(dict(x=i + 0.40, y=y, text=f"{size} bp", showarrow=False,
                                    font=dict(size=9, color=C["muted"]), xanchor="left"))
        annotations.append(dict(x=i, y=1.02, text=lane["name"], showarrow=False,
                                font=dict(size=11, color=C["text"]), yanchor="bottom"))
    for size in (100, 200, 500, 1000, 2000, 5000):
        if min_bp * 0.5 <= size <= max_bp * 1.2:
            annotations.append(dict(x=-0.6, y=y_of(size), text=f"{size}", showarrow=False,
                                    font=dict(size=9, color=C["muted"])))
    fig.update_layout(shapes=shapes, annotations=annotations,
                      xaxis=dict(range=[-0.9, len(lanes) - 0.2], showticklabels=False),
                      yaxis=dict(range=[-0.02, 1.06], showticklabels=False,
                                 title="migration →"))
    return fig


def fig_radar(values: Dict[str, float], title: str = "Safety profile", height: int = 400) -> go.Figure:
    cats = list(values.keys())
    vals = [values[c] for c in cats]
    fig = go.Figure(go.Scatterpolar(r=vals + [vals[0]], theta=cats + [cats[0]],
                                    fill="toself", fillcolor="rgba(255,122,89,0.25)",
                                    line=dict(color=C["plasma"], width=2)))
    fig.update_layout(**PLOT_LAYOUT, height=height,
                      polar=dict(bgcolor="rgba(0,0,0,0)",
                                 radialaxis=dict(gridcolor=C["line"], range=[0, 1],
                                                 tickfont=dict(size=9, color=C["muted"])),
                                 angularaxis=dict(gridcolor=C["line"],
                                                  tickfont=dict(size=10, color=C["text"]))),
                      title=dict(text=title, font=dict(size=14, color=C["muted"])))
    return fig


def fig_histogram(values: Sequence[float], title: str = "", xlabel: str = "",
                  height: int = 340, color: str = C["teal"]) -> go.Figure:
    fig = _fig(title, height)
    fig.add_trace(go.Histogram(x=list(values), marker_color=color, nbinsx=24))
    return _axes(fig, xlabel, "patients")


def fig_scatter(x: Sequence[float], y: Sequence[float], title: str = "",
                xlabel: str = "", ylabel: str = "", color: str = C["gold"],
                labels: Optional[Sequence[str]] = None, height: int = 400) -> go.Figure:
    fig = _fig(title, height)
    fig.add_trace(go.Scatter(x=list(x), y=list(y), mode="markers",
                             text=list(labels) if labels else None,
                             marker=dict(size=10, color=color, line=dict(color=C["bg"], width=1))))
    return _axes(fig, xlabel, ylabel)


# --------------------------------------------------------------------------- #
# SVG biology
# --------------------------------------------------------------------------- #
def svg_bacterium(width: int = 460, immune: bool = False, array_spacers: int = 7,
                  new_spacers: int = 0, label: str = "") -> str:
    """Rod-shaped bacterium with nucleoid, CRISPR array and (optional) phage."""
    spacer_w = max(6, int(300 / max(1, array_spacers + 1)))
    ladder = []
    x = 90
    ladder.append(f'<rect x="{x}" y="86" width="10" height="18" fill="{C["green"]}" rx="2"/>')
    x += 14
    for i in range(array_spacers + new_spacers):
        colour = C["gold"] if i % 2 == 0 else C["violet"]
        if new_spacers and i < new_spacers:
            colour = C["green"]
        ladder.append(f'<rect x="{x}" y="86" width="{spacer_w}" height="18" fill="{colour}" '
                      f'rx="2" opacity="0.92"/>')
        x += spacer_w + 3
    ring = C["green"] if immune else C["muted"]
    return f"""
<svg viewBox="0 0 460 190" width="100%" height="{int(width*0.41)}" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <linearGradient id="cellg" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#16324a" stop-opacity="0.95"/>
      <stop offset="100%" stop-color="#0d1b2c" stop-opacity="0.95"/>
    </linearGradient>
  </defs>
  <rect x="20" y="20" width="420" height="150" rx="72" fill="url(#cellg)" stroke="{ring}" stroke-width="2.5"/>
  <ellipse cx="230" cy="95" rx="150" ry="46" fill="#0b1524" opacity="0.75"/>
  <ellipse cx="230" cy="95" rx="150" ry="46" fill="none" stroke="{C['line']}" stroke-dasharray="4 5"/>
  {' '.join(ladder)}
  <text x="230" y="128" text-anchor="middle" fill="{C['muted']}" font-size="11"
        font-family="Inter, sans-serif">nucleoid · CRISPR array {array_spacers + new_spacers} spacers</text>
  <text x="230" y="160" text-anchor="middle" fill="{C['text']}" font-size="12" font-weight="700"
        font-family="Inter, sans-serif">{html.escape(label)}</text>
  {''.join(f'<circle cx="{52 + i*22}" cy="{38 + (i%2)*8}" r="4" fill="{C["teal"]}" opacity="0.55"/>' for i in range(16))}
  <text x="36" y="24" fill="{C['muted']}" font-size="10" font-family="Inter, sans-serif">ribosomes</text>
</svg>"""


def svg_phage(width: int = 220, acr: bool = False) -> str:
    return f"""
<svg viewBox="0 0 220 190" width="100%" height="{int(width*0.86)}" xmlns="http://www.w3.org/2000/svg">
  <polygon points="110,18 152,44 152,92 110,118 68,92 68,44" fill="#233252"
           stroke="{C['plasma']}" stroke-width="2.4"/>
  <polygon points="110,30 140,50 140,84 110,104 80,84 80,50" fill="#16233c"/>
  <rect x="104" y="118" width="12" height="42" fill="{C['muted']}"/>
  <path d="M84 168 L104 160 M136 168 L116 160 M110 172 L110 160" stroke="{C['violet']}"
        stroke-width="3" fill="none"/>
  <text x="110" y="146" text-anchor="middle" fill="{C['muted']}" font-size="10"
        font-family="Inter, sans-serif">tail fibres</text>
  {f'<text x="110" y="14" text-anchor="middle" fill="{C["gold"]}" font-size="11" font-weight="700">+ anti-CRISPR genes</text>' if acr else ''}
</svg>"""


def svg_human_cell(width: int = 420, edited_percent: int = 0, label: str = "") -> str:
    return f"""
<svg viewBox="0 0 420 210" width="100%" height="{int(width*0.5)}" xmlns="http://www.w3.org/2000/svg">
  <ellipse cx="210" cy="105" rx="200" ry="98" fill="#0e1a2b" stroke="{C['line']}" stroke-width="2"/>
  <ellipse cx="210" cy="100" rx="118" ry="66" fill="#16233c" stroke="{C['violet']}" stroke-width="1.6"/>
  {''.join(f'<path d="M{150+i*17} {70+ (i%3)*22} q10 -14 22 0 q-10 16 -22 0" stroke="{C["muted"]}" fill="none" stroke-width="1.1"/>' for i in range(8))}
  <text x="210" y="48" text-anchor="middle" fill="{C['muted']}" font-size="11"
        font-family="Inter, sans-serif">nucleus · 46 chromosomes</text>
  <text x="210" y="196" text-anchor="middle" fill="{C['text']}" font-size="12" font-weight="700"
        font-family="Inter, sans-serif">{html.escape(label)}</text>
  <rect x="120" y="176" width="180" height="8" rx="4" fill="#1b2740"/>
  <rect x="120" y="176" width="{int(1.8*edited_percent)}" height="8" rx="4" fill="{C['green']}"/>
  <text x="308" y="184" fill="{C['muted']}" font-size="10" font-family="Inter, sans-serif">
     {edited_percent}% alleles edited</text>
</svg>"""


def svg_cas_domains(enzyme: Dict[str, Any], width: int = 700) -> str:
    """Domain map for Cas9-type effectors with the cut position marked."""
    fam = enzyme.get("family", "Cas9")
    if "Cas12" in fam:
        doms = [("RuvC", 0.03, 0.22, C["plasma"]), ("WED", 0.22, 0.38, C["teal"]),
                ("PI", 0.38, 0.52, C["violet"]), ("REC", 0.52, 0.72, C["blue"]),
                ("Nuc", 0.72, 0.97, C["gold"])]
    elif "Cas13" in fam:
        doms = [("HEPN1", 0.05, 0.28, C["plasma"]), ("HEPN2", 0.30, 0.55, C["plasma"]),
                ("REC", 0.55, 0.78, C["blue"]), ("crRNA-binding", 0.78, 0.97, C["teal"])]
    else:
        doms = [("RuvC-I", 0.01, 0.06, C["plasma"]), ("REC1", 0.06, 0.16, C["teal"]),
                ("Bridge helix", 0.16, 0.21, C["gold"]), ("REC2", 0.21, 0.30, C["blue"]),
                ("REC3", 0.30, 0.55, C["blue"]), ("HNH", 0.55, 0.69, C["plasma"]),
                ("RuvC-II", 0.69, 0.76, C["plasma"]), ("RuvC-III", 0.76, 0.84, C["plasma"]),
                ("PI (PAM readout)", 0.84, 0.99, C["violet"])]
    blocks = []
    for name, a, b, colour in doms:
        blocks.append(
            f'<rect x="{20 + a*660:.1f}" y="34" width="{(b-a)*660:.1f}" height="34" rx="6" '
            f'fill="{colour}" opacity="0.85"/>'
            f'<text x="{20 + (a+ (b-a)/2)*660:.1f}" y="56" text-anchor="middle" '
            f'font-size="10" fill="#0b0f18" font-family="Inter, sans-serif" font-weight="700">'
            f'{html.escape(name)}</text>')
    cut = enzyme.get("nts_cut", 17)
    x_cut = 20 + 660 * 0.5 + (0 if enzyme.get("pam_side") == "5" else 60)
    return f"""
<svg viewBox="0 0 700 120" width="100%" height="{int(width*0.17)}" xmlns="http://www.w3.org/2000/svg">
  <text x="20" y="18" fill="{C['muted']}" font-size="11" font-family="Inter, sans-serif">
    {html.escape(enzyme['name'])} — {enzyme['size_aa']} aa · {html.escape(enzyme['system'])} ·
    PAM {html.escape(enzyme['pam'])} ({enzyme['pam_side']}') · spacer {enzyme['spacer_len']} nt
  </text>
  {''.join(blocks)}
  <line x1="{x_cut}" y1="24" x2="{x_cut}" y2="86" stroke="{C['red']}" stroke-width="2" stroke-dasharray="4 4"/>
  <text x="{x_cut}" y="100" text-anchor="middle" fill="{C['red']}" font-size="11"
        font-family="Inter, sans-serif">catalytic cut · {enzyme['cut_geometry']} · Δ{enzyme.get('overhang',0)} nt</text>
</svg>"""


def sequence_diff_html(original: str, edited: str, context: int = 26) -> str:
    """Show the first differing region of two alleles, with the difference painted."""
    import difflib

    if original == edited:
        return ('<div class="seqbox">No difference — the edited allele is identical to '
                'the reference.</div>')
    matcher = difflib.SequenceMatcher(None, original, edited, autojunk=False)
    diffs = [op for op in matcher.get_opcodes() if op[0] != "equal"]
    tag, i1, i2, j1, j2 = diffs[0]
    lo = max(0, i1 - context)

    def paint(seq: str, a: int, b: int, pre_start: int) -> str:
        pre = seq[pre_start:a]
        core = seq[a:b]
        post = seq[b:min(len(seq), b + context)]
        core_html = "".join(
            f'<span style="background:{C["gold"]};color:#0b0f18;font-weight:700">'
            f'{html.escape(ch)}</span>' for ch in core)
        return (f'<span style="color:{C["muted"]}">{html.escape(pre)}</span>' + core_html
                + f'<span style="color:{C["muted"]}">{html.escape(post)}</span>')

    ref_line = paint(original, i1, i2, min(lo, i1))
    alt_line = paint(edited, j1, j2, min(lo, j1))
    kind = {"replace": "substitution", "delete": "deletion",
            "insert": "insertion"}.get(tag, tag)
    return ('<div class="seqbox">'
            f'<div style="color:{C["muted"]};font-size:0.72rem;margin-bottom:4px">'
            f'{kind} · reference positions {i1 + 1}-{i2} ({i2 - i1} nt) → '
            f'edited positions {j1 + 1}-{j2} ({j2 - j1} nt)</div>'
            f'<div><span class="ruler">ref</span>  {ref_line}</div>'
            f'<div><span class="ruler">edit</span>  {alt_line}</div>'
            "</div>")
