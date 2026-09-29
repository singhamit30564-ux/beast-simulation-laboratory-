"""Visual identity: palette, CSS, cards and HTML sequence rendering."""
from __future__ import annotations

import html
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

#: Core palette used by every figure and card in the app.
C = {
    "bg": "#080b12",
    "panel": "#111827",
    "panel2": "#0e1524",
    "line": "#1f2b44",
    "text": "#e8eefc",
    "muted": "#93a3bd",
    "gold": "#d4af37",
    "teal": "#3fb6a8",
    "plasma": "#ff7a59",
    "violet": "#7f8cff",
    "green": "#48d597",
    "red": "#ff5470",
    "blue": "#4aa8ff",
}

BASE_COLORS = {"A": "#4aa8ff", "C": "#48d597", "G": "#ffb454", "T": "#ff7a8a", "N": "#93a3bd"}

CSS = f"""
<style>
  html, body, [class*="css"] {{ font-family: 'Inter', system-ui, sans-serif; }}
  .stApp {{ background:
      radial-gradient(1100px 600px at 12% -10%, rgba(212,175,55,0.10), transparent 60%),
      radial-gradient(900px 500px at 95% 0%, rgba(63,182,168,0.10), transparent 55%),
      {C['bg']}; }}
  section[data-testid="stSidebar"] {{ background: {C['panel2']}; border-right: 1px solid {C['line']}; }}
  h1, h2, h3 {{ letter-spacing: -0.01em; }}
  code, .mono {{ font-family: 'JetBrains Mono', ui-monospace, monospace; }}

  .beast-hero {{
    border: 1px solid {C['line']};
    background: linear-gradient(120deg, rgba(212,175,55,0.10), rgba(127,140,255,0.07) 45%, rgba(63,182,168,0.10));
    border-radius: 18px; padding: 22px 26px; margin-bottom: 14px;
  }}
  .beast-hero h1 {{ margin: 0 0 6px 0; font-size: 2.05rem; color: {C['text']}; }}
  .beast-hero p {{ margin: 0; color: {C['muted']}; font-size: 1.0rem; }}
  .beast-badges {{ margin-top: 12px; }}
  .badge {{
    display: inline-block; padding: 3px 11px; border-radius: 999px; font-size: 0.74rem;
    font-weight: 700; letter-spacing: 0.04em; text-transform: uppercase;
    margin-right: 7px; margin-bottom: 5px; border: 1px solid transparent;
  }}
  .card {{
    border: 1px solid {C['line']}; border-left: 3px solid {C['gold']};
    background: {C['panel']}; border-radius: 12px; padding: 14px 16px; margin-bottom: 10px;
  }}
  .card h4 {{ margin: 0 0 6px 0; font-size: 0.98rem; color: {C['text']}; }}
  .card p {{ margin: 0; color: {C['muted']}; font-size: 0.87rem; line-height: 1.45; }}
  .card.teal {{ border-left-color: {C['teal']}; }}
  .card.violet {{ border-left-color: {C['violet']}; }}
  .card.plasma {{ border-left-color: {C['plasma']}; }}
  .card.green {{ border-left-color: {C['green']}; }}
  .card.red {{ border-left-color: {C['red']}; }}
  .kv {{ display: flex; flex-wrap: wrap; gap: 8px; }}
  .kv div {{
    background: {C['panel2']}; border: 1px solid {C['line']}; border-radius: 9px;
    padding: 7px 11px; font-size: 0.82rem; color: {C['muted']};
  }}
  .kv div b {{ color: {C['text']}; }}
  .seqbox {{
    font-family: 'JetBrains Mono', monospace; font-size: 0.82rem; line-height: 1.65;
    background: #0a0f1a; border: 1px solid {C['line']}; border-radius: 10px;
    padding: 12px 14px; overflow-x: auto; white-space: pre;
  }}
  .ruler {{ color: {C['muted']}; font-size: 0.7rem; }}
  .legend span {{ margin-right: 12px; font-size: 0.78rem; color: {C['muted']}; }}
  .callout {{
    border-radius: 10px; padding: 11px 14px; font-size: 0.87rem; margin: 8px 0;
    border: 1px solid {C['line']}; background: {C['panel2']}; color: {C['muted']};
  }}
  .callout.warn {{ border-left: 3px solid {C['plasma']}; }}
  .callout.good {{ border-left: 3px solid {C['green']}; }}
  .callout.info {{ border-left: 3px solid {C['blue']}; }}
  .metricrow {{ display: flex; gap: 10px; flex-wrap: wrap; }}
  .metric {{
    flex: 1 1 140px; background: {C['panel']}; border: 1px solid {C['line']};
    border-radius: 12px; padding: 11px 13px;
  }}
  .metric .l {{ font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.06em; color: {C['muted']}; }}
  .metric .v {{ font-size: 1.5rem; font-weight: 800; color: {C['text']}; }}
  .metric .s {{ font-size: 0.72rem; color: {C['muted']}; }}
  .arrayrow {{ display: flex; flex-wrap: wrap; align-items: center; gap: 3px; margin: 6px 0 12px 0; }}
  .repeat-block {{
    background: repeating-linear-gradient(45deg, #d4af37, #d4af37 4px, #b3922d 4px, #b3922d 8px);
    color: #12161f; border-radius: 4px; font-size: 0.6rem; font-weight: 800;
    padding: 5px 4px; text-align: center; min-width: 16px;
  }}
  .spacer-block {{
    background: linear-gradient(180deg, #1c2a44, #16213a); border: 1px solid #2c3f61;
    color: {C['text']}; border-radius: 4px; font-size: 0.62rem; padding: 5px 6px;
    font-family: 'JetBrains Mono', monospace;
  }}
  .spacer-block.new {{ border-color: {C['green']}; box-shadow: 0 0 0 1px rgba(72,213,151,0.35); }}
  @media (max-width: 640px) {{
    .beast-hero {{ padding: 16px; }}
    .beast-hero h1 {{ font-size: 1.55rem; }}
    [data-testid="stHorizontalBlock"] {{ flex-wrap: wrap; }}
    [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {{
      width: 100% !important; flex: 1 1 100% !important; min-width: 0 !important;
    }}
    .stButton button, .stDownloadButton button {{ min-height: 44px; }}
    .card, .callout {{ overflow-wrap: anywhere; }}
  }}
</style>
"""


def badge(text: str, color: Optional[str] = None, bg: Optional[str] = None) -> str:
    color = color or C["gold"]
    bg = bg or "rgba(212,175,55,0.14)"
    return f'<span class="badge" style="color:{color};background:{bg};border-color:{color}33">{html.escape(text)}</span>'


def hero(title: str, subtitle: str, badges: Iterable[str] = ()) -> str:
    b = "".join(badge(x) for x in badges)
    return (f'<div class="beast-hero"><h1>{title}</h1><p>{subtitle}</p>'
            f'<div class="beast-badges">{b}</div></div>')


def card(title: str, body: str, kind: str = "") -> str:
    return (f'<div class="card {kind}"><h4>{title}</h4><p>{body}</p></div>')


def kv(pairs: Dict[str, Any]) -> str:
    items = "".join(f"<div><b>{html.escape(str(k))}</b> · {html.escape(str(v))}</div>"
                    for k, v in pairs.items())
    return f'<div class="kv">{items}</div>'


def callout(text: str, kind: str = "info") -> str:
    return f'<div class="callout {kind}">{text}</div>'


def metrics(items: Sequence[Tuple[str, str, str]]) -> str:
    """``items`` = [(label, value, sub), ...]"""
    inner = "".join(f'<div class="metric"><div class="l">{html.escape(l)}</div>'
                    f'<div class="v">{v}</div><div class="s">{html.escape(s)}</div></div>'
                    for l, v, s in items)
    return f'<div class="metricrow">{inner}</div>'


# --------------------------------------------------------------------------- #
# sequence rendering
# --------------------------------------------------------------------------- #
def render_sequence(seq: str, start: int = 1, width: int = 60,
                    highlights: Sequence[Tuple[int, int, str, str]] = (),
                    wrap: bool = True, show_ruler: bool = True) -> str:
    """Colourised, monospace sequence block.

    ``highlights`` entries are ``(start, end, colour, label)`` using 1-based
    inclusive coordinates in the same frame as ``seq``.
    """
    seq = seq.upper()
    spans = sorted(highlights, key=lambda h: h[0])
    lines: List[str] = []
    ruler_lines: List[str] = []

    def wrap_line(chunk_start: int, text: str, chunk_len: int) -> str:
        out: List[str] = []
        for offset, base in enumerate(text):
            pos = chunk_start + offset  # 1-based
            colour = BASE_COLORS.get(base, C["muted"])
            style = f"color:{colour}"
            for (hs, he, hc, _lab) in spans:
                if hs <= pos <= he:
                    style = f"color:#0b0f18;background:{hc};border-radius:2px;font-weight:700"
                    break
            out.append(f'<span style="{style}">{base}</span>')
        idx = "".join(str((chunk_start + k - 1) % 10) for k in range(chunk_len))
        ruler_lines.append(f'<span class="ruler">{chunk_start:>6}</span>  <span class="ruler">{idx}</span>')
        return "".join(out)

    if wrap:
        for i in range(0, len(seq), width):
            chunk = seq[i:i + width]
            lines.append(f'<span class="ruler">{start + i:>6}</span>  '
                         + wrap_line(start + i, chunk, len(chunk)))
    else:
        lines.append(f'<span class="ruler">{start:>6}</span>  ' + wrap_line(start, seq, len(seq)))
    body = "\n".join(lines)
    legend = ""
    if spans:
        legend = '<div class="legend" style="margin-bottom:6px">' + "".join(
            f'<span><span style="display:inline-block;width:10px;height:10px;background:{hc};'
            f'border-radius:2px;margin-right:4px"></span>{html.escape(lab)} '
            f'({hs}–{he})</span>' for hs, he, hc, lab in spans) + "</div>"
    return f'<div class="seqbox">{legend}{body}</div>'


def crispr_array_html(repeat: str, spacers: Sequence[str], new_count: int = 0,
                      leader: str = "") -> str:
    """Draw a CRISPR array as leader + repeat/spacer ladder."""
    out = ['<div class="arrayrow">']
    if leader:
        out.append(f'<div class="spacer-block" style="background:#0f1a12;border-color:{C["green"]}">'
                   f'LEADER ({len(leader)} bp)</div>')
    out.append('<div class="repeat-block" title="direct repeat">R</div>')
    for i, sp in enumerate(spacers):
        cls = "spacer-block new" if i < new_count else "spacer-block"
        label = f"S{len(spacers) - i}"
        out.append(f'<div class="{cls}" title="{html.escape(sp)}">{label}·{len(sp)}nt</div>')
        out.append('<div class="repeat-block" title="direct repeat">R</div>')
    out.append("</div>")
    return "".join(out)


def genome_browser_html(seq: str, features: Sequence[Dict[str, Any]], width: int = 900) -> str:
    """Horizontal track showing feature blocks over a genome bar."""
    n = max(1, len(seq))
    bars = []
    for f in features:
        left = 100 * f["start"] / n
        w = max(0.4, 100 * (f["end"] - f["start"]) / n)
        bars.append(
            f'<div title="{html.escape(str(f.get("label","")))}" style="position:absolute;'
            f'left:{left:.3f}%;width:{w:.3f}%;height:{f.get("height", 16)}px;'
            f'top:{f.get("top", 4)}px;background:{f.get("color", C["gold"])};'
            f'border-radius:3px;opacity:0.9"></div>')
    return (f'<div style="position:relative;width:100%;height:44px;background:{C["panel2"]};'
            f'border:1px solid {C["line"]};border-radius:8px;overflow:hidden">'
            + "".join(bars) + "</div>")
