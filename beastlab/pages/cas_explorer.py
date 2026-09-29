"""Cas explorer — the effector encyclopedia and PAM-coverage bench."""
from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd
import streamlit as st
from ..presentation import table, chart

from .. import theme, viz
from ..crispr import build_guide_rna, pam_coverage, scan_pams
from ..data.cas_enzymes import (BASE_EDITORS, EFFECTOR_BY_ID, NUCLEASES,
                                PAM_COVERAGE_SETS, POINT_MUTATION_CLASSES,
                                PRIME_EDITORS, RNA_EFFECTORS)
from ..data.targets import all_loci, get_locus

ALL = NUCLEASES + RNA_EFFECTORS

CLASS_DIAGRAM = f"""
<div class="card teal">
<h4>Class 1 — multi-protein machines (the ancestors)</h4>
<p><b>Type I</b> — Cascade (Cas8/Cas7/Cas5/Cas6) recruits <b>Cas3</b>, a helicase-nuclease that
chews the target processively. PAM 5'-AAG-3' (E. coli I-E) or 5'-CC-3' (I-F). <i>E. coli,
P. aeruginosa, T. thermophilus.</i><br>
<b>Type III</b> — Csm/Cmr complexes bind the <i>transcript</i> of a target and cut
transcriptionally active DNA, plus cyclic-oligonucleotide signalling. No strict PAM.
<i>S. thermophilus CRISPR2, M. tuberculosis.</i><br>
<b>Type IV</b> — no nuclease domain at all; partnerships with other systems.</p>
</div>
<div class="card violet">
<h4>Class 2 — single-protein effectors (the tools)</h4>
<p><b>Type II</b> — Cas9 + tracrRNA + crRNA. Blunt cut 3 bp upstream of the PAM.
<i>SpCas9, SaCas9, Nme2Cas9, CjCas9 …</i><br>
<b>Type V</b> — Cas12a/b/f/i/j. Single crRNA, staggered cut, T-rich PAMs, some process their own
arrays. <i>AsCas12a, LbCas12a, Un1Cas12f1 …</i><br>
<b>Type VI</b> — Cas13a/b/d. RNA targets only, collateral RNase activity (SHERLOCK/DETECTR).</p>
</div>
"""


def _effector_table(effectors: List[Dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame([{
        "id": e["id"], "name": e["name"], "family": e.get("family", ""),
        "species": e.get("species", ""), "type": e.get("system", ""),
        "size (aa)": e.get("size_aa", ""),
        "PAM": e.get("pam", ""), "PAM side": (e.get("pam_side", "") + "'")
            if e.get("pam_side") != "none" else "—",
        "spacer": e.get("spacer_len", ""), "cut": e.get("cut_geometry", ""),
        "seed": e.get("seed_len", ""),
        "tracrRNA": "yes" if e.get("tracrrna") else "no",
        "collateral": "yes" if e.get("collateral") else "no",
        "temp": e.get("temperature", ""),
    } for e in effectors])


def render() -> None:
    st.markdown(theme.hero(
        "🔬 Cas explorer",
        "Every effector in the toolbox: what it recognises, how big it is, where it cuts, and "
        "what it does that its relatives cannot.",
        ["30+ effectors", "editors included", "PAM coverage bench"],
    ), unsafe_allow_html=True)

    tab_eff, tab_compare, tab_pam, tab_editors, tab_class = st.tabs(
        ["Effector library", "Compare side by side", "PAM coverage bench",
         "Base & prime editors", "Classification"])

    # ---------------------------------------------------------------- library
    with tab_eff:
        c1, c2, c3, c4 = st.columns([1.3, 1, 1, 1.2])
        with c1:
            families = sorted({e.get("family", "").split(" (")[0] for e in ALL})
            fam = st.multiselect("Family", families, default=[])
        with c2:
            min_aa, max_aa = st.slider("Size (aa)", 400, 1700, (400, 1700), step=25)
        with c3:
            side = st.multiselect("PAM side", ["3", "5", "none"], default=[])
        with c4:
            query = st.text_input("Search", placeholder="species, PAM, id…")

        rows = []
        for e in ALL:
            if fam and e.get("family", "").split(" (")[0] not in fam:
                continue
            if not (min_aa <= e.get("size_aa", 0) <= max_aa):
                continue
            if side and e.get("pam_side") not in side:
                continue
            hay = " ".join(str(e.get(k, "")) for k in ("id", "name", "species", "pam", "notes")).lower()
            if query and query.lower() not in hay:
                continue
            rows.append(e)
        st.caption(f"{len(rows)} effectors match.")
        table(_effector_table(rows), width="stretch", hide_index=True, height=420)

        st.markdown("### Effector detail")
        chosen = st.selectbox("Effector", [e["id"] for e in ALL],
                              format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
        e = EFFECTOR_BY_ID[chosen]
        st.markdown(viz.svg_cas_domains(e), unsafe_allow_html=True)
        a, b = st.columns([1.15, 1])
        with a:
            st.markdown(theme.kv({
                "species": e.get("species", ""), "system": e.get("system", ""),
                "size": f"{e.get('size_aa','?')} aa", "PAM": e.get("pam", "—"),
                "PAM side": e.get("pam_side"), "spacer": f"{e.get('spacer_len')} nt",
                "cut": f"{e.get('cut_geometry')} · Δ{e.get('overhang', 0)} nt",
                "seed": f"{e.get('seed_len')} nt", "tracrRNA": e.get("tracrrna"),
                "collateral": e.get("collateral"), "temperature": e.get("temperature"),
            }), unsafe_allow_html=True)
            st.markdown(theme.card("Assessment", e.get("notes", ""), "teal"), unsafe_allow_html=True)
            st.caption(f"Reference: {e.get('ref','—')}")
        with b:
            ps = st.text_input("Guide spacer (protospacer) to build the RNA for",
                               value="GTCACCAATCCTGTCCCTAG", max_chars=40,
                               help="Any 20-28 nt sequence; this is only to render the RNA.")
            g = build_guide_rna(ps, e)
            st.markdown(theme.card(f"{g['name']} · {g['len']} nt",
                                   f"<span class='mono' style='word-break:break-all'>{g['rna']}</span><br>"
                                   f"<span style='color:{theme.C['muted']}'>{g['note']}</span>", "violet"),
                        unsafe_allow_html=True)
            if g.get("scaffold"):
                st.caption("Scaffold: " + g["scaffold"])
            st.markdown(theme.callout(
                "Scaffolds differ between orthologs — swapping a Cas9 without swapping its "
                "guide architecture is one of the commonest reasons an experiment simply "
                "does not work.", "info"), unsafe_allow_html=True)

    # --------------------------------------------------------------- compare
    with tab_compare:
        picks = st.multiselect("Pick up to 4 effectors to compare",
                               [e["id"] for e in ALL], default=["SpCas9", "SaCas9",
                                                                "AsCas12a", "Un1Cas12f1"])
        if picks:
            df = _effector_table([EFFECTOR_BY_ID[p] for p in picks[:4]])
            # transpose so each column is one effector; mixed dtypes must be cast to
            # text before they reach Arrow (int + str in one column otherwise fails)
            cmp_df = df.set_index("name").T.astype(str)
            table(cmp_df, width="stretch")
            st.markdown(theme.callout(
                "The trade-off triangle: <b>size</b> (delivery), <b>PAM</b> (where you can cut) "
                "and <b>cut geometry</b> (blunt vs staggered — which decides whether the target "
                "survives to be re-cut).", "info"), unsafe_allow_html=True)

    # ------------------------------------------------------------ PAM bench
    with tab_pam:
        st.caption("How much of a real sequence does each PAM actually reach? This is computed "
                   "on the shipped genome library.")
        c1, c2, c3 = st.columns([1.2, 1.4, 1])
        with c1:
            locus_id = st.selectbox("Target", [l["id"] for l in all_loci()],
                                    format_func=lambda i: get_locus(i)["name"])
        with c2:
            pams = st.multiselect("PAM sets", list(PAM_COVERAGE_SETS.keys()),
                                  default=["SpCas9 (NGG)", "SpCas9-NG (NG)", "SpRY (NRN)",
                                           "SpRY (NYN)"])
        with c3:
            window = st.slider("Coverage window (bp)", 10, 200, 50, step=10)

        locus = get_locus(locus_id)
        if pams and locus:
            rows = []
            for name in pams:
                cov = pam_coverage(locus["sequence"], PAM_COVERAGE_SETS[name], window=window)
                cov["pam"] = name
                rows.append(cov)
            chart(viz.fig_pam_coverage(rows), width="stretch")
            table(pd.DataFrame([{
                "PAM set": r["pam"], "sites in locus": r["n_pam_sites"],
                "per kb": round(r["pam_per_kb"], 1),
                "median gap (bp)": r["median_distance_bp"],
                f"coverage within {window} bp": f"{100*r[f'coverage_within_{window}_bp']:.1f}%",
            } for r in rows]), width="stretch", hide_index=True)
            st.markdown(theme.callout(
                f"Same {locus['length']:,} bp of {locus['name']}. NGG sites are frequent but not "
                "uniform; near-PAMless variants reach almost everything — and pay for it with a "
                "bigger off-target search space (see Guide design).", "info"),
                unsafe_allow_html=True)
            with st.expander("Every PAM site in this locus"):
                site_rows = []
                for name in pams:
                    sites = scan_pams(locus["sequence"], {
                        "id": name, "pam": PAM_COVERAGE_SETS[name],
                        "pam_side": "3", "spacer_len": 20, "nts_cut": 17, "ts_cut": 17,
                        "cut_geometry": "blunt"})
                    for s in sites[:80]:
                        site_rows.append({"PAM set": name, "strand": s["strand"],
                                          "position": s["start_1"],
                                          "protospacer": s["protospacer"],
                                          "PAM": s["pam"],
                                          "cut at": s["cut_site_1"]})
                table(pd.DataFrame(site_rows), width="stretch", hide_index=True)

    # --------------------------------------------------------------- editors
    with tab_editors:
        st.markdown("#### Base editors — chemistry instead of scissors")
        table(pd.DataFrame([{
            "id": b["id"], "family": b["family"], "conversion": b["conversion"],
            "Cas domain": b["cas"], "PAM": b["pam"],
            "window (protospacer positions)": f"{b['window'][0]}–{b['window'][1]}",
            "motif": b["preferred_motif"], "bystander risk": b["bystander"],
            "indel rate": f"{100*b['indel_rate']:.1f}%",
            "max efficiency": f"{100*b['max_efficiency']:.0f}%",
        } for b in BASE_EDITORS]), width="stretch", hide_index=True)
        st.markdown(theme.callout(
            "Cytosine editors (CBE) deaminate C→U, which replication reads as T. Adenine editors "
            "(ABE) deaminate A→inosine, read as G. The window is a property of the protein "
            "architecture: ABE8e is fast but wide (bystander risk), ABE9 is narrow and precise.",
            "info"), unsafe_allow_html=True)

        st.markdown("#### Prime editors — search and replace")
        table(pd.DataFrame([{
            "id": p["id"], "cas/RT": p["cas"], "PAM": p["pam"],
            "max efficiency": f"{100*p['max_efficiency']:.0f}%",
            "indel rate": f"{100*p['indel_rate']:.1f}%", "notes": p["notes"],
        } for p in PRIME_EDITORS]), width="stretch", hide_index=True)
        st.markdown("#### The 12 point mutations — and which editor can make them")
        table(pd.DataFrame(POINT_MUTATION_CLASSES), width="stretch", hide_index=True)

    # ---------------------------------------------------------- class tree
    with tab_class:
        st.markdown(CLASS_DIAGRAM, unsafe_allow_html=True)
        counts = {
            "Type II (Cas9)": sum(1 for e in NUCLEASES if e["system"].startswith("II")),
            "Type V (Cas12)": sum(1 for e in NUCLEASES if e["system"].startswith("V")),
            "Type VI (Cas13)": len(RNA_EFFECTORS),
        }
        chart(viz.fig_bar(list(counts), list(counts.values()),
                                    "Effectors in this database by type",
                                    theme.C["gold"], 300), width="stretch")
        st.caption("Class 1 systems (I, III, IV) are multi-protein machines that we do not "
                   "repackage as tools; their biology lives in the Immunity arena where it belongs.")
