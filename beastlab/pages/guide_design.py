"""Guide design bench — real PAM scans, ranked guides, off-target honesty."""
from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd
import streamlit as st

from .. import theme, viz
from ..crispr import (build_guide_rna, design_guides, find_offtargets,
                      pam_coverage, protospacer_for_guide, recommend_pbs, scan_pams)
from ..data.cas_enzymes import EFFECTOR_BY_ID, NUCLEASES, PRIME_EDITORS
from ..data.targets import all_loci, coding_strand, get_locus
from ..sequtils import clean, gc_content


def _guide_frame(guides: List[Dict[str, Any]]) -> pd.DataFrame:
    return pd.DataFrame([{
        "rank": i + 1,
        "guide": g["guide"],
        "PAM": g["pam"],
        "strand": g["strand"],
        "cut at": g["cut_site_1"],
        "GC%": round(100 * g["gc"], 1),
        "on-target": g["on_target_score"],
        "specificity": g["specificity_score"],
        "off-targets": g["n_offtargets"],
        "seed-perfect off-targets": g["n_seed_perfect_offtargets"],
        "warnings": "; ".join(g["warnings"]) or "—",
    } for i, g in enumerate(guides)])


def render() -> None:
    st.markdown(theme.hero(
        "🎯 Guide design bench",
        "Scan real sequence for real PAMs, rank the candidates, and then face the off-target "
        "table. The score is a transparent heuristic — the app shows every component that "
        "produced it.",
        ["IUPAC PAM scan", "off-target search", "donor & pegRNA design"],
    ), unsafe_allow_html=True)

    # ------------------------------------------------------------ target ----
    st.markdown("## 1 · Target")
    c1, c2 = st.columns([1.3, 1.7])
    with c1:
        locus_ids = [l["id"] for l in all_loci()]
        use_paste = st.toggle("Paste a custom sequence instead", value=False)
        if not use_paste:
            locus_id = st.selectbox("Genome library", locus_ids,
                                    format_func=lambda i: get_locus(i)["name"])
            locus = get_locus(locus_id)
            sequence = coding_strand(locus)
            st.caption(f"{locus['locus']} · {locus['strand']} strand gene · "
                       f"{locus['length']:,} bp · {locus['source']}")
        else:
            pasted = st.text_area("Sequence (ACGT, any length)", height=120,
                                  placeholder="paste a gene, an amplicon, a plasmid…")
            sequence = clean(pasted) if pasted else ""
            locus = None
    with c2:
        if sequence:
            st.markdown(theme.kv({
                "length": f"{len(sequence):,} bp",
                "GC": f"{100*gc_content(sequence):.1f}%",
                "N content": f"{100*sequence.count('N')/len(sequence):.1f}%",
                "source": locus["source"] if locus else "user input",
            }), unsafe_allow_html=True)
        st.caption("Everything downstream is computed on this exact sequence: PAM positions, "
                   "cut coordinates, off-target hits.")

    if not sequence:
        st.info("Choose a locus or paste a sequence to continue.")
        return

    # --------------------------------------------------------- effector -----
    st.markdown("## 2 · Effector")
    e1, e2, e3 = st.columns([1.4, 1, 1])
    with e1:
        effector_id = st.selectbox("Nuclease", [e["id"] for e in NUCLEASES], index=0,
                                   format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
        enzyme = EFFECTOR_BY_ID[effector_id]
    with e2:
        max_mm = st.slider("Off-target mismatch budget", 0, 5, 3,
                           help="How many mismatches a site may carry and still be counted "
                                "as a potential off-target. Fewer is stricter.")
    with e3:
        top_n = st.slider("Guides to score deeply", 5, 60, 20, step=5)

    region = st.slider("Target window (bp in this sequence)",
                       1, max(2, len(sequence)), (1, len(sequence)))
    region_0 = (region[0] - 1, region[1])

    run = st.button("🔍 Design and rank guides", type="primary")
    if run:
        with st.spinner("scanning PAMs and searching the off-target space…"):
            guides = design_guides(sequence, enzyme, region=region_0,
                                   max_offtarget_mismatches=max_mm, offtarget_top_n=top_n)
        st.session_state["guides"] = guides
        st.session_state["guides_key"] = (locus["id"] if locus else "custom", effector_id, region)

    guides = st.session_state.get("guides") or []
    if not guides:
        st.info("No guides designed yet — press the button. If nothing is found, the PAM is "
                "simply absent from this window: try a PAM-flexible variant (SpG, SpRY) or "
                "Cas12a for AT-rich sequence.")
        return

    st.markdown("## 3 · Ranked guides")
    st.plotly_chart(viz.fig_guide_scatter(guides), width="stretch")
    st.dataframe(_guide_frame(guides), width="stretch", hide_index=True, height=320)
    st.download_button("⬇️ Download the guide table (CSV)",
                       _guide_frame(guides).to_csv(index=False), "guides.csv")

    st.markdown("## 4 · Guide detail")
    g1, g2 = st.columns([1.5, 1])
    with g1:
        pick = st.selectbox("Guide", list(range(len(guides))),
                            format_func=lambda i: f"#{i+1}  {guides[i]['guide']}  "
                                                  f"({guides[i]['pam']}, {guides[i]['strand']})")
        g = guides[pick]
        site = g["site"]
        start = max(0, site["start"] - 25)
        end = min(len(sequence), site["end"] + 25)
        ctx = sequence[start:end]
        hl = [(site["start"] - start + 1, site["end"] - start, theme.C["violet"], "protospacer"),
              (site["pam_start"] - start + 1, site["pam_end"] - start, theme.C["gold"], "PAM")]
        st.markdown(theme.render_sequence(ctx, start=start + 1, width=110, highlights=hl),
                    unsafe_allow_html=True)
        st.markdown(theme.kv({
            "cut at": g["cut_site_1"],
            "strand": g["strand"],
            "geometry": f"{enzyme['cut_geometry']} · Δ{enzyme.get('overhang',0)} nt",
            "Tm of spacer": f"{g['tm']:.1f} °C",
            "homopolymer": g["homopolymer"],
        }), unsafe_allow_html=True)
        if g["warnings"]:
            for w in g["warnings"]:
                st.markdown(theme.callout("⚠️ " + w, "warn"), unsafe_allow_html=True)
    with g2:
        rna = build_guide_rna(g["guide"], enzyme)
        st.markdown(theme.card(f"Guide RNA · {rna['len']} nt",
                               f"<span class='mono' style='word-break:break-all'>{rna['rna']}</span>",
                               "violet"), unsafe_allow_html=True)
        comp = g["components"]
        st.plotly_chart(viz.fig_bar(list(comp.keys()), [100 * v for v in comp.values()],
                                    "Score components (%)", theme.C["teal"], 300,
                                    horizontal=True), width="stretch")
        st.markdown(theme.callout(
            "Score weights: GC 30%, position preferences 20%, homopolymer 15%, U6 poly-T 10%, "
            "self-complementarity 10%, specificity 15%. <b>This is a heuristic, not the "
            "Doench/Azimuth Rule Set 2 model</b> — use it to compare guides, not to predict a "
            "percentage.", "info"), unsafe_allow_html=True)

    # -------------------------------------------------------- off-targets ---
    st.markdown("## 5 · Off-target sites (≤ %d mismatches)" % max_mm)
    ots = find_offtargets(sequence, g["guide"], enzyme, max_mm)
    if not ots:
        st.markdown(theme.callout("No near-matches inside this sequence — but note that a "
                                  "1600 bp locus is not a genome. Real off-target searches run "
                                  "against the whole genome (and against the intended target "
                                  "for a therapeutic guide).", "good"), unsafe_allow_html=True)
    else:
        st.dataframe(pd.DataFrame([{
            "mismatches": o["mismatches"],
            "seed mismatches": o["seed_mismatches"],
            "strand": o["strand"],
            "position": f"{o['start']+1}-{o['end']}",
            "site": o["sequence"],
            "mismatch map": o["mismatch_string"],
            "PAM": o["pam"] or "—",
            "engages?": "yes" if o["pam_valid"] else "no PAM",
        } for o in ots]), width="stretch", hide_index=True)
        st.caption("`.` = match, `|` = mismatch. A hit with zero seed mismatches and a valid PAM "
                   "is the kind of site that shows up as a real off-target cut in a GUIDE-seq "
                   "experiment.")

    # ------------------------------------------------------------ donor -----
    st.markdown("## 6 · Donor design (if you need precise repair)")
    d1, d2, d3 = st.columns([1, 1, 1])
    with d1:
        edit = st.text_input("Desired edit (bases to write at the cut)", value="",
                             placeholder="e.g. GAG→GCG, or the corrected base")
        arm = st.slider("Homology arm (bp each side)", 20, 200, 60, step=10)
    with d2:
        donor_type = st.selectbox("Donor format", ["ssODN", "dsDNA", "plasmid", "aav",
                                                    "lentivirus"])
    with d3:
        offset = st.slider("Edit offset from the cut", -10, 10, 0)
    if edit:
        from ..crispr import design_hdr_donor
        donor = design_hdr_donor(sequence, g["cut_site_1"], clean(edit), arm=arm,
                                 edit_offset=offset)
        st.markdown(theme.render_sequence(
            donor["donor"], width=110,
            highlights=[(len(donor["left_arm"]) + 1,
                         len(donor["left_arm"]) + len(clean(edit)), theme.C["green"], "edit"),
                        (1, len(donor["left_arm"]), theme.C["blue"], "left arm"),
                        (len(donor["left_arm"]) + len(clean(edit)) + 1, len(donor["donor"]),
                         theme.C["teal"], "right arm")]),
            unsafe_allow_html=True)
        st.markdown(theme.kv({
            "donor": f"{donor['length']} nt", "format": donor_type,
            "left arm Tm": f"{donor['left_tm']:.1f} °C",
            "right arm Tm": f"{donor['right_tm']:.1f} °C",
            "GC": f"{100*donor['gc']:.1f}%",
        }), unsafe_allow_html=True)
        st.markdown(theme.callout(donor["silent_pam_mutation"], "info"), unsafe_allow_html=True)

    # ----------------------------------------------------------- pegRNA -----
    if PRIME_EDITORS and st.toggle("Design a pegRNA instead (prime editing)"):
        pe = st.selectbox("Prime editor", [p["id"] for p in PRIME_EDITORS], index=1,
                          format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
        p1, p2, p3 = st.columns(3)
        with p1:
            pbs_len = st.slider("PBS length", 6, 18, 13)
        with p2:
            rtt_len = st.slider("RT template length", 8, 40, 16)
        with p3:
            motif = st.toggle("3' protective motif", value=True)
        pe_edit = st.text_input("Sequence to install", value="", key="pe_edit",
                                placeholder="bases to write after the nick")
        from ..crispr import design_pegRNA
        try:
            peg = design_pegRNA(sequence, g["guide"], enzyme, clean(pe_edit or "A"),
                                edit_offset=len(g["guide"]) // 2, pbs_len=pbs_len,
                                rtt_len=rtt_len, add_motif=motif)
            st.markdown(theme.kv({
                "PBS": peg["pbs"], "PBS Tm": f"{peg['pbs_tm']:.1f} °C",
                "PBS GC": f"{100*peg['pbs_gc']:.0f}%",
                "RT template": peg["rtt"],
                "nick at": peg["nick_1"],
                "pegRNA length": len(peg["pegRNA"]),
            }), unsafe_allow_html=True)
            st.markdown(theme.card("pegRNA architecture",
                                   f"<span class='mono' style='word-break:break-all'>{peg['pegRNA']}</span>",
                                   "green"), unsafe_allow_html=True)
            st.markdown(theme.callout(peg["note"], "info"), unsafe_allow_html=True)
            tm_table = pd.DataFrame([{"PBS length": L, "Tm (°C)": round(tm, 1), "PBS": p}
                                     for L, tm, p in recommend_pbs(peg["rtt"], 45.0)])
            st.dataframe(tm_table, width="stretch", hide_index=True)
        except ValueError as exc:
            st.warning(f"Could not place the guide in this sequence: {exc}")

    # --------------------------------------------------------- PAM sweep ----
    st.markdown("## 7 · PAM coverage of this sequence")
    rows = []
    from ..data.cas_enzymes import PAM_COVERAGE_SETS
    for name, pam in PAM_COVERAGE_SETS.items():
        cov = pam_coverage(sequence, pam, window=50)
        cov["pam"] = name
        rows.append(cov)
    st.plotly_chart(viz.fig_pam_coverage(rows), width="stretch")
    st.caption("Same sequence, different PAM requirements. This is the quantitative version of "
               "«there is no NGG near my cut site».")
