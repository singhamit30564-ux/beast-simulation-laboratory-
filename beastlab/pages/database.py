"""Database — the reference layer: genomes, bacteria, phages, effectors, trials."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
import streamlit as st

from .. import theme, viz
from ..data.cas_enzymes import (ALL_DNA_EFFECTORS, BASE_EDITORS, EFFECTOR_BY_ID,
                                NUCLEASES, PAM_COVERAGE_SETS, POINT_MUTATION_CLASSES,
                                PRIME_EDITORS, RNA_EFFECTORS)
from ..data.clinical import DISEASES, TRIALS
from ..data.microbes import ANTI_CRISPR, BACTERIA, BACTERIA_BY_ID, PHAGES
from ..data.targets import GENOME_DIR, load_records, summary_table
from ..immunity import parse_crispr_array
from ..sequtils import (clean, complexity, find_motifs, gc_content, orfs, translate)

ROOT = Path(__file__).resolve().parents[2]
REFERENCES = ROOT / "data" / "references.md"

FAMILY_COLOR = {"Cas9": theme.C["gold"], "Cas12": theme.C["teal"], "Cas13": theme.C["violet"],
                "Cas12f": theme.C["teal"], "Cascade": theme.C["blue"], "Csm": theme.C["plasma"],
                "Cmr": theme.C["plasma"], "CBE": theme.C["green"], "ABE": theme.C["green"],
                "PE": theme.C["violet"]}


def _family_of(e: Dict[str, Any]) -> str:
    fam = e.get("family", "")
    for key in ("Cas12f", "Cas12j", "Cas12", "Cas13", "Cas9", "Cascade", "Csm", "Cmr", "CBE",
                "ABE", "PE"):
        if fam.startswith(key):
            return key
    return fam or "—"


def _records_for(filename: str) -> Dict[str, Dict[str, Any]]:
    return load_records(filename)


def _genome_tab() -> None:
    st.markdown("#### Shipped genomes")
    files = sorted(p.name for p in GENOME_DIR.glob("*.fasta")) if GENOME_DIR.exists() else []
    if not files:
        st.error(f"No FASTA files found in {GENOME_DIR}.")
        return
    c1, c2, c3 = st.columns([1.2, 1.2, 1])
    with c1:
        filename = st.selectbox("FASTA file", files)
    records = _records_for(filename)
    if not records:
        st.warning("That file could not be parsed into records.")
        return
    with c2:
        rec_id = st.selectbox("Record", list(records))
    rec = records[rec_id]
    seq = clean(rec["sequence"])
    with c3:
        motif = st.text_input("Search motif (IUPAC ok)", value="ATG")

    st.markdown(theme.metrics([
        ("Length", f"{len(seq):,} bp", filename),
        ("GC", f"{100*gc_content(seq):.1f}%", "of the record"),
        ("Sequence complexity", f"{complexity(seq):.2f}", "distinct 4-mers / total 4-mers"),
        ("ORFs ≥ 30 aa", f"{len(orfs(seq, 30))}", "six frames, no splicing"),
    ]), unsafe_allow_html=True)
    st.caption(rec.get("header", rec_id))

    hits = find_motifs(seq, motif, both_strands=True) if motif else []
    features: List[Dict[str, Any]] = [
        {"start": s, "end": e, "label": f"ORF {e-s} bp", "color": theme.C["violet"],
         "top": 8, "height": 22} for (s, e, _aa) in orfs(seq, 50)[:60]]
    features += [{"start": pos, "end": pos + len(motif), "label": f"{motif} ×{i+1}",
                  "color": theme.C["gold"], "top": 30, "height": 8}
                 for i, (pos, _m) in enumerate(hits[:80])]
    st.markdown(theme.genome_browser_html(seq, features), unsafe_allow_html=True)
    st.caption(f"Track: {len(features) - len(hits[:80])} ORFs (violet above) and "
               f"{len(hits)} {motif} motifs (gold below, both strands).")

    if hits:
        st.markdown("##### Motif positions")
        half = 24
        rows = []
        for pos, matched in hits[:60]:
            rows.append({"position (1-based)": pos + 1, "matched": matched,
                         "context": seq[max(0, pos - 6):pos + len(matched) + 6]})
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True, height=220)
    else:
        st.info(f"No {motif} motif in this record on either strand.")

    # sequence viewer
    st.markdown("##### Sequence")
    start = st.slider("Start at (bp)", 1, max(1, len(seq) - 60), 1, step=1)
    chunk = seq[start - 1:start - 1 + 600]
    hl = [(h[0] + 1 - (start - 1), h[0] + len(motif) - (start - 1), theme.C["gold"], motif)
          for h in hits if start - 1 <= h[0] < start - 1 + 600]
    st.markdown(theme.render_sequence(chunk, start=start, width=100, highlights=hl[:40]),
                unsafe_allow_html=True)

    with st.expander("Translation (frame +1) and ORF table"):
        protein = translate(seq, 0)
        st.markdown(theme.card("Frame +1 translation",
                               f"<span class='mono' style='word-break:break-all'>"
                               f"{protein[:600]}{'…' if len(protein) > 600 else ''}</span>",
                               "violet"), unsafe_allow_html=True)
        orf_rows = [{"start": s + 1, "end": e, "length (aa)": (e - s) // 3,
                     "first 40 aa": aa[:40]} for (s, e, aa) in orfs(seq, 60)]
        if orf_rows:
            st.dataframe(pd.DataFrame(orf_rows).sort_values("length (aa)", ascending=False),
                         width="stretch", hide_index=True, height=240)

    st.download_button("⬇️ Download this record as FASTA",
                       f">{rec_id}\n{seq}\n", f"{rec_id}.fasta")


def _bacteria_tab() -> None:
    c1, c2 = st.columns([1, 2])
    with c1:
        bid = st.selectbox("Bacterium", [b["id"] for b in BACTERIA],
                           format_func=lambda i: BACTERIA_BY_ID[i]["name"])
        b = BACTERIA_BY_ID[bid]
        st.markdown(theme.kv({
            "domain": b["domain"], "phylum": b["phylum"], "Gram": b["gram"],
            "morphology": b["morphology"], "genome": f"{b['genome_mb']} Mb",
            "GC": f"{100*b['gc']:.0f}%", "genes": f"{b['genes']:,}",
            "accession": b["accession"], "habitat": b["habitat"],
        }), unsafe_allow_html=True)
        st.markdown(theme.callout(b["interesting"], "info"), unsafe_allow_html=True)
    with c2:
        st.markdown(f"##### CRISPR–Cas systems in {b['short']}")
        for s in b["systems"]:
            kind = "gold" if s.get("repeat_source") == "sequenced" else "teal"
            st.markdown(theme.card(
                f"{s['name']} · type {s['type']}{('-' + s['subtype'].split('-')[-1]) if '-' in s.get('subtype','') else ''}",
                f"<b>Effector:</b> {s['effector']}<br>"
                f"<b>cas genes:</b> {', '.join(s['cas_genes'])}<br>"
                f"<b>PAM:</b> {s['pam']}<br>"
                f"<b>Protospacer:</b> {s['protospacer_len']} nt<br>"
                f"<b>Repeat:</b> <span class='mono'>{s['repeat'] or '—'}</span> "
                f"({s['repeat_len']} nt, <i>{s['repeat_source']}</i>)<br>"
                f"<b>Activity:</b> {s['activity']}<br>"
                f"<span class='muted'>{s['notes']}</span><br>"
                f"<span class='muted'>Ref: {s['ref']}</span>", kind), unsafe_allow_html=True)
            if s.get("repeat_source") == "sequenced" and (ROOT / "data" / "genomes" /
                                                          "ecoli_k12.fasta").exists():
                recs = _records_for("ecoli_k12.fasta")
                arr_rec = next((r for k, r in recs.items() if "crispr1" in k), None)
                if arr_rec:
                    parsed = parse_crispr_array(arr_rec["sequence"], repeat=s["repeat"])
                    st.markdown(theme.crispr_array_html(
                        parsed["repeat"], parsed["spacers"],
                        leader=parsed.get("leader", "") or ""), unsafe_allow_html=True)
                    st.caption(f"Parsed live from the shipped FASTA: {parsed['repeats']} repeats, "
                               f"{len(parsed['spacers'])} spacers, repeat {len(parsed['repeat'])} nt, "
                               f"leader {parsed.get('leader_length', 0)} nt.")
        st.markdown(theme.callout(
            "Repeats tagged <i>sequenced</i> were read out of a FASTA in this repository; "
            "<i>literature</i> repeats come from the cited paper and <i>type-level</i> entries "
            "are the canonical repeat of the type, not this strain's own array — treat those as "
            "illustrative.", "warn"), unsafe_allow_html=True)


def _phages_tab() -> None:
    st.markdown("#### Phages used in the immunity simulations")
    st.dataframe(pd.DataFrame([{
        "phage": p["name"], "host": p["host"], "family": p["family"],
        "genome": p["genome_type"], "size (bp)": p["genome_size"],
        "GC": f"{100*p['gc']:.0f}%", "genes": p["genes"], "accession": p["accession"],
        "latency (min)": p["latency_min"], "burst": p["burst_size"],
        "anti-CRISPR": ", ".join(p.get("anti_crispr", [])) or "—",
    } for p in PHAGES]), width="stretch", hide_index=True)
    st.markdown("#### Anti-CRISPR proteins")
    st.dataframe(pd.DataFrame(ANTI_CRISPR), width="stretch", hide_index=True)
    st.markdown(theme.callout(
        "Anti-CRISPRs are the reason 'fully immune' is never fully immune: a phage that carries "
        "an Acr gene can shut the defence down before it fires. In the immunity arena this is the "
        "anti-CRISPR switch on the invader.", "info"), unsafe_allow_html=True)
    st.markdown("#### Phage detail")
    pid = st.selectbox("Phage", [p["id"] for p in PHAGES],
                       format_func=lambda i: next(p["name"] for p in PHAGES if p["id"] == i))
    p = next(x for x in PHAGES if x["id"] == pid)
    st.markdown(theme.kv({
        "host": p["host"], "family": p["family"], "genome": p["genome_type"],
        "size": f"{p['genome_size']:,} bp", "GC": f"{100*p['gc']:.1f}%",
        "genes": p["genes"], "accession": p["accession"], "reference": p["ref"],
    }), unsafe_allow_html=True)
    st.markdown(p["notes"])
    st.caption("Phage parameters used by the simulation: virulence "
               f"{p['virulence']}, latency {p['latency_min']} min, burst {p['burst_size']}.")


def _cas_tab() -> None:
    st.markdown(f"#### Effector catalogue — {len(EFFECTOR_BY_ID)} entries")
    kind = st.radio("Show", ["All", "DNA nucleases", "RNA targeting", "Base editors",
                             "Prime editors"], horizontal=True)
    picks = {"All": NUCLEASES + RNA_EFFECTORS + BASE_EDITORS + PRIME_EDITORS,
             "DNA nucleases": ALL_DNA_EFFECTORS, "RNA targeting": RNA_EFFECTORS,
             "Base editors": BASE_EDITORS, "Prime editors": PRIME_EDITORS}[kind]

    rows = []
    for e in picks:
        rows.append({
            "id": e["id"], "name": e["name"], "family": _family_of(e),
            "Cas effector": e.get("cas", e.get("effector", "")),
            "PAM": e.get("pam", "—"),
            "pam side": e.get("pam_side", "—"),
            "spacer": e.get("spacer_len", e.get("protospacer_len", "—")),
            "cut": e.get("cut_geometry", "—"),
            "overhang": e.get("overhang", "—"),
            "PAM-proximal window": str(e.get("window", "—")),
            "ref": e.get("ref", ""),
        })
    df = pd.DataFrame(rows)
    st.dataframe(df, width="stretch", hide_index=True, height=380)
    st.download_button("⬇️ Download the catalogue (CSV)", df.to_csv(index=False),
                       "cas_catalogue.csv")

    st.markdown("##### Effector detail")
    eid = st.selectbox("Effector", [e["id"] for e in picks],
                       format_func=lambda i: EFFECTOR_BY_ID[i]["name"])
    e = EFFECTOR_BY_ID[eid]
    c1, c2 = st.columns([1.4, 1])
    with c1:
        st.markdown(viz.svg_cas_domains(e), unsafe_allow_html=True)
        st.markdown(theme.kv({k: str(v) for k, v in e.items() if k not in ("notes", "ref")}),
                    unsafe_allow_html=True)
    with c2:
        st.markdown(theme.card(e["name"], e.get("notes", ""), "gold"), unsafe_allow_html=True)
        st.caption(f"Reference: {e.get('ref','—')}")
        if "conversion" in e:
            st.markdown(theme.metrics([
                ("Conversion", e["conversion"], "chemistry"),
                ("Window", str(e.get("window", "—")), "protospacer positions"),
                ("Max efficiency", f"{100*e.get('max_efficiency', 0):.0f}%", "as modelled"),
                ("Indel rate", f"{100*e.get('indel_rate', 0):.1f}%", "DSB byproducts"),
            ]), unsafe_allow_html=True)
        if e["id"] in PAM_COVERAGE_SETS:
            st.caption("PAM density for this enzyme is charted in the PAM coverage bench.")

    st.markdown("#### Point mutations and the editors that can revert them")
    st.dataframe(pd.DataFrame(POINT_MUTATION_CLASSES), width="stretch", hide_index=True)


def _registry_tab() -> None:
    st.markdown("#### Loci shipped with the app")
    st.dataframe(pd.DataFrame(summary_table()), width="stretch", hide_index=True)
    st.markdown("#### Clinical programmes")
    st.dataframe(pd.DataFrame([{
        "product": t["product"], "modality": t["modality"], "indication": t["indication"],
        "phase": t["phase"], "target": t["target"], "ref": t["ref"],
    } for t in TRIALS]), width="stretch", hide_index=True)
    st.markdown("#### Disease models")
    st.dataframe(pd.DataFrame([{
        "key": k, "disease": v["name"], "programme": v["market"],
        "endpoint": v["endpoint"], "reference": v["ref"],
    } for k, v in DISEASES.items()]), width="stretch", hide_index=True)

    st.markdown("#### Provenance and limitations")
    if REFERENCES.exists():
        text = REFERENCES.read_text()
        st.caption(f"{REFERENCES.name} · {len(text.splitlines())} lines · "
                   f"{len(text):,} characters")
        with st.expander("Open the provenance file", expanded=False):
            st.markdown(text)
        st.download_button("⬇️ Download references.md", text, "references.md")
    else:
        st.warning("references.md was not found next to the shipped genomes.")

    st.markdown(theme.callout(
        "Every sequence in this app is either a real record (with its accession) or clearly "
        "labelled as a type-level illustration. Every number produced by the simulations is a "
        "model output with the assumptions printed next to it — none of it is a prediction you "
        "should act on clinically.", "warn"), unsafe_allow_html=True)


def render() -> None:
    st.markdown(theme.hero(
        "🗄️ Database",
        "The reference layer: the genomes that ship with the app, the bacteria and their CRISPR "
        "systems, the phages and their counter-defences, the full effector catalogue, and the "
        "provenance for every number.",
        [f"{len(NUCLEASES)} nucleases", f"{len(BACTERIA)} bacteria", f"{len(PHAGES)} phages",
         f"{len(TRIALS)} programmes"],
    ), unsafe_allow_html=True)

    tabs = st.tabs(["Genomes", "Bacteria & CRISPR", "Phages & anti-CRISPR", "Cas catalogue",
                    "Loci, trials & provenance"])
    with tabs[0]:
        _genome_tab()
    with tabs[1]:
        _bacteria_tab()
    with tabs[2]:
        _phages_tab()
    with tabs[3]:
        _cas_tab()
    with tabs[4]:
        _registry_tab()
