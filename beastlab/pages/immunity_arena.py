"""Immunity arena — real arrays, real repeat/spacer architecture, real arms races."""
from __future__ import annotations

import zlib
from typing import Any, Dict, List, Tuple

import pandas as pd
import streamlit as st
from ..presentation import table, chart

from .. import theme, viz
from ..data.microbes import ANTI_CRISPR, BACTERIA, BACTERIA_BY_ID, PHAGES
from ..data.targets import load_records
from ..immunity import (candidate_protospacers, immunity_efficiency, interference_scan,
                        moist_curve, parse_crispr_array, rebuild_array, seed_matched_fraction,
                        simulate_adaptation, simulate_phage_challenge)
from ..sequtils import clean, gc_content, random_sequence

PHAGE_BY_ID = {p["id"]: p for p in PHAGES}
REAL_ARRAY_RECORD = "ecoli_k12_crispr1_locus"


def _stable_seed(text: str, mod: int = 99991) -> int:
    """Process-independent seed (Python's hash() is randomised per run)."""
    return zlib.crc32(text.encode()) % mod


def _invader_genome(phage: Dict[str, Any]) -> Tuple[str, str]:
    """Real genome when the FASTA ships with the app, otherwise a seeded synthetic one."""
    file = phage.get("file")
    if file:
        recs = load_records(file)
        for rec in recs.values():
            if "NC_" in rec["id"] or rec["id"] in ("NC_001422.1", "phix174"):
                return clean(rec["sequence"]), f"real genome · {rec['id']} · {rec.get('source','')}"
    seed = _stable_seed(str(phage["id"]))
    n = int(phage["genome_size"])
    return (random_sequence(n, gc=float(phage["gc"]), seed=seed),
            f"synthetic genome, {n:,} bp at {100*phage['gc']:.0f}% GC (seeded from the phage id)")


def _array_for(host: Dict[str, Any], system: Dict[str, Any],
               invader: str) -> Tuple[Dict[str, Any], str]:
    """Real array from the shipped FASTA when possible, otherwise a modelled history."""
    real = load_records("ecoli_k12.fasta").get(REAL_ARRAY_RECORD) if system.get("array_file") else None
    repeat = system.get("repeat") or "GTTTTAGAGCTGTGTTGTTTCGAATGGTTCCAAAAC"
    if real and host["id"] == "ecoli_k12":
        arr = parse_crispr_array(real["sequence"], repeat=repeat)
        return arr, (f"Real array, parsed live from `ecoli_k12.fasta` "
                     f"({real['id']}, {system['array_region']}).")
    hist = candidate_protospacers(invader, system, limit=8)
    spacers = [c["protospacer"] for c in hist[:6]]
    arr = {
        "repeat": repeat, "spacers": spacers, "repeats": len(spacers) + 1,
        "spacer_objects": [], "leader": "", "leader_length": 0,
        "array_start": 0, "array_end": 0,
        "spacer_lengths": [len(s) for s in spacers],
        "array_sequence": rebuild_array(repeat, spacers),
    }
    return arr, (f"Modelled exposure history: six spacers derived from this invader's own "
                 f"{system['pam']} sites, using the published {system['name']} repeat "
                 f"({system['repeat_source']}).")


def _genome_track(invader: str, hits: List[Dict[str, Any]],
                  new_spacers: List[str]) -> str:
    feats: List[Dict[str, Any]] = []
    for h in hits:
        feats.append({"start": h["start"], "end": h["end"], "label":
                      f"spacer {h['spacer_index']} targets here (PAM {h['pam'] or 'none'})",
                      "color": theme.C["gold"] if h["interference"] else theme.C["red"],
                      "top": 6, "height": 20})
    inv = clean(invader)
    for sp in new_spacers:
        pos = inv.find(sp)
        if pos >= 0:
            feats.append({"start": pos, "end": pos + len(sp), "label": "newly acquired spacer",
                          "color": theme.C["green"], "top": 29, "height": 9})
    return theme.genome_browser_html(invader, feats)


def render() -> None:
    st.markdown(theme.hero(
        "🛡️ Immunity arena",
        "A bacterium is not a bag of enzymes — it is a lineage with a memory of past infections, "
        "written into its own chromosome. Read a real array, give it a real invader, and watch "
        "what the defence does.",
        ["real array from U00096.3", "6 phages", "9 anti-CRISPRs", "arms-race ODE"],
    ), unsafe_allow_html=True)

    # -------------------------------------------------------------- set-up --
    c1, c2, c3 = st.columns([1.3, 1.3, 1])
    with c1:
        host_id = st.selectbox("Host bacterium", [b["id"] for b in BACTERIA],
                               format_func=lambda i: BACTERIA_BY_ID[i]["name"])
        host = BACTERIA_BY_ID[host_id]
    with c2:
        sys_idx = st.selectbox("CRISPR–Cas system",
                               list(range(len(host["systems"]))),
                               format_func=lambda i: host["systems"][i]["name"])
        system = host["systems"][sys_idx]
    with c3:
        phage_id = st.selectbox("Invading phage", [p["id"] for p in PHAGES],
                                format_func=lambda i: PHAGE_BY_ID[i]["name"])
        phage = PHAGE_BY_ID[phage_id]

    invader, invader_note = _invader_genome(phage)
    array, array_note = _array_for(host, system, invader)
    base_spacers = [clean(s) for s in array["spacers"]]

    # ------------------------------------------------ adaptation (up front) --
    with st.expander("🧬 Adaptation controls — let Cas1–Cas2 write a new spacer", expanded=False):
        a1, a2, a3 = st.columns(3)
        with a1:
            n_new = st.slider("New spacers to acquire", 1, 5, 2)
        with a2:
            mode = st.selectbox("Acquisition mode", ["naive", "primed"],
                                help="Primed adaptation is biased towards the neighbourhood of a "
                                     "protospacer that Cascade is already bound to.")
        with a3:
            seed = st.number_input("Seed", 0, 9999, _stable_seed(host_id + phage_id) % 9999,
                                   key="adapt_seed")
        primed_target = (base_spacers[0] if base_spacers else None)
        adapt = simulate_adaptation(invader, system, array, n_new=int(n_new), mode=mode,
                                    existing_target=primed_target, seed=int(seed))
        new_spacers = [clean(s) for s in adapt["new_spacers"]]
        include_new = st.toggle(
            "Count the newly acquired spacers as part of the defence", value=bool(new_spacers),
            help="Acquisition is a population event and it takes time to spread, but you can "
                 "watch what happens once the memory is written.")

    spacers = (new_spacers if include_new else []) + base_spacers
    hits = interference_scan(invader, spacers, system)
    active = [h for h in hits if h["interference"]]

    acr_present = bool(phage.get("anti_crispr"))
    acr_strength = st.slider("Anti-CRISPR strength", 0.0, 0.99,
                             min(0.95, 0.90) if acr_present else 0.0, step=0.01,
                             help="How completely the Acr protein disables the system. Only "
                                  "meaningful for phages that carry anti-CRISPR genes.")
    efficiency = immunity_efficiency(len(active), system, acr_present=acr_present,
                                    acr_strength=acr_strength)

    st.markdown(theme.kv({
        "host": f"{host['name']} · {host['genome_mb']} Mb · {100*host['gc']:.0f}% GC",
        "system": f"type {system['type']} ({system.get('subtype','—')}) · {system['effector']}",
        "PAM": f"{system['pam']} ({system.get('pam_iupac','—') or 'none'})",
        "repeat": f"{system['repeat'] or '—'} ({system['repeat_len']} nt, {system['repeat_source']})",
        "invader": f"{phage['name']} · {len(invader):,} bp · {100*gc_content(invader):.1f}% GC",
    }), unsafe_allow_html=True)
    st.caption(f"{array_note} {invader_note}")

    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.markdown(theme.metrics([("Spacers in the array", str(len(spacers)),
                                    f"{array['repeats']} repeats")]), unsafe_allow_html=True)
    with m2:
        st.markdown(theme.metrics([("Spacers matching this phage", str(len(active)),
                                    f"{len(hits)} exact matches found")]), unsafe_allow_html=True)
    with m3:
        st.markdown(theme.metrics([("Per-cell immunity", f"{100*efficiency:.1f}%",
                                    "chance an infection is aborted")]), unsafe_allow_html=True)
    with m4:
        st.markdown(theme.metrics([("Anti-CRISPR", ", ".join(phage.get("anti_crispr", [])) or "none",
                                    "carried by this phage")]), unsafe_allow_html=True)

    tabs = st.tabs(["① The array", "② Adaptation", "③ Interference",
                    "④ Arms race", "⑤ Escape playground", "⑥ Anti-CRISPR"])

    # -------------------------------------------------------------- array ---
    with tabs[0]:
        st.markdown("### Repeat–spacer architecture")
        st.markdown(theme.crispr_array_html(array["repeat"], spacers,
                                            new_count=0,
                                            leader=array.get("leader", "") or ""),
                    unsafe_allow_html=True)
        if st.toggle("Show the array as raw sequence", value=False):
            st.markdown(theme.card("Locus",
                                   f"<span class='mono' style='word-break:break-all'>"
                                   f"{array['array_sequence']}</span>", "violet"),
                        unsafe_allow_html=True)
        rows = []
        for i, sp in enumerate(spacers):
            matching = [h for h in hits if h["spacer_index"] == i]
            rows.append({
                "spacer": f"S{len(spacers)-i}", "length": len(sp),
                "GC": f"{100*gc_content(sp):.0f}%",
                "sequence (5'→3')": sp,
                "targets this invader": "yes" if any(h["interference"] for h in matching)
                                        else ("no PAM" if matching else "no"),
            })
        table(pd.DataFrame(rows), width="stretch", hide_index=True, height=280)
        chart(viz.fig_array_map(array), width="stretch")
        st.markdown(theme.callout(
            "Evolutionary reading: the newest spacers sit next to the leader (5' end) because "
            "Cas1–Cas2 always inserts there. A real array is a chronological record — the "
            "leader-distal spacers are the oldest immunity.", "info"), unsafe_allow_html=True)
        st.markdown(theme.callout(system["activity"], "info"), unsafe_allow_html=True)

    # --------------------------------------------------------- adaptation ---
    with tabs[1]:
        st.markdown("### Adaptation — writing a new memory")
        st.caption("Controls for this bench sit in the **Adaptation controls** panel above, so that "
                   "the acquired spacers flow into the interference and arms-race numbers.")
        st.markdown(theme.metrics([
            ("Candidates in this genome", str(adapt.get("candidates", 0)),
             f"{system['pam']} sites long enough for a {system['protospacer_len']} nt spacer"),
            ("New spacers", str(len(new_spacers)), f"mode: {adapt.get('mode', mode)}"),
            ("Counted in the defence", "yes" if include_new else "no",
             "see the toggle in the controls panel"),
        ]), unsafe_allow_html=True)
        if new_spacers:
            ordered = new_spacers + base_spacers
            st.markdown(theme.crispr_array_html(array["repeat"], ordered,
                                                new_count=len(new_spacers),
                                                leader=array.get("leader", "") or ""),
                        unsafe_allow_html=True)
            table(pd.DataFrame([{
                "event": e["event"], "new spacer": e["spacer"], "PAM": e["pam"],
                "position in the invader": e["invader_position"], "GC": f"{100*e['gc']:.0f}%",
            } for e in adapt["events"]]), width="stretch", hide_index=True)
            st.markdown(theme.callout(adapt.get("note", ""), "good"), unsafe_allow_html=True)
            if st.toggle("Where do the new spacers point in the invader genome?"):
                st.markdown(_genome_track(invader, hits, new_spacers), unsafe_allow_html=True)
                st.caption("Gold = a spacer in the array that hits this genome (red if the PAM "
                           "is missing); green = the spacers just acquired.")
            st.session_state["arena_new_spacers"] = new_spacers
        else:
            st.warning("No candidate protospacers were found in this invader for this system. "
                       "That is a real result: an AT-rich genome with a GC-rich PAM has few "
                       "targetable sites.")

    # -------------------------------------------------------- interference --
    with tabs[2]:
        st.markdown("### Interference — does the array fire?")
        if hits:
            table(pd.DataFrame([{
                "spacer": f"S{len(spacers)-h['spacer_index']}",
                "strand": h["strand"],
                "position": f"{h['start']+1}–{h['end']}",
                "PAM found": h["pam"] or "—",
                "PAM site": h["pam_side"],
                "engages": "YES" if h["interference"] else "no",
                "note": h["note"],
            } for h in hits]), width="stretch", hide_index=True)
        else:
            st.markdown(theme.callout(
                "None of the spacers in this array matches this invader — the population has "
                "never met this phage. <b>Three ways forward:</b> open the adaptation controls "
                "and let Cas1–Cas2 write a spacer against this very genome (then re-read this "
                "table), pick a different invader, or accept the result — E. coli K-12's "
                "CRISPR-1 spacers were acquired against other invaders, and the locus is "
                "transcriptionally silenced at 37 °C anyway.", "warn"), unsafe_allow_html=True)
        st.markdown("#### Match quality for the best spacers")
        probe = active[0]["protospacer"] if active else (
            max(spacers, key=len)[:system["protospacer_len"]] if spacers else "")
        st.markdown(theme.kv({
            "exact matches": str(len(hits)),
            "engaging (PAM-checked)": str(len(active)),
            "per-cell immunity": f"{100*immunity_efficiency(len(active), system):.1f}%",
            "with anti-CRISPR": f"{100*immunity_efficiency(len(active), system, True, acr_strength):.1f}%"
                                if acr_present else "n/a — this phage carries none",
        }), unsafe_allow_html=True)
        st.caption("This model requires an exact protospacer match plus a valid PAM (a "
                   "mismatch-tolerant version with seed scoring would engage more sites — the "
                   "approximation is stated rather than hidden).")

        st.markdown("#### Poisson infection at a given MOI")
        moi = st.slider("MOI (phage per bacterium)", 0.01, 50.0, 5.0, step=0.01)
        curve = moist_curve(moi, 1.0, efficiency)
        chart(viz.fig_moi(curve), width="stretch")
        st.markdown(theme.kv({
            "P(sensitive cell infected)": f"{100*curve['p_infected_sensitive']:.1f}%",
            "P(immune cell survives)": f"{100*curve['p_survive_if_immune']:.1f}%",
            "protection factor": f"{curve['protection_factor']:.1f}×",
        }), unsafe_allow_html=True)
        st.caption("At high MOI even a 97%-effective single spacer leaks: immunity is "
                   "probabilistic, which is why real defences stack several spacers.")

    # ------------------------------------------------------------ arms race -
    with tabs[3]:
        st.markdown("### The arms race")
        st.caption("A deterministic six-compartment model: sensitive bacteria, immune bacteria, "
                   "CRISPR-suppressed (anti-CRISPR) bacteria, wild-type phage, escaped phage and "
                   "anti-CRISPR phage.")
        r1, r2, r3, r4 = st.columns(4)
        with r1:
            hours = st.slider("Duration (h)", 4, 72, 24)
        with r2:
            moi0 = st.slider("Starting MOI (phage per bacterium)", 0.01, 20.0, 1.0, step=0.01)
        with r3:
            immune_frac = st.slider("Fraction of the host population that is immune", 0.0, 1.0,
                                    0.2 if len(active) else 0.0, step=0.01,
                                    help="Even a perfect spacer does not help a population in "
                                         "which nobody carries it yet.")
        with r4:
            escape = st.number_input("Phage escape rate", 1e-7, 1e-2, 2e-6, format="%.1e")

        r5, r6, r7, r8 = st.columns(4)
        with r5:
            acquire = st.number_input("Spacer acquisition rate", 1e-6, 1e-1, 1e-4, format="%.1e")
        with r6:
            immunity_cost = st.slider("Cost of immunity", 0.0, 0.2, 0.03, step=0.01)
        with r7:
            initial_b = st.number_input("Initial bacteria (CFU/mL)", 1e4, 1e9, 1e7, format="%.0e")
        with r8:
            acr_fraction = st.slider("Phage carrying anti-CRISPR", 0.0, 1.0,
                                     0.3 if acr_present else 0.0, step=0.05)

        initial_p = float(moi0) * float(initial_b)
        st.caption(f"Starting phage density implied by that MOI: {initial_p:.2e} PFU/mL.")
        spa_frac = min(immune_frac, immune_frac * acr_fraction)
        run = st.button("▶ Run the infection", type="primary")
        if run or "arena_challenge" in st.session_state:
            if run:
                st.session_state["arena_challenge"] = simulate_phage_challenge(
                    immune_fraction=float(immune_frac - spa_frac),
                    spa_fraction=float(spa_frac),
                    immunity_efficiency_per_cell=float(efficiency if efficiency else 0.0),
                    acr_fraction=float(acr_fraction),
                    acquisition_rate=float(acquire),
                    escape_rate=float(escape),
                    initial_bacteria=float(initial_b),
                    initial_phage=initial_p,
                    hours=float(hours),
                    immunity_cost=float(immunity_cost),
                    burst_size=float(phage["burst_size"]),
                    latency_h=float(phage["latency_min"]) / 60.0,
                )
            res = st.session_state["arena_challenge"]
            chart(viz.fig_arms_race(res["series"], res["events"]),
                            width="stretch")
            f = res["final"]
            st.markdown(theme.metrics([
                ("Bacteria", f"{f['bacteria']:.2e}", "CFU/mL at the end"),
                ("Phage", f"{f['phage']:.2e}", "PFU/mL at the end"),
                ("Immune fraction", f"{f['immune_percent']:.1f}%", "of surviving bacteria"),
                ("Escaped phage", f"{f['escape_phage_percent']:.2f}%", "of the phage population"),
            ]), unsafe_allow_html=True)
            st.markdown(theme.card("Outcome", res["outcome"],
                                   "green" if "bacteria win" in res["outcome"] else "plasma"),
                        unsafe_allow_html=True)
            if res["events"]:
                st.markdown("#### Event log")
                for ev in res["events"]:
                    st.markdown(f"- **{ev['t']:.2f} h** — {ev['event']}: {ev['detail']}")
            with st.expander("Assumptions and parameters"):
                st.markdown(theme.kv(res["parameters"]), unsafe_allow_html=True)
                for a in res["assumptions"]:
                    st.markdown(f"- {a}")
        else:
            st.info("Press **Run the infection**. The default MOI is 1 phage per bacterium at "
                    "t = 0.")

    # ----------------------------------------------------- escape playground -
    with tabs[4]:
        st.markdown("### Escape playground")
        st.caption("A single base change can defeat a spacer, but only if it lands in the right "
                   "place. Mutate the target and watch the immunity status flip.")
        if not spacers:
            st.info("This array has no spacers to work with.")
        else:
            if active:
                base_hit = active[0]
                ctx_seq = clean(invader)[max(0, base_hit["start"] - 4):
                                         base_hit["end"] + 4]
            else:
                base_hit = None
                ctx_seq = clean(invader)[:int(system.get("protospacer_len", 32)) + 8]
            spacer = base_hit["spacer"] if base_hit else spacers[0]
            e1, e2, e3 = st.columns([1.4, 1, 1])
            with e1:
                st.markdown(theme.render_sequence(ctx_seq, start=1, width=80),
                            unsafe_allow_html=True)
            with e2:
                pos = st.slider("Position to mutate (1-based)", 1, len(ctx_seq), 5)
                bases = ["A", "C", "G", "T"]
                new_base = st.selectbox("New base", bases,
                                        index=bases.index(ctx_seq[pos - 1]) if pos <= len(ctx_seq)
                                        else 0)
            with e3:
                if st.button("🧬 Mutate and re-test"):
                    mutated = ctx_seq[:pos - 1] + new_base + ctx_seq[pos:]
                    mhits = interference_scan(mutated, [spacer], system)
                    engaging = any(h["interference"] for h in mhits)
                    st.session_state["escape_result"] = (ctx_seq, mutated, engaging, mhits)
            if "escape_result" in st.session_state:
                old, mutated, engaging, mhits = st.session_state["escape_result"]
                st.markdown(viz.sequence_diff_html(old, mutated, context=10),
                            unsafe_allow_html=True)
                if engaging:
                    st.markdown(theme.callout(
                        "Immunity survives this mutation: an escape variant needs the change to "
                        "hit the protospacer (and, for most types, to leave the PAM intact).",
                        "good"), unsafe_allow_html=True)
                else:
                    st.markdown(theme.callout(
                        "Immunity lost. The array still carries the spacer, but this phage "
                        "variant is invisible to it — this is what the model calls an escape "
                        "mutant, and it is why escape phage appear in the arms race.", "warn"),
                        unsafe_allow_html=True)
            if base_hit:
                st.markdown(theme.kv({
                    "spacer": spacer,
                    "protospacer": base_hit["protospacer"],
                    "PAM": base_hit["pam"] or "none",
                    "seed match (last 8 nt)": f"{100*seed_matched_fraction(spacer, base_hit['protospacer']):.0f}%",
                }), unsafe_allow_html=True)
            st.markdown(theme.callout(
                "Where mutations hurt: the PAM and the PAM-proximal seed are the most "
                "escape-sensitive positions. A mismatch in the PAM-proximal seed usually "
                "kills targeting outright; mismatches at the PAM-distal end are often "
                "tolerated — which is exactly why cross-reactivity between related phages "
                "is common.", "info"), unsafe_allow_html=True)

    # --------------------------------------------------------- anti-CRISPR --
    with tabs[5]:
        st.markdown("### Anti-CRISPR proteins — the counter-defence")
        table(pd.DataFrame(ANTI_CRISPR), width="stretch", hide_index=True)
        st.markdown(theme.callout(
            "Anti-CRISPRs are why 'fully immune' is never fully immune. They are also a "
            "laboratory tool: an Acr protein can switch a CRISPR system off on demand, and "
            "Acr-based 'off switches' are being engineered for safer editing.", "info"),
            unsafe_allow_html=True)
        st.markdown("#### Which phages in this arena carry one")
        for p in PHAGES:
            carried = p.get("anti_crispr")
            st.markdown(f"- **{p['name']}** — "
                        + (f"carries {', '.join(carried)}" if carried else
                           "no known anti-CRISPR genes"), unsafe_allow_html=True)
        st.caption("Compiled from Bondy-Denomy 2013 Nature 493:429; Rauch 2017 Cell 168:150; "
                   "Watters 2018 Science 362:236; Meeske 2020 Nature 578:381.")
