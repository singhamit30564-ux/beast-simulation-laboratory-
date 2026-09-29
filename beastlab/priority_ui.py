"""Button-gated priority pack UI. No torch imports on ordinary page loads."""
import streamlit as st
from .io import parse_dna, export_center


@st.cache_resource(show_spinner=False)
def efficiency_resource():
    from .ml.efficiency import load
    return load()


@st.cache_resource(show_spinner=False)
def dqn_resource():
    from .rl.dqn import load
    return load()


def guide_tools(sequence, source=None):
    with st.expander("Custom gRNA · GRCh38 scan · ML & RL · Export center", expanded=True):
        st.caption("SpCas9 sandbox: 20-nt spacer, adjacent NGG. Independent of the effector chosen in the ranking bench below.")
        guide_text = st.text_input("Custom gRNA / RNA spacer (20 nt)", value="ACGTACGTACGTACGTACGT", key="priority_guide")
        try:
            guide = parse_dna(guide_text, guide=True)
        except ValueError as exc:
            st.error(str(exc)); return
        from .offtarget import REGIONS, reference_or_fallback, scan
        name = st.selectbox("GRCh38 reference window (NOT whole genome)", list(REGIONS))
        st.caption("Ensembl REST public reference only; your guide is never sent. ≤3 substitutions, both strands, adjacent NGG, 12-nt PAM-proximal seed. No bulges/variants/chromatin. Exact matches are included; intended site is not auto-assigned. Risk is a heuristic, not measured cleavage.")
        if st.button("Scan vs GRCh38", key="scan_grch38"):
            with st.spinner("Fetching bounded public reference; offline fallback available…"):
                seq, meta, error = reference_or_fallback(name)
                hits = scan(seq, guide, offset=meta["start"], chromosome=meta["chromosome"])
            st.session_state["priority_scan"] = ((guide, name), {"hits": hits, "reference": meta, "fallback_reason": error})
        if st.button("Scan bundled phiX174 offline", key="scan_offline"):
            from .data.targets import GENOME_DIR
            from .sequtils import parse_fasta
            seq = parse_fasta((GENOME_DIR / "phix174.fasta").read_text())[0][1]
            st.session_state["priority_scan"] = ((guide, name), {"hits": scan(seq, guide, chromosome="phiX174"),
                "reference": {"source": "Bundled NCBI phiX174; NOT GRCh38", "start": 1, "end": len(seq)}, "fallback_reason": "Offline selected"})
        results = {}
        saved = st.session_state.get("priority_scan")
        if saved and saved[0] == (guide, name):
            results["offtarget"] = saved[1]
            if saved[1]["fallback_reason"]:
                st.warning("Offline phiX174 demonstration — NOT a human safety screen. " + saved[1]["fallback_reason"])
            st.write(saved[1]["reference"])
            st.write(f"{len(saved[1]['hits'])} matches within this bounded reference only.")
            from .presentation import table
            table(saved[1]["hits"])
            st.caption("What this does not claim: whole-GRCh38 coverage or absence of off-targets elsewhere. phiX174 is scanned linearly; origin-spanning circular sites are omitted.")
        st.markdown("#### Educational MLP · 24 → 32 → 16 → 1")
        st.warning("Simulated training data — illustrative model, not clinical-grade")
        pam = st.selectbox("Model PAM category (illustrative)", ["NGG", "NAG", "OTHER"])
        if st.button("Predict Efficiency"):
            from .ml.efficiency import predict
            model, meta = st.session_state.get("retrained_efficiency") or efficiency_resource()
            st.session_state["priority_prediction"] = ((guide, pam), predict(model, meta, guide, pam))
        saved = st.session_state.get("priority_prediction")
        if saved and saved[0] == (guide, pam):
            result = saved[1]
            results["efficiency"] = result
            st.metric("Illustrative cut efficiency", f"{result['efficiency_percent']:.1f}%")
            st.write("Illustrative confidence band (not clinical): %.1f–%.1f%%" % tuple(result["illustrative_band"]))
            st.caption(result["band_method"])
        with st.expander("Training page · simulated data / loss curve / retrain"):
            st.caption("Seed 42; 5,000 simulated guides (4,000 train / 1,000 held-out). 24 features: GC fraction + 16 overlapping dinucleotide fractions + 4 position-weighted quarters + 3 PAM categories. No experimental data. CPU-only; optional retraining affects this session's predictor, not the shipped RL reward model.")
            if st.button("Show training loss"):
                from .ml.efficiency import ROOT
                import json
                st.session_state["training_meta"] = json.loads((ROOT / "metadata.json").read_text())
            if st.button("Retrain MLP · seed=42"):
                from .ml.efficiency import train
                with st.spinner("Training on 5,000 synthetic guides in RAM…"):
                    model, meta = train()
                st.session_state["retrained_efficiency"] = (model, meta)
                st.session_state["training_meta"] = meta
                st.session_state.pop("priority_prediction", None)
            if "training_meta" in st.session_state:
                import pandas as pd
                loss = pd.DataFrame(st.session_state["training_meta"]["loss"])
                if st.session_state.get("lite_mode"):
                    st.dataframe(loss.tail(50), hide_index=True)
                else:
                    st.line_chart(loss.set_index("epoch"))
        st.markdown("#### 🧠 RL guide optimizer · sequence-score sandbox")
        st.caption("Assumptions: shipped MLP reward with NGG fixed; 80 mutations, max 10 steps or score >90%. Mutations do NOT preserve complementarity to your target. Do not use the result as a target-valid guide. Educational simulation, not biological optimization.")
        if st.button("🧠 Auto-Optimize Guide"):
            from .rl.dqn import benchmark, scorer
            with st.spinner("Comparing DQN, random and greedy on the same starting guide…"):
                model, _ = efficiency_resource()
                runs = benchmark(dqn_resource(), scorer(model), guide)
            st.session_state["priority_rl"] = (guide, runs)
        saved = st.session_state.get("priority_rl")
        if saved and saved[0] == guide:
            from .presentation import table
            runs = saved[1]
            results["optimization"] = runs
            st.code(runs[0]["final_guide"], language=None)
            table(runs[0]["steps"])
            st.markdown("**Benchmark · actual results, one start / seed 42**")
            table([{k: r[k] for k in ("policy", "before", "after")} for r in runs])
            st.caption("Not a population benchmark. Same 10-step budget, unequal compute: greedy evaluates 80 candidates/step. No promise DQN wins; no target-binding or safety constraint.")
        export_center({"sequence": sequence, "guide": guide, "pam": pam, "sequence_source": source}, results,
                      "User input / bundled library; source metadata per result; simulated ML/RL", "guide_tools")
