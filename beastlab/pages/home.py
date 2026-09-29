"""Mission control — the landing page and the map of the laboratory."""
from __future__ import annotations


import streamlit as st
from ..presentation import table, chart

from .. import theme, viz
from ..data.cas_enzymes import (BASE_EDITORS, NUCLEASES, PRIME_EDITORS,
                                RNA_EFFECTORS)
from ..data.clinical import TRIALS
from ..data.microbes import BACTERIA, PHAGES
from ..data.targets import all_loci


def _total_real_bp() -> int:
    return sum(l["length"] for l in all_loci())


def render() -> None:
    st.markdown(theme.hero(
        "🧬 BEAST Simulation Laboratory",
        "A digital genome-surgery suite: run CRISPR immunity battles inside bacteria, "
        "design guides against real human and bacterial sequence, then watch how a cell "
        "repairs — or fails to repair — the cut.",
        ["real sequences", "8 workbenches", "no black boxes"],
    ), unsafe_allow_html=True)

    loci = all_loci()
    st.markdown(theme.metrics([
        ("Cas effectors", str(len(NUCLEASES) + len(RNA_EFFECTORS)),
         "DNA + RNA targeting, natural and engineered"),
        ("Editors", str(len(BASE_EDITORS) + len(PRIME_EDITORS)),
         "cytosine / adenine base editors, prime editors"),
        ("Hosts", str(len(BACTERIA)), "bacteria + archaea with real CRISPR systems"),
        ("Phages", str(len(PHAGES)), "invaders, one with a fully embedded genome"),
        ("Genome library", f"{_total_real_bp() / 1000:.1f} kb",
         f"{len(loci)} real loci (NCBI / Ensembl GRCh38)"),
        ("Trials on file", str(len(TRIALS)), "approved and in-flight programmes"),
    ]), unsafe_allow_html=True)

    st.markdown("### The two arenas")

    left, right = st.columns(2)
    with left:
        st.markdown(theme.card(
            "🦠 Immunity arena — the bacterial arms race",
            "Bacteria keep a molecular photo album of every phage that ever attacked them. "
            "Challenge a real host genome with a real phage genome: watch Cas1-Cas2 write a "
            "new spacer into the array at the leader end, watch Cascade find the target, and "
            "then watch the phage mutate its way out. Population dynamics included.",
            "teal"), unsafe_allow_html=True)
        st.markdown(viz.svg_bacterium(immune=True, array_spacers=7, new_spacers=2,
                                      label="E. coli K-12 · type I-E · new spacer acquired"),
                    unsafe_allow_html=True)
    with right:
        st.markdown(theme.card(
            "✂️ Genome surgery — cut, edit, repair",
            "Pick a cell type, a delivery route and an editing modality, then run the "
            "experiment: DSB repair competes between c-NHEJ, MMEJ and HDR, indels appear, "
            "base editors convert single letters without cutting, prime editors write new "
            "sequence at the nick. Allele tables, a simulated gel and a safety radar at the end.",
            "violet"), unsafe_allow_html=True)
        st.markdown(viz.svg_human_cell(edited_percent=62, label="CD34+ HSPC · BCL11A enhancer"), unsafe_allow_html=True)

    st.markdown("### How to drive it")
    steps = st.columns(4)
    guide = [
        ("1 · Choose a target", "Guide design", "Point the bench at a real locus — sickle HBB, "
         "the BCL11A enhancer, PCSK9, TTR, CEP290, the E. coli CRISPR array — or paste your own sequence."),
        ("2 · Choose a weapon", "Cas explorer", "NGG-cutting SpCas9 or a near-PAMless SpRY, "
         "a staggered-cut Cas12a, a 529 aa miniature Cas12f, or an editor that never cuts at all."),
        ("3 · Run the experiment", "Genome surgery", "Set the cell type, the dose and the "
         "delivery route, then read the allele table, the gel and the collateral-damage radar."),
        ("4 · Understand the repair", "Repair lab", "Open the decision tree: why HDR barely "
         "happens in a hepatocyte, why microhomology decides the size of the deletion."),
    ]
    for col, (title, page, body) in zip(steps, guide):
        with col:
            st.markdown(
                f'<div class="card"><h4>{theme.badge("step", theme.C["gold"])} {title}</h4>'
                f'<p>{body}</p><p style="margin-top:8px;color:{theme.C["teal"]};'
                f'font-weight:700">→ {page}</p></div>', unsafe_allow_html=True)

    st.markdown("### What is real, and what is modelled")

    a, b = st.columns(2)
    with a:
        st.markdown(theme.card(
            "Measured / retrieved data",
            "• Every genome in the library: NCBI RefSeq or Ensembl GRCh38, with coordinates printed "
            "on screen.<br>"
            "• The E. coli K-12 CRISPR array is parsed out of the real chromosome window at run time.<br>"
            "• PAM sequences, spacer lengths, cut geometries and protein sizes of every effector.<br>"
            "• Repeat sequences of S. thermophilus CRISPR1/2/3, S. pyogenes, E. coli, P. aeruginosa.<br>"
            "• Clinical programme numbers: patient counts, doses, response rates.",
            "green"), unsafe_allow_html=True)
    with b:
        st.markdown(theme.card(
            "Simulation / heuristic",
            "• The on-target guide score is a transparent heuristic — not Doench Rule Set 2.<br>"
            "• Indel sizes are drawn from an empirical spectrum; microhomology deletions *are* "
            "computed from sequence.<br>"
            "• Phage-bacteria dynamics use textbook parameters; the arms race is qualitative.<br>"
            "• Patient cohort projections illustrate dose-response logic, nothing more.<br>"
            "• Off-target scoring uses the published MIT/Hsu position weights.",
            "plasma"), unsafe_allow_html=True)

    st.markdown("### Pick a question")
    st.caption("A quick router — every suggestion opens a live simulation, not a text answer.")

    q = st.radio("What do you want to find out?",
                 ["Why does a phage sometimes beat CRISPR?",
                  "Which Cas enzyme should I use, and can it even target my site?",
                  "Why is precise repair (HDR) so hard, and how do I get more of it?",
                  "How do you cure sickle-cell disease without a donor?",
                  "What does base editing do that nucleases cannot?"],
                 horizontal=False, label_visibility="collapsed")
    advice = {
        "Why does a phage sometimes beat CRISPR?": (
            "Immunity arena", "Run the challenge: escape mutations in the protospacer or PAM, "
            "plus anti-CRISPR proteins from real phages, are the two ways a phage wins. "
            "The escape playground lets you mutate the target base by base and watch immunity fail."),
        "Which Cas enzyme should I use, and can it even target my site?": (
            "Cas explorer → Guide design", "Compare PAM density (NGG every ~8-16 bp vs near-PAMless "
            "SpRY vs a T-rich Cas12a PAM), then design and rank guides with real off-target scans."),
        "Why is precise repair (HDR) so hard, and how do I get more of it?": (
            "Repair lab", "Open the decision tree, then try the levers: synchronise the cells, "
            "tether or protect the donor, inhibit 53BP1, and watch the HDR slice grow."),
        "How do you cure sickle-cell disease without a donor?": (
            "Therapeutics → Genome surgery", "Casgevy does not repair HBB at all: it breaks the "
            "BCL11A enhancer so fetal haemoglobin returns. Simulate the editing fraction and see "
            "what HbF and crisis rate follow."),
        "What does base editing do that nucleases cannot?": (
            "Genome surgery (modality = base editing)", "Convert a single base with no double-strand "
            "break: far fewer indels and large deletions, but you are locked into A>G or C>T and "
            "into a 4-8 nt window with bystander positions."),
    }
    page, text = advice[q]
    st.markdown(theme.callout(f"<b>{page}</b> — {text}", "info"), unsafe_allow_html=True)

    st.markdown("### The vocabulary, in one breath")
    st.markdown(theme.kv({
        "PAM": "a short motif next to the protospacer that the effector must see before it engages",
        "Protospacer": "the sequence in the invader / genome that matches the guide",
        "Spacer": "the same sequence as stored in the host's CRISPR array",
        "crRNA / sgRNA": "the guide RNA; sgRNA = crRNA fused to the tracrRNA scaffold",
        "Cascade": "the multi-protein type I surveillance complex (no Cas9 involved)",
        "Cas3": "the helicase-nuclease Cascade recruits to chew up the target",
        "R-loop": "the RNA:DNA hybrid invasion that tests a candidate target",
        "Seed": "the PAM-proximal guide nucleotides that dominate specificity",
        "NHEJ / MMEJ / HDR": "the three fates of a double-strand break",
    }), unsafe_allow_html=True)
