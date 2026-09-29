"""Per-bench, seed-stable mascot rotation without global RNG state."""
import hashlib
import streamlit as st

NOTES = {
 "Mission control": ["Start with a question, not a button.", "A simulation is a map, not the village.", "Small screens deserve big questions."],
 "Immunity arena": ["Bacteria remember invaders in molecular snapshots.", "Defense has a cost; survival is a trade-off.", "A spacer is a memory, not a magic shield."],
 "Genome surgery": ["PAM is the doorman — no PAM, no entry.", "A cut is only the beginning; repair writes the ending.", "Two strands, one responsibility: check your assumptions."],
 "Repair lab": ["The cell holds the repair tools, not the scissors.", "An indel is not automatically a knockout.", "A repair distribution is not one cell's destiny."],
 "Guide design": ["A guide without off-target checking is a promise without proof.", "A higher score is a hypothesis, not proof.", "Changing a guide can change its target — never forget the match."],
 "Cas explorer": ["Different scissors need different doorways.", "Size matters when your delivery vehicle is small.", "No single Cas is best for every question."],
 "Therapeutics": ["HbF can help without correcting HBB itself.", "A classroom verdict is not a patient's prognosis.", "Delivery, safety and durability belong in every cure story."],
 "Database & sources": ["A source is the beginning of trust, not its end.", "Record the assembly before comparing coordinates.", "A hash proves consistency, not truth."],
}


def render(bench):
    key = "note_" + bench
    with st.container(border=True):
        st.markdown("### 🧑🏽‍🔬 Dr. Titan · Field notes")
        notes = NOTES.get(bench, NOTES["Mission control"])
        base = int(hashlib.sha256(f"42:{bench}".encode()).hexdigest()[:8], 16)
        index = st.session_state.get(key, 0)
        st.write(notes[(base + index) % len(notes)])
        if st.button("Next field note", key=key + "_next"):
            st.session_state[key] = index + 1
            st.rerun()
