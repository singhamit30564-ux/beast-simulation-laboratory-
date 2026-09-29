"""Strict, bounded user input and in-memory, hash-verifiable exports."""
import csv
import hashlib
import io
import json

ASSUMPTIONS = "Educational simulation; sequence-only models and declared toy assumptions."
LIMITATIONS = "Not clinical advice, experimental validation, a whole-genome safety screen, or proof of cure."


def parse_dna(text, *, guide=False, max_length=10000):
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        raise ValueError("Enter a sequence.")
    if lines[0].startswith(">"):
        lines = lines[1:]
    if any(line.startswith(">") for line in lines):
        raise ValueError("Use one FASTA record only.")
    seq = "".join("".join(lines).split()).upper()
    if guide:
        seq = seq.replace("U", "T")
    if not seq or set(seq) - set("ACGT"):
        raise ValueError("Only A, C, G, T DNA is accepted; RNA U is allowed in guides. No ambiguous bases or digits.")
    if len(seq) > max_length:
        raise ValueError(f"Maximum {max_length:,} bases for phone-friendly computation.")
    if guide and len(seq) != 20:
        raise ValueError("SpCas9 guide must be exactly 20 bases, without PAM.")
    return seq


def canonical(payload):
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def export_bundle(inputs, results, source):
    payload = {"schema": "beastlab.v2", "inputs": inputs, "results": results,
               "source": source, "assumptions": ASSUMPTIONS, "what_this_does_not_claim": LIMITATIONS}
    digest = hashlib.sha256(canonical(payload)).hexdigest()
    envelope = {"sha256": digest, "hash_scope": "canonical payload JSON (sorted keys, compact separators, ASCII)",
                "payload": payload}
    js = json.dumps(envelope, indent=2).encode()
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(["sha256", "hash_scope", "payload_json"])
    writer.writerow([digest, envelope["hash_scope"], canonical(payload).decode()])
    from reportlab.pdfgen import canvas
    from reportlab.lib.utils import simpleSplit
    output = io.BytesIO()
    pdf = canvas.Canvas(output, pagesize=(595, 842))
    pdf.setTitle("BEAST Laboratory - educational simulation")
    y = 805
    for line in ("BEAST Laboratory | Educational simulation\n" + json.dumps(envelope, indent=2)).splitlines():
        for part in simpleSplit(line, "Courier", 8, 515):
            if y < 40:
                pdf.showPage()
                y = 805
            pdf.setFont("Courier", 8)
            pdf.drawString(40, y, part)
            y -= 11
    pdf.save()
    return {"json": js, "csv": stream.getvalue().encode(), "pdf": output.getvalue()}


def export_center(inputs, results, source, key):
    import streamlit as st
    with st.expander("User I/O · Export center — PDF / CSV / JSON"):
        st.caption(ASSUMPTIONS + " What this does not claim: " + LIMITATIONS)
        st.caption("SHA-256 verifies the canonical payload, not scientific validity or authorship. CSV contains a lossless JSON payload; PDF prints the same payload and hash.")
        if st.button("Prepare exports in RAM", key=key + "_prepare"):
            st.session_state[key + "_exports"] = (canonical([inputs, results, source]), export_bundle(inputs, results, source))
        saved = st.session_state.get(key + "_exports")
        if saved and saved[0] == canonical([inputs, results, source]):
            for ext, data in saved[1].items():
                st.download_button(f"Download {ext.upper()}", data, f"beast-{key}.{ext}",
                                   mime={"json": "application/json", "csv": "text/csv", "pdf": "application/pdf"}[ext], key=key + ext)
        st.caption("Zero Data Retention — computed in RAM. Inputs remain in this live session until cleared or expired; no user sequences are written to disk or sent to APIs. Public reference sequences alone are cached.")
