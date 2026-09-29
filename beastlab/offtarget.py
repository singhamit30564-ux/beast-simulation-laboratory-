"""Bounded public reference fetching; exact adjacent PAM, both-strand near matches."""
import json
import hashlib
from pathlib import Path
from urllib.request import Request, urlopen
import streamlit as st
from .sequtils import revcomp, parse_fasta

REGIONS = {"HBB": ("11", 5225464, 5227071), "BCL11A enhancer": ("2", 60495000, 60495500)}


@st.cache_data(ttl=86400, max_entries=4, show_spinner=False)
def fetch_reference(name):
    chrom, start, end = REGIONS[name]
    url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1?coord_system_version=GRCh38"
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "BEAST-education/2"})
    with urlopen(req, timeout=8) as response:
        obj = json.loads(response.read(100000))
    seq = obj["seq"].upper()
    if len(seq) != end - start + 1 or set(seq) - set("ACGTN"):
        raise ValueError("Reference response failed length/alphabet validation")
    return seq, {"source": "Ensembl REST", "assembly": "GRCh38", "chromosome": chrom,
                 "start": start, "end": end, "url": url, "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest()}


def reference_or_fallback(name):
    try:
        seq, meta = fetch_reference(name)
        return seq, meta, None
    except (OSError, ValueError, KeyError) as exc:
        path = Path(__file__).resolve().parents[1] / "data/genomes/phix174.fasta"
        seq = parse_fasta(path.read_text())[0][1]
        return seq, {"source": "Bundled NCBI phiX174 (offline demonstration, NOT human)",
                     "assembly": "NC_001422.1", "chromosome": "phiX174", "start": 1, "end": len(seq)}, type(exc).__name__


def scan(seq, guide, *, offset=1, chromosome="local", max_mm=3):
    """20-nt SpCas9 guides; includes exact matches (not presumed intended)."""
    from .io import parse_dna
    guide = parse_dna(guide, guide=True)
    hits = []
    for strand, strand_seq in (("+", seq.upper()), ("-", revcomp(seq))):
        for i in range(len(strand_seq) - 22):
            site, pam = strand_seq[i:i+20], strand_seq[i+20:i+23]
            if pam[1:] != "GG" or set(site + pam) - set("ACGT"):
                continue
            mm = [j for j in range(20) if site[j] != guide[j]]
            if len(mm) > max_mm:
                continue
            start = i if strand == "+" else len(seq) - i - 20
            seed = sum(j >= 8 for j in mm)  # last 12 bases, PAM-proximal
            hits.append({"position": f"{chromosome}:{offset+start}-{offset+start+19}",
                         "strand": strand, "mismatches": len(mm), "seed mismatches (last 12)": seed,
                         "PAM": pam, "site": site,
                         "risk (heuristic)": "high" if seed == 0 else "moderate" if seed == 1 else "lower",
                         "identity": "exact match; intended site not assigned" if not mm else "near match"})
    return sorted(hits, key=lambda row: (row["mismatches"], row["seed mismatches (last 12)"]))
