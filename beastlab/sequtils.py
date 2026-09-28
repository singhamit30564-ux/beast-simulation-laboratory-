"""Low level nucleic-acid utilities.

Everything in this module is deterministic and dependency-free so that the
simulation core can be unit-tested without Streamlit, plotly or numpy.
"""
from __future__ import annotations

import math
import random
from typing import Dict, Iterator, List, Optional, Tuple

#: IUPAC nucleotide code -> set of concrete bases.
IUPAC: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
}

#: Standard genetic code (DNA codons).
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_1TO3: Dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Stop",
}


# --------------------------------------------------------------------------- #
# basic string operations
# --------------------------------------------------------------------------- #
def clean(seq: str) -> str:
    """Upper-case a nucleotide string, mapping U->T and dropping whitespace/digits."""
    return "".join(c for c in str(seq).upper().replace("U", "T") if c in IUPAC)


def revcomp(seq: str) -> str:
    """Reverse complement of a (possibly degenerate) nucleotide sequence."""
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(clean(seq)))


def gc_content(seq: str) -> float:
    """Fraction of G+C bases (0.0 for an empty sequence)."""
    s = clean(seq)
    if not s:
        return 0.0
    return sum(b in "GC" for b in s) / len(s)


def melting_temp(seq: str, na_conc_mM: float = 50.0) -> float:
    """Approximate duplex melting temperature (salt-adjusted, nearest-neighbour free).

    Uses the classic Wallace rule for short oligos (<= 13 nt) and the
    GC%-based formula of Marmur-Doty otherwise.  Indicative only — the app
    quotes it for PBS / donor-arm design sanity checks.
    """
    s = clean(seq)
    n = len(s)
    if n == 0:
        return 0.0
    gc = sum(b in "GC" for b in s)
    if n <= 13:
        return 2.0 * (n - gc) + 4.0 * gc
    salt_corr = 16.6 * math.log10(max(na_conc_mM, 1) / 1000.0) + 16.6
    return 64.9 + 41.0 * (gc - 16.4) / n + salt_corr - 16.6


def hamming(a: str, b: str) -> int:
    """Number of positions at which two equal-length strings differ."""
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def mismatch_positions(a: str, b: str) -> List[int]:
    """1-based positions at which two aligned sequences differ."""
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def wrap(seq: str, width: int = 60) -> str:
    """Hard-wrap a sequence for display / FASTA output."""
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def longest_homopolymer(seq: str) -> Tuple[str, int]:
    """Return the longest run of a single base (base, length)."""
    best, cur, best_base = 1, 1, ""
    s = clean(seq)
    if not s:
        return "", 0
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            cur += 1
        else:
            cur = 1
        if cur > best:
            best, best_base = cur, s[i]
    return (best_base or s[0]), best


# --------------------------------------------------------------------------- #
# motif searching
# --------------------------------------------------------------------------- #
def _expand(pattern: str) -> List[str]:
    """Expand a short IUPAC pattern into all concrete sequences."""
    pools = [IUPAC[c] for c in clean(pattern)]
    out: List[str] = [""]
    for pool in pools:
        out = [prefix + b for prefix in out for b in pool]
    return out


def find_motifs(seq: str, pattern: str, both_strands: bool = False) -> List[Tuple[int, str]]:
    """Find IUPAC `pattern` in `seq`.

    Returns a list of ``(0-based start, matched_text)`` tuples.  With
    ``both_strands=True`` the reverse-complement of the pattern is searched too.
    """
    seq = clean(seq)
    targets = _expand(pattern)
    hits: List[Tuple[int, str]] = []
    seen = set()
    for t in targets:
        start = seq.find(t)
        while start != -1:
            if (start, t) not in seen:
                seen.add((start, t))
                hits.append((start, t))
            start = seq.find(t, start + 1)
    if both_strands:
        rc = revcomp(pattern)
        for t in _expand(rc):
            start = seq.find(t)
            while start != -1:
                if (start, t) not in seen:
                    seen.add((start, t))
                    hits.append((start, t))
                start = seq.find(t, start + 1)
    return sorted(hits)


def pattern_regex(pattern: str) -> str:
    """Compile an IUPAC pattern into a regular-expression string."""
    return "".join(f"[{IUPAC[c]}]" for c in clean(pattern))


# --------------------------------------------------------------------------- #
# translation and reading frames
# --------------------------------------------------------------------------- #
def translate(seq: str, frame: int = 0, stop_symbol: str = "*") -> str:
    """Translate a nucleotide sequence in the given reading frame (0/1/2)."""
    s = clean(seq)[frame:]
    return "".join(
        CODON_TABLE.get(s[i:i + 3], "X") if len(s[i:i + 3]) == 3 else ""
        for i in range(0, len(s) - len(s) % 3, 3)
    ).replace("*", stop_symbol)


def orfs(seq: str, min_aa: int = 30) -> List[Tuple[int, int, str]]:
    """Very small six-frame ORF finder -> (start, end, protein) on the plus strand."""
    s = clean(seq)
    out: List[Tuple[int, int, str]] = []
    for frame in range(3):
        aa = translate(s, frame, stop_symbol="*")
        start_aa: Optional[int] = None
        for i, res in enumerate(aa):
            if res == "M" and start_aa is None:
                start_aa = i
            if res == "*" and start_aa is not None:
                if i - start_aa >= min_aa:
                    out.append((frame + start_aa * 3, frame + i * 3 + 3, aa[start_aa:i]))
                start_aa = None
    return out


# --------------------------------------------------------------------------- #
# reading / writing sequence records
# --------------------------------------------------------------------------- #
def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Parse FASTA text into ``[(header, sequence), ...]``."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(clean(line))
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def to_fasta(header: str, seq: str, width: int = 70) -> str:
    """Serialise one record as FASTA text."""
    return f">{header}\n{wrap(clean(seq), width)}\n"


# --------------------------------------------------------------------------- #
# sequence surgery helpers
# --------------------------------------------------------------------------- #
def mutate(seq: str, edits: Dict[int, str]) -> str:
    """Apply 1-based substitutions ``{position: base}`` and return the new sequence."""
    s = list(clean(seq))
    for pos, base in edits.items():
        if 1 <= pos <= len(s):
            s[pos - 1] = base.upper()
    return "".join(s)


def microhomologies(seq: str, cut: int, min_len: int = 3, max_len: int = 12,
                    window: int = 40) -> List[Tuple[int, int, int, str]]:
    """Find microhomologies flanking a cut site (MMEJ substrate search).

    Returns ``(left_start, right_start, length, repeat)`` with 0-based
    coordinates on `seq`.  Only perfect direct repeats on the *same* strand are
    reported, which is the class of repeat that MMEJ can anneal.
    """
    s = clean(seq)
    lo, hi = max(0, cut - window), min(len(s), cut + window)
    left_zone = s[lo:cut]
    right_zone = s[cut:hi]
    found: List[Tuple[int, int, int, str]] = []
    for length in range(max_len, min_len - 1, -1):
        for i in range(len(left_zone) - length + 1):
            repeat = left_zone[i:i + length]
            j = right_zone.find(repeat)
            if j != -1:
                found.append((lo + i, cut + j, length, repeat))
    # keep longest repeat per left start, sorted by repeat length
    best: Dict[int, Tuple[int, int, int, str]] = {}
    for hit in found:
        if hit[0] not in best or hit[2] > best[hit[0]][2]:
            best[hit[0]] = hit
    return sorted(best.values(), key=lambda h: (-h[2], h[0]))


def random_sequence(n: int, gc: float = 0.5, seed: Optional[int] = None) -> str:
    """Random DNA of length `n` with target GC content."""
    rng = random.Random(seed)
    return "".join(
        rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n)
    )


def complexity(seq: str, k: int = 4) -> float:
    """Simple k-mer complexity score (unique k-mers / possible k-mers)."""
    s = clean(seq)
    if len(s) < k:
        return 0.0
    kmer = {s[i:i + k] for i in range(len(s) - k + 1)}
    return len(kmer) / max(1, len(s) - k + 1)


def slice_with_context(seq: str, start: int, end: int, context: int = 10,
                       upper: bool = True) -> str:
    """Return ``seq[start:end]`` flanked by lowercase context (0-based, half-open)."""
    lo, hi = max(0, start - context), min(len(seq), end + context)
    left, core, right = seq[lo:start], seq[start:end], seq[end:hi]
    if not upper:
        return left.lower() + core.lower() + right.lower()
    return left.lower() + core.upper() + right.lower()


def chunk_positions(length: int, chunk: int = 10) -> Iterator[Tuple[int, int]]:
    """Yield (start, end) blocks covering ``range(length)``."""
    for i in range(0, length, chunk):
        yield i, min(i + chunk, length)


def kmer_counts(seq: str, k: int = 2) -> Dict[str, int]:
    """Count k-mers (used for dinucleotide/codon usage summaries)."""
    s = clean(seq)
    counts: Dict[str, int] = {}
    for i in range(len(s) - k + 1):
        word = s[i:i + k]
        counts[word] = counts.get(word, 0) + 1
    return counts
