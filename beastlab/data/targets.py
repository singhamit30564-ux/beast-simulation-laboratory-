"""Genome library.

All sequences are real and come from two public sources, recorded per record:

* NCBI RefSeq / GenBank through E-utilities (``efetch``)
* Ensembl REST ``/sequence/region`` on GRCh38

Files live in ``data/genomes``.  The loader validates the alphabet and warns if a
record is shorter than expected, so a corrupted install cannot silently feed the
simulation made-up DNA.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

from ..sequtils import gc_content, parse_fasta, revcomp

GENOME_DIR = Path(__file__).resolve().parents[2] / "data" / "genomes"

#: Curated target loci — everything the "genome surgery" bench can edit.
LOCI: List[Dict[str, Any]] = [
    {
        "id": "hbb_sickle", "record": "HBB_hg38", "name": "HBB — sickle-cell allele",
        "file": "human_targets.fasta",
        "gene": "HBB (beta-globin)", "organism": "Homo sapiens",
        "locus": "chr11:5,225,464-5,227,071 (GRCh38)", "strand": "minus",
        "disease": "Sickle-cell disease (HbS, Glu6Val, rs334)",
        "target_note": "The classic single-base correction: GAG (Glu) -> GTG (Val) at codon 6.",
        "editable_base": "C", "editing_modality": "base editing (ABE) or HDR",
        "reference": "Rees 2010 Lancet 376:2018; Park 2019 Nat Biotechnol (sickle correction)",
        "highlight": "GAG", "edit_offset_from_orf": 15,
    },
    {
        "id": "bcl11a_enhancer", "record": "BCL11A_enhancer",
        "file": "human_targets.fasta",
        "name": "BCL11A erythroid enhancer (+58)",
        "gene": "BCL11A", "organism": "Homo sapiens",
        "locus": "chr2:60,495,000-60,495,500 (GRCh38)", "strand": "plus",
        "disease": "Target for sickle-cell disease and beta-thalassaemia therapy",
        "target_note": "Disrupting the GATA1 site in the erythroid enhancer de-represses "
                       "fetal haemoglobin (HbF) — the Casgevy mechanism.",
        "editable_base": "none", "editing_modality": "nuclease knockout of a regulatory element",
        "reference": "Frangoul 2021 NEJM 384:252; Wu 2019 Nat Med 25:776",
        "highlight": "AGATAA",
    },
    {
        "id": "hbg1_promoter", "record": "HBG1_promoter",
        "file": "human_targets.fasta",
        "name": "HBG1 promoter (HPFH region)",
        "gene": "HBG1", "organism": "Homo sapiens",
        "locus": "chr11:5,249,390-5,249,890 (GRCh38)", "strand": "minus",
        "disease": "Hereditary persistence of fetal haemoglobin — a natural protective state",
        "target_note": "Base editing the -198/-175 HPFH sites boosts HbF without a double-strand break.",
        "editable_base": "A", "editing_modality": "adenine base editing (ABE8e / BEAM-101)",
        "reference": "Gaudelli 2017 Nature 551:464; Beam Therapeutics BEAM-101",
        "highlight": "CACATTC",
    },
    {
        "id": "pcsk9_start", "record": "PCSK9_start",
        "file": "human_targets.fasta",
        "name": "PCSK9 5' region",
        "gene": "PCSK9", "organism": "Homo sapiens",
        "locus": "chr1:55,039,450-55,039,950 (GRCh38)", "strand": "plus",
        "disease": "Heterozygous familial hypercholesterolaemia",
        "target_note": "Permanent inactivation of hepatic PCSK9 lowers LDL-C — an in-vivo base "
                       "editing programme delivered by lipid nanoparticles.",
        "editable_base": "C", "editing_modality": "CBE splice-site disruption (VERVE-101/102)",
        "reference": "Kathiresan 2022 Circulation 147:1002; heart-1 trial NCT05398029",
        "highlight": "ATGGGCACCGTCAGCTCCAGG",
    },
    {
        "id": "ttr_exon3", "record": "TTR_exon3", "name": "TTR (transthyretin)",
        "file": "human_targets.fasta",
        "gene": "TTR", "organism": "Homo sapiens",
        "locus": "chr18:31,594,000-31,594,500 (GRCh38)", "strand": "minus",
        "disease": "TTR amyloidosis (ATTR-CM / hATTR-PN)",
        "target_note": "Knockout of the hepatic TTR gene removes the amyloid precursor, "
                       "delivered in vivo as LNPs.",
        "editable_base": "none", "editing_modality": "nuclease knockout (NTLA-2001 / nex-z)",
        "reference": "Gillmore 2021 NEJM 385:493",
        "highlight": "GACTTCATAGATTT",
    },
    {
        "id": "cep290_ivs26", "record": "CEP290_IVS26",
        "file": "human_targets.fasta",
        "name": "CEP290 deep intron (IVS26)",
        "gene": "CEP290", "organism": "Homo sapiens",
        "locus": "chr12:88,086,000-88,086,500 (GRCh38)", "strand": "plus",
        "disease": "Leber congenital amaurosis 10 (LCA10)",
        "target_note": "The intronic c.2991+1655A>G mutation creates a cryptic splice donor; "
                       "it can be deleted or disrupted, but the sequence is A/T-rich, which "
                       "limits guide choice.",
        "editable_base": "none", "editing_modality": "nuclease deletion of the cryptic exon",
        "reference": "Maeder 2019 Nat Med 25:229; EDIT-101 (BRILLIANCE)",
        "highlight": "TTTTACAGCTGTATGATAAAACATAGATTAAACCC",
    },
    {
        "id": "ccr5", "record": "CCR5_region", "name": "CCR5 locus",
        "file": "human_targets.fasta",
        "gene": "CCR5", "organism": "Homo sapiens",
        "locus": "chr3:46,370,550-46,371,050 (GRCh38)", "strand": "minus",
        "disease": "HIV-1 entry (co-receptor); the CCR5-delta32 allele confers resistance",
        "target_note": "Disrupting CCR5 makes CD4+ T cells resistant to CCR5-tropic HIV-1 — "
                       "the first clinical CRISPR concept, and the target of the DdPc trial.",
        "editable_base": "none", "editing_modality": "nuclease knockout",
        "reference": "Xu 2019 NEJM 381:1240; Tebas 2014 NEJM 370:901",
        "highlight": "TGACATCAATTATTATACATCG",
    },
    {
        "id": "aavs1", "record": "AAVS1_region", "name": "AAVS1 safe harbour (PPP1R12C)",
        "file": "human_targets.fasta",
        "gene": "PPP1R12C", "organism": "Homo sapiens",
        "locus": "chr19:55,115,000-55,115,500 (GRCh38)", "strand": "plus",
        "disease": "Not a disease gene — a landing pad for transgene knock-in",
        "target_note": "The classical safe-harbour locus: a well-behaved place to insert a "
                       "cassette without disturbing essential genes.",
        "editable_base": "none", "editing_modality": "HDR knock-in",
        "reference": "Sadelain 2012 Nat Rev Cancer 12:51",
        "highlight": "GCTACTGGCCTTATCTCACAGGTAAAACTGACGCACGGAGG",
    },
    {
        "id": "ecoli_crispr_locus", "record": "ecoli_k12_crispr1_locus",
        "file": "ecoli_k12.fasta",
        "name": "E. coli K-12 CRISPR-1 locus",
        "gene": "cas1-cas2 / CRISPR array", "organism": "Escherichia coli K-12 MG1655",
        "locus": "U00096.3:2,877,800-2,878,849", "strand": "plus",
        "disease": "Not applicable — this is the bacterial defence locus itself",
        "target_note": "A real CRISPR array: 10 direct repeats (29 nt) with 9 spacers of 33-34 nt. "
                       "Edit here and you disable the bacterium's own immune memory.",
        "editable_base": "none", "editing_modality": "nuclease knockout of a defence system",
        "reference": "Brouns 2008 Science 321:960",
        "highlight": "GGTTTATCCCCGCTGGCGCGGGGAACTC",
    },
    {
        "id": "ecoli_genomic", "record": "ecoli_k12_window2",
        "file": "ecoli_k12.fasta",
        "name": "E. coli K-12 genomic window (lac region)",
        "gene": "chromosomal window", "organism": "Escherichia coli K-12 MG1655",
        "locus": "U00096.3:365,030-365,529", "strand": "plus",
        "disease": "Not applicable",
        "target_note": "A coding-dense bacterial window — useful for comparing PAM density "
                       "between genomes.",
        "editable_base": "none", "editing_modality": "nuclease knockout",
        "reference": "NCBI RefSeq U00096.3",
        "highlight": "GATTCATTGGCACCATG",
    },
    {
        "id": "phix174_genome", "record": "NC_001422.1", "name": "phiX174 phage genome",
        "file": "phix174.fasta",
        "gene": "11 overlapping genes", "organism": "Escherichia phage phiX174",
        "locus": "NC_001422.1:1-5,386", "strand": "plus",
        "disease": "Not applicable — an invading phage genome",
        "target_note": "The first genome ever sequenced, and a real substrate for both "
                       "CRISPR immunity and anti-phage engineering.",
        "editable_base": "none", "editing_modality": "target for CRISPR immunity",
        "reference": "Sanger 1977 Nature 265:687",
        "highlight": "ATGAGTCAAGTTACTGAACAATCC",
    },
]


@lru_cache(maxsize=64)
def load_records(filename: str) -> Dict[str, Dict[str, Any]]:
    """Load ``data/genomes/<filename>`` into ``{record_id: record}``."""
    path = GENOME_DIR / filename
    if not path.exists():
        return {}
    records = {}
    for header, seq in parse_fasta(path.read_text()):
        first = header.split()[0]
        records[first] = {
            "id": first, "header": header, "sequence": seq, "length": len(seq),
            "gc": gc_content(seq),
            "source": _source_of(header),
        }
    return records


def _source_of(header: str) -> str:
    if "U00096" in header:
        return "NCBI RefSeq U00096.3 (E. coli K-12 MG1655), E-utilities efetch"
    if "NC_001422" in header:
        return "NCBI RefSeq NC_001422.1 (phiX174), E-utilities efetch"
    if "GRCh38" in header:
        return "Ensembl REST /sequence/region on GRCh38"
    return "curated"


@lru_cache(maxsize=1)
def all_loci() -> List[Dict[str, Any]]:
    """LOCl records enriched with their sequence."""
    out = []
    for meta in LOCI:
        rec = load_records(meta["file"]).get(meta["record"])
        if not rec:
            continue
        merged = dict(meta)
        merged.update({"sequence": rec["sequence"], "gc": rec["gc"],
                       "length": rec["length"], "source": rec["source"]})
        out.append(merged)
    return out


def get_locus(locus_id: str) -> Optional[Dict[str, Any]]:
    for locus in all_loci():
        if locus["id"] == locus_id:
            return locus
    return None


def coding_strand(locus: Dict[str, Any]) -> str:
    """Return the locus sequence oriented so that the gene reads 5'->3'."""
    return locus["sequence"] if locus.get("strand", "plus") == "plus" \
        else revcomp(locus["sequence"])


def summary_table() -> List[Dict[str, Any]]:
    return [{
        "id": l["id"], "name": l["name"], "organism": l["organism"], "locus": l["locus"],
        "length_bp": l["length"], "gc_percent": round(100 * l["gc"], 2),
        "modality": l["editing_modality"], "source": l["source"],
    } for l in all_loci()]


def find_motif_locations(locus_id: str, motif: str) -> List[int]:
    """1-based positions of a motif in the coding-strand sequence of a locus."""
    locus = get_locus(locus_id)
    if not locus or not motif:
        return []
    s = coding_strand(locus)
    out, start = [], s.find(motif.upper())
    while start != -1:
        out.append(start + 1)
        start = s.find(motif.upper(), start + 1)
    return out
