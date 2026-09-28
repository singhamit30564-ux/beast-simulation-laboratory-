"""Bacteria, archaea, their CRISPR-Cas systems and the phages that hunt them.

Repeat sequences are marked with ``repeat_source``:

* ``literature``  – exact repeat published in the cited paper
* ``sequenced``   – taken directly from the genome window shipped in ``data/genomes``
* ``type-level``  – no repeat sequence is claimed; only the system type/PAM, which
  is what the record actually documents

The UI surfaces this field so that nothing is passed off as curated when it is
not.
"""
from __future__ import annotations

from typing import Any, Dict, List

BACTERIA: List[Dict[str, Any]] = [
    {
        "id": "ecoli_k12",
        "name": "Escherichia coli K-12 MG1655",
        "short": "E. coli K-12",
        "domain": "Bacteria",
        "phylum": "Pseudomonadota (γ-proteobacteria)",
        "gram": "negative",
        "morphology": "rod, ~2 µm x 0.5 µm",
        "habitat": "lower gut of warm-blooded animals; the workhorse of molecular biology",
        "genome_mb": 4.64, "gc": 0.508, "genes": 4140,
        "accession": "U00096.3",
        "systems": [
            {
                "name": "CRISPR1 (type I-E, CASS2)",
                "type": "I", "subtype": "I-E",
                "cas_genes": ["cas3", "cse1(cas8e/casA)", "cse2(casB)", "cas7", "cas5", "cas6e", "cas1", "cas2"],
                "effector": "Cascade (Cas8e-Cse2-Cas7₆-Cas5-Cas6e) + Cas3 helicase-nuclease",
                "pam": "5'-AAG-3' (protospacer flanked by 3'-TTC-5')",
                "pam_iupac": "AAG", "pam_side": "5",
                "protospacer_len": 32,
                "repeat": "GGTTTATCCCCGCTGGCGCGGGGAACTC",
                "repeat_len": 28,
                "repeat_source": "sequenced",
                "array_file": "ecoli_k12.fasta",
                "array_region": "U00096.3:2,877,800-2,878,849",
                "spacers_observed": 9,
                "activity": "silenced by H-NS at 37 °C; active at 30 °C — a textbook example of a "
                            "defence system that evolution has switched down rather than deleted",
                "notes": "No Cas9. Interference is a Cascade/Cas3 R-loop that is degraded processively "
                         "from the PAM-proximal end.",
                "ref": "Brouns 2008 Science 321:960; Jore 2011 Nat Struct Mol Biol 18:529; "
                       "Datsenko 2012 PNAS 109:E2579",
            }
        ],
        "interesting": "The K-12 lab strain carries a CRISPR-1 array with only a handful of spacers — "
                       "and it is transcriptionally repressed. E. coli is a bacterium that largely gave up on CRISPR.",
        "presentation": "gram-negative, no CRISPR defence active at 37 °C",
    },
    {
        "id": "sthermophilus",
        "name": "Streptococcus thermophilus DGCC7710",
        "short": "S. thermophilus",
        "domain": "Bacteria",
        "phylum": "Bacillota (Firmicutes)",
        "gram": "positive",
        "morphology": "ovoid cocci in chains, ~0.8 µm",
        "habitat": "dairy starter culture for yogurt and cheese",
        "genome_mb": 1.86, "gc": 0.390, "genes": 1900,
        "accession": "CP000419 / LMD-9",
        "systems": [
            {
                "name": "CRISPR1 (type II-A)",
                "type": "II", "subtype": "II-A",
                "cas_genes": ["cas9", "cas1", "cas2", "csn2"],
                "effector": "Cas9 + sgRNA (crRNA:tracrRNA)",
                "pam": "5'-NNAGAAW-3'",
                "pam_iupac": "NNAGAAW", "pam_side": "3",
                "protospacer_len": 30,
                "repeat": "GTTTTTGTACTCTCAAGATTTAAGTAACTGTACAAC",
                "repeat_len": 36, "repeat_source": "literature",
                "spacers_observed": 32,
                "activity": "the most active natural spacer-acquisition system known — "
                            "phage challenge produces bacteriophage-insensitive mutants (BIMs) within one passage",
                "notes": "Adaptation is polarized: new spacers are inserted at the leader-proximal end.",
                "ref": "Barrangou 2007 Science 315:1709; Horvath 2008 J Bacteriol 190:1401",
            },
            {
                "name": "CRISPR3 (type II-A)",
                "type": "II", "subtype": "II-A",
                "cas_genes": ["cas9", "cas1", "cas2", "csn2"],
                "effector": "Cas9 (St3Cas9) + sgRNA",
                "pam": "5'-NGGNG-3'",
                "pam_iupac": "NGGNG", "pam_side": "3",
                "protospacer_len": 30,
                "repeat": "GTTTTAGAGCTGTGTTGTTTCGAATGGTTCCAAAAC",
                "repeat_len": 36, "repeat_source": "literature",
                "spacers_observed": 8,
                "activity": "secondary acquisition locus; requires Csn2 for new spacer integration",
                "notes": "cas9 alone is sufficient for interference, but not for adaptation.",
                "ref": "Sapranauskas 2011 Nucleic Acids Res 39:9275; Garneau 2010 Nature 468:67",
            },
            {
                "name": "CRISPR2 (type III-A)",
                "type": "III", "subtype": "III-A",
                "cas_genes": ["cas10/csm1", "csm2", "csm3", "csm4", "csm5", "csm6", "cas1", "cas2", "cas6"],
                "effector": "Csm complex (targets transcriptionally active DNA, plus cOA signalling via Csm6)",
                "pam": "no strict PAM (transcription-dependent)",
                "pam_iupac": "", "pam_side": "none",
                "protospacer_len": 37,
                "repeat": "GATATAAACCTAATTACCTCGAGAGGGGACGGAAAC",
                "repeat_len": 36, "repeat_source": "literature",
                "spacers_observed": 5,
                "activity": "no spacer acquisition observed in the lab",
                "notes": "Type III systems cut the DNA of transcribed targets and also degrade RNA.",
                "ref": "Horvath 2008 J Bacteriol 190:1401",
            },
            {
                "name": "CRISPR4 (type I-E)",
                "type": "I", "subtype": "I-E",
                "cas_genes": ["cas3", "casA(cse1)", "casB(cse2)", "cas7", "cas5", "cas6e", "cas1", "cas2"],
                "effector": "St-Cascade + Cas3",
                "pam": "5'-AAG-3'-like (I-E family)",
                "pam_iupac": "AAG", "pam_side": "5",
                "protospacer_len": 33,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 6,
                "activity": "functional Cascade/Cas3 reconstituted in vitro",
                "notes": "The 61-nt crRNA bound by St-Cascade was the first type I crRNA structure solved.",
                "ref": "Sashital 2011; Mulepati & Bailey 2011; Jore 2011",
            },
        ],
        "interesting": "Yogurt is where CRISPR immunity was discovered: challenging S. thermophilus with "
                       "phages produced resistant mutants that had literally written pieces of the phage "
                       "genome into their own chromosome (Barrangou 2007).",
        "presentation": "gram-positive chain former, two type II-A systems plus I-E and III-A",
    },
    {
        "id": "spyogenes",
        "name": "Streptococcus pyogenes SF370",
        "short": "S. pyogenes",
        "domain": "Bacteria",
        "phylum": "Bacillota (Firmicutes)",
        "gram": "positive",
        "morphology": "cocci in chains, ~1 µm",
        "habitat": "human pharynx and skin — the source of SpCas9",
        "genome_mb": 1.85, "gc": 0.385, "genes": 1700,
        "accession": "NC_002737.2",
        "systems": [
            {
                "name": "CRISPR1 (type II-A) — origin of SpCas9",
                "type": "II", "subtype": "II-A",
                "cas_genes": ["cas9", "cas1", "cas2", "csn2"],
                "effector": "Cas9 + tracrRNA:crRNA",
                "pam": "5'-NGG-3'",
                "pam_iupac": "NGG", "pam_side": "3",
                "protospacer_len": 20,
                "repeat": "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC",
                "repeat_len": 36, "repeat_source": "literature",
                "spacers_observed": 6,
                "activity": "low in the lab; the array is short and largely static",
                "notes": "The tracrRNA of this locus was fused to the crRNA to make the first single-guide RNA.",
                "ref": "Jinek 2012 Science 337:816; Deltcheva 2011 Nature 471:602",
            }
        ],
        "interesting": "The single-guide RNA that made CRISPR easy is an artificial fusion of two RNAs "
                       "this bacterium uses naturally.",
        "presentation": "gram-positive pathogen, the birthplace of SpCas9",
    },
    {
        "id": "paeruginosa",
        "name": "Pseudomonas aeruginosa PA14",
        "short": "P. aeruginosa",
        "domain": "Bacteria",
        "phylum": "Pseudomonadota (γ-proteobacteria)",
        "gram": "negative",
        "morphology": "rod, ~2 µm, single polar flagellum",
        "habitat": "soil, water and hospitals — an opportunistic pathogen with big CRISPR arrays",
        "genome_mb": 6.54, "gc": 0.663, "genes": 5900,
        "accession": "NC_008463.1",
        "systems": [
            {
                "name": "CRISPR1 (type I-F)",
                "type": "I", "subtype": "I-F",
                "cas_genes": ["csy1", "csy2", "csy3", "cas6f(csy4)", "cas1", "cas2-cas3 fusion"],
                "effector": "Csy complex (9 subunits + crRNA) + Cas2/3 fusion helicase-nuclease",
                "pam": "5'-CC-3'",
                "pam_iupac": "CC", "pam_side": "5",
                "protospacer_len": 32,
                "repeat": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
                "repeat_len": 28, "repeat_source": "literature",
                "spacers_observed": 40,
                "activity": "constitutive and highly active; the first type I-F system shown to block phage infection",
                "notes": "P. aeruginosa phages fight back with anti-CRISPR (AcrIF) proteins.",
                "ref": "Cady 2011 J Bacteriol 193:3433; Bondy-Denomy 2013 Nature 493:429",
            }
        ],
        "interesting": "The battlefield where anti-CRISPR proteins were discovered: half of the "
                       "P. aeruginosa phages tested encode proteins that shut CRISPR-Cas down.",
        "presentation": "gram-negative opportunist with a 40-spacer array and crafty phages",
    },
    {
        "id": "m tuberculosis",
        "name": "Mycobacterium tuberculosis H37Rv",
        "short": "M. tuberculosis",
        "domain": "Bacteria",
        "phylum": "Actinomycetota",
        "gram": "positive (acid-fast, mycolic acid wall)",
        "morphology": "slow-growing rod, ~2 µm",
        "habitat": "human lung macrophages; the leading bacterial cause of death from a single agent",
        "genome_mb": 4.41, "gc": 0.657, "genes": 4000,
        "accession": "NC_000962.3",
        "systems": [
            {
                "name": "CRISPR1 (type III-A)",
                "type": "III", "subtype": "III-A",
                "cas_genes": ["cas10/csm1", "csm2", "csm3", "csm4", "csm5", "csm6", "cas1", "cas2", "cas6"],
                "effector": "Csm complex — transcription-dependent DNA interference + cOA signalling",
                "pam": "no strict PAM",
                "pam_iupac": "", "pam_side": "none",
                "protospacer_len": 36,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 0,
                "activity": "arrays are small and largely degenerate; clinical isolates often lack systems entirely",
                "notes": "Type III-A effectors recognise the nascent transcript of an actively "
                         "transcribed target — a defence tuned to expression, not just sequence.",
                "ref": "Makarova 2020 Nat Rev Microbiol 18:67; Grissa 2007 BMC Bioinformatics 8:172",
            }
        ],
        "interesting": "The CRISPR arrays of M. tuberculosis isolates were used as a genotyping "
                       "fingerprint (spoligotyping) long before the immune function was understood.",
        "presentation": "acid-fast pathogen, type III-A system, no strict PAM",
    },
    {
        "id": "n meningitidis",
        "name": "Neisseria meningitidis MC58",
        "short": "N. meningitidis",
        "domain": "Bacteria",
        "phylum": "Pseudomonadota (β-proteobacteria)",
        "gram": "negative",
        "morphology": "diplococcus, ~1 µm",
        "habitat": "human nasopharynx — also the source of NmeCas9",
        "genome_mb": 2.27, "gc": 0.515, "genes": 2100,
        "accession": "NC_003112.2",
        "systems": [
            {
                "name": "CRISPR1 (type II-C) — origin of Nme2Cas9",
                "type": "II", "subtype": "II-C",
                "cas_genes": ["cas9", "cas1", "cas2"],
                "effector": "Cas9 (Nme1/Nme2) + tracrRNA:crRNA",
                "pam": "5'-NNNNCC-3' (Nme2) / 5'-NNNNGATT-3' (Nme1)",
                "pam_iupac": "NNNNCC", "pam_side": "3",
                "protospacer_len": 22,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 12,
                "activity": "active; type II-C systems keep a short RNA duplex and are strictly guided",
                "notes": "Type II-C Cas9s are small, high-fidelity and recognise long PAMs — "
                         "excellent for allele-specific editing.",
                "ref": "Zhang 2013 Cell Rep 5:698; Esvelt 2013 Nat Methods 10:1116",
            }
        ],
        "interesting": "N. meningitidis Cas9s need an 8 nt PAM — less convenient, but far more "
                       "selective in a human genome.",
        "presentation": "gram-negative diplococcus, compact type II-C Cas9",
    },
    {
        "id": "cj jejuni",
        "name": "Campylobacter jejuni",
        "short": "C. jejuni",
        "domain": "Bacteria",
        "phylum": "Campylobacterota",
        "gram": "negative",
        "morphology": "curved rod, ~1.5 µm, microaerophilic",
        "habitat": "poultry gut — a leading cause of foodborne enteritis",
        "genome_mb": 1.64, "gc": 0.305, "genes": 1650,
        "accession": "NC_002163.1",
        "systems": [
            {
                "name": "CRISPR1 (type II-C) — origin of CjCas9",
                "type": "II", "subtype": "II-C",
                "cas_genes": ["cas9", "cas1", "cas2"],
                "effector": "Cas9 (984 aa) + tracrRNA:crRNA",
                "pam": "5'-NNNNRYAC-3'",
                "pam_iupac": "NNNNRYAC", "pam_side": "3",
                "protospacer_len": 22,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 4,
                "activity": "low GC genome, AT-rich; arrays change during poultry colonisation",
                "notes": "CjCas9 is the smallest widely used Cas9 (~984 aa) and therefore fits in AAV.",
                "ref": "Kim 2017 Nat Commun 8:14500",
            }
        ],
        "interesting": "The C. jejuni genome is only ~30% G+C, which explains why its Cas9 evolved a "
                       "long pyrimidine-rich PAM instead of a G-rich one.",
        "presentation": "AT-rich curved rod, tiny CjCas9",
    },
    {
        "id": "saureus",
        "name": "Staphylococcus aureus NCTC 8325",
        "short": "S. aureus",
        "domain": "Bacteria",
        "phylum": "Bacillota (Firmicutes)",
        "gram": "positive",
        "morphology": "grape-like cocci, ~1 µm",
        "habitat": "human nares and skin; a major cause of drug-resistant infection",
        "genome_mb": 2.82, "gc": 0.329, "genes": 2900,
        "accession": "NC_007795.1",
        "systems": [
            {
                "name": "CRISPR (type II-A) — origin of SaCas9",
                "type": "II", "subtype": "II-A",
                "cas_genes": ["cas9", "cas1", "cas2", "csn2"],
                "effector": "Cas9 (SaCas9, 1053 aa)",
                "pam": "5'-NNGRRT-3'",
                "pam_iupac": "NNGRRT", "pam_side": "3",
                "protospacer_len": 21,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 3,
                "activity": "active; the 21-nt spacer and NNGRRT PAM limit off-target space",
                "notes": "SaCas9 packaging fits one AAV; used in the first in-body CRISPR therapy attempts.",
                "ref": "Ran 2015 Nature 520:186",
            }
        ],
        "interesting": "SaCas9 is ~300 amino acids smaller than SpCas9 — the difference that makes "
                       "single-vector AAV delivery possible.",
        "presentation": "gram-positive coccus, AAV-friendly SaCas9",
    },
    {
        "id": "tthemophilus",
        "name": "Thermus thermophilus HB8",
        "short": "T. thermophilus",
        "domain": "Bacteria",
        "phylum": "Deinococcota (Thermus)",
        "gram": "negative",
        "morphology": "rod, ~2 µm, thermophilic (65-72 °C)",
        "habitat": "hot springs — its proteins are stable enough for structural biology",
        "genome_mb": 1.85, "gc": 0.697, "genes": 2200,
        "accession": "NC_006461.1",
        "systems": [
            {
                "name": "CRISPR1 (type I-E)",
                "type": "I", "subtype": "I-E",
                "cas_genes": ["cas3", "cse1", "cse2", "cas7", "cas5", "cas6e", "cas1", "cas2"],
                "effector": "Tt-Cascade + Cas3",
                "pam": "5'-AAG-3'",
                "pam_iupac": "AAG", "pam_side": "5",
                "protospacer_len": 32,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 15,
                "activity": "the structural model system for type I — the first full Cascade structure",
                "notes": "Type IV/III-B systems are also present; thermostability made the first "
                         "atomic models of Cascade possible.",
                "ref": "Jore 2011 Nat Struct Mol Biol 18:529; Lintner 2011 J Biol Chem 286:21643",
            },
            {
                "name": "CRISPR2 (type III-B)",
                "type": "III", "subtype": "III-B",
                "cas_genes": ["cmr1", "cmr2", "cmr3", "cmr4", "cmr5", "cmr6", "cas1", "cas2", "cas6"],
                "effector": "Cmr complex — RNA-targeting, co-transcriptional",
                "pam": "none (RNA target)",
                "pam_iupac": "", "pam_side": "none",
                "protospacer_len": 36,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 10,
                "activity": "degrades invader RNA in vitro",
                "notes": "Type III-B Cmr complexes were the first CRISPR effectors shown to cut RNA.",
                "ref": "Hale 2009 Cell 139:945",
            },
        ],
        "interesting": "Hot-spring bacteria need heat-stable defence proteins; those same proteins "
                       "gave us the structural snapshots of CRISPR immunity.",
        "presentation": "thermophile with both type I-E and type III-B systems",
    },
    {
        "id": "haloferax",
        "name": "Haloferax volcanii DS2",
        "short": "H. volcanii",
        "domain": "Archaea",
        "phylum": "Euryarchaeota (Halobacteria)",
        "gram": "negative (archaeal, no peptidoglycan layer)",
        "morphology": "irregular disc/plate-shaped cell, ~2-3 µm",
        "habitat": "hypersaline lakes such as the Dead Sea",
        "genome_mb": 4.01, "gc": 0.653, "genes": 3800,
        "accession": "NC_013967",
        "systems": [
            {
                "name": "CRISPR locus (type I-B)",
                "type": "I", "subtype": "I-B",
                "cas_genes": ["cas3", "cas8b(cas1)", "cas7", "cas5", "cas6", "cas1", "cas2", "cas4"],
                "effector": "Cascade (I-B) + Cas3",
                "pam": "5'-TTC-3' (3' of the protospacer on the target strand)",
                "pam_iupac": "TTC", "pam_side": "5",
                "protospacer_len": 30,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 20,
                "activity": "one of the most experimentally tractable archaeal systems "
                            "— spacer acquisition and primed adaptation were both dissected here",
                "notes": "Archaea get half the CRISPR world's attention: the first in vivo spacer "
                         "acquisition assay was built in this organism.",
                "ref": "Barrangou 2007; Li 2014 Nucleic Acids Res 42:2483; Maier 2019 Nature 571:219",
            }
        ],
        "interesting": "The archaeal kingdom invented CRISPR too — every major CRISPR discovery "
                       "has been mirrored in Archaea, often first.",
        "presentation": "halophilic archaeon, salinity-tolerant type I-B defence",
    },
    {
        "id": "llactis",
        "name": "Lactococcus lactis IL1403",
        "short": "L. lactis",
        "domain": "Bacteria",
        "phylum": "Bacillota (Firmicutes)",
        "gram": "positive",
        "morphology": "ovoid cocci in pairs/chains, ~0.5-1 µm",
        "habitat": "milk — a workhorse of industrial fermentation",
        "genome_mb": 2.37, "gc": 0.354, "genes": 2300,
        "accession": "NC_002662.1",
        "systems": [
            {
                "name": "CRISPR1 (type II-A) and CRISPR2/3 (type III-A)",
                "type": "II + III", "subtype": "II-A / III-A",
                "cas_genes": ["cas9", "cas1", "cas2", "csn2", "cas10", "csm2-6"],
                "effector": "Cas9 for II-A; Csm complex for III-A",
                "pam": "5'-NGG-3' (II-A)",
                "pam_iupac": "NGG", "pam_side": "3",
                "protospacer_len": 20,
                "repeat": "", "repeat_len": 0, "repeat_source": "type-level",
                "spacers_observed": 30,
                "activity": "active; industrial strains are selected for phage-resistant arrays",
                "notes": "Dairy fermentations fail when phages win — CRISPR arrays are under constant "
                         "selection in starter cultures.",
                "ref": "Millen 2012; Barrangou & Horvath 2012",
            }
        ],
        "interesting": "Cheese-making is an arms race: phages can wipe out a starter culture in hours, "
                       "and CRISPR is the bacterium's answer.",
        "presentation": "dairy fermenter with layered II-A and III-A defences",
    },
]

BACTERIA_BY_ID: Dict[str, Dict[str, Any]] = {b["id"]: b for b in BACTERIA}

# --------------------------------------------------------------------------- #
# Phages and mobile elements
# --------------------------------------------------------------------------- #
PHAGES: List[Dict[str, Any]] = [
    {
        "id": "phix174",
        "name": "Escherichia phage phiX174",
        "host": "Escherichia coli C (and relatives)",
        "family": "Microviridae", "genome_type": "ssDNA (circular)",
        "genome_size": 5386, "gc": 0.445, "genes": 11,
        "accession": "NC_001422.1",
        "file": "phix174.fasta",
        "notes": "The first DNA genome ever sequenced (Sanger 1977). Lytically infects E. coli in "
                 "~20 minutes; its single-stranded genome is a classic substrate for CRISPR "
                 "spacer-acquisition assays.",
        "virulence": 0.85, "latency_min": 20, "burst_size": 150,
        "ref": "Sanger 1977 Nature 265:687",
    },
    {
        "id": "lambda",
        "name": "Enterobacteria phage lambda",
        "host": "Escherichia coli K-12",
        "family": "Siphoviridae (Caudoviricetes)", "genome_type": "dsDNA (linear, cos ends)",
        "genome_size": 48502, "gc": 0.498, "genes": 73,
        "accession": "NC_001416.1",
        "notes": "The classic temperate phage: integrates as a prophage, excises on induction. "
                 "A clean model for the lysogeny decision that CRISPR can interrupt.",
        "virulence": 0.55, "latency_min": 45, "burst_size": 100,
        "ref": "Lederberg 1950; Ptashne 2004",
    },
    {
        "id": "t7",
        "name": "Escherichia phage T7",
        "host": "Escherichia coli (with the appropriate receptor)",
        "family": "Podoviridae (Caudoviricetes)", "genome_type": "dsDNA (linear, direct repeats)",
        "genome_size": 39937, "gc": 0.482, "genes": 56,
        "accession": "NC_001604.1",
        "notes": "Fast lytic phage; its tightly packed genome leaves few PAMs, which changes how "
                 "a targeting guide must be chosen.",
        "virulence": 0.95, "latency_min": 25, "burst_size": 200,
        "ref": "Dunn 1983 J Mol Biol 166:477",
    },
    {
        "id": "st2972",
        "name": "Streptococcus phage 2972",
        "host": "Streptococcus thermophilus DGCC7710",
        "family": "Siphoviridae (cos-type)", "genome_type": "dsDNA",
        "genome_size": 34704, "gc": 0.418, "genes": 43,
        "accession": "NC_007019.1",
        "notes": "The phage used in the original 2007 yogurt experiment. Challenging S. thermophilus "
                 "with 2972 produced bacteriophage-insensitive mutants carrying new CRISPR spacers.",
        "virulence": 0.9, "latency_min": 60, "burst_size": 120,
        "ref": "Barrangou 2007 Science 315:1709",
    },
    {
        "id": "jbd30",
        "name": "Pseudomonas phage JBD30 (and relatives)",
        "host": "Pseudomonas aeruginosa PA14",
        "family": "Siphoviridae", "genome_type": "dsDNA",
        "genome_size": 45600, "gc": 0.575, "genes": 60,
        "accession": "JX275870",
        "notes": "A real phage with a real counter-defence: it carries anti-CRISPR genes (acrIF1, "
                 "acrIF2). A bacterial population that is fully immune can still be infected.",
        "virulence": 0.8, "latency_min": 50, "burst_size": 100, "anti_crispr": ["AcrIF1", "AcrIF2"],
        "ref": "Bondy-Denomy 2013 Nature 493:429",
    },
    {
        "id": "ms2",
        "name": "Escherichia virus MS2",
        "host": "Escherichia coli (F+ strains)",
        "family": "Leviviridae", "genome_type": "ssRNA (+) linear",
        "genome_size": 3569, "gc": 0.516, "genes": 4,
        "accession": "NC_001417.2",
        "notes": "An RNA phage — invisible to the DNA-targeting CRISPR types I, II and V. "
                 "Only type VI (Cas13) or type III-B (Cmr) can fight it.",
        "virulence": 0.6, "latency_min": 15, "burst_size": 500,
        "ref": "Fiers 1976 Nature 260:500",
    },
]

#: Anti-CRISPR proteins — real counter-defence genes carried by (pro)phages.
ANTI_CRISPR: List[Dict[str, str]] = [
    {"name": "AcrIF1", "target": "type I-F (Csy complex)", "mechanism": "Binds the Csy complex backbone and blocks target recruitment.",
     "origin": "P. aeruginosa phage JBD5", "ref": "Bondy-Denomy 2013 Nature 493:429"},
    {"name": "AcrIF2", "target": "type I-F (Csy complex)", "mechanism": "Mimics DNA and occupies the complex surface.",
     "origin": "P. aeruginosa phage JBD30", "ref": "Bondy-Denomy 2013 Nature 493:429"},
    {"name": "AcrIE1", "target": "type I-E (Cascade)", "mechanism": "Binds Cascade's Cas8e subunit; blocks Cas3 recruitment.",
     "origin": "P. aeruginosa prophage", "ref": "Pawluk 2014 mBio 5:e01010"},
    {"name": "AcrIIA2", "target": "type II-A (SpCas9)", "mechanism": "Binds the Cas9-sgRNA complex and prevents DNA binding.",
     "origin": "Listeria monocytogenes prophage", "ref": "Rauch 2017 Cell 168:150"},
    {"name": "AcrIIA4", "target": "type II-A (SpCas9)", "mechanism": "Mimics the PAM-interacting domain; the best-studied Cas9 inhibitor.",
     "origin": "Listeria monocytogenes prophage", "ref": "Rauch 2017 Cell 168:150"},
    {"name": "AcrIIC1", "target": "type II-C (NmeCas9, and some II-A)", "mechanism": "Broad-spectrum: binds the HNH domain and blocks cleavage while DNA binding continues.",
     "origin": "N. meningitidis prophage", "ref": "Harrington 2017 Cell 170:1224"},
    {"name": "AcrIIC3", "target": "type II-C (NmeCas9)", "mechanism": "Forces Cas9 to dimerise and prevents target binding.",
     "origin": "N. meningitidis prophage", "ref": "Harrington 2017 Cell 170:1224"},
    {"name": "AcrVA1", "target": "type V-A (Cas12a)", "mechanism": "A nuclease that destroys the crRNA — Cas12a is left without a guide.",
     "origin": "Moraxella bovis prophage", "ref": "Watters 2018 Science 362:236"},
    {"name": "AcrVIB1", "target": "type VI-B (Cas13b)", "mechanism": "Binds Cas13b and blocks its collateral RNase activity.",
     "origin": "bacterial prophage", "ref": "Meeske 2020 Nature 578:381"},
]
