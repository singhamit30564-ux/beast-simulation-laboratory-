"""Curated CRISPR effector / editor database.

Every entry records what the simulation needs to compute a cut, a base edit or a
prime edit.  Values are taken from the primary literature and from the review
tables cited in :mod:`beastlab.data.references`.  Where a value is a consensus or
an approximation it is flagged in ``notes`` (``size_aa`` is the canonical
isoform length and can differ by a few residues between annotations).

Field guide
-----------
pam            IUPAC PAM recognised by the effector
pam_side       "3" -> PAM sits 3' of the protospacer (Cas9-type)
               "5" -> PAM sits 5' of the protospacer (Cas12a-type)
               "none" -> no DNA PAM (Cas13 = RNA target, PFS instead)
spacer_len     guide-derived protospacer length used by the simulation
cut_geometry   "blunt" | "staggered"
nts_cut        cut on the non-target strand, counted in nt from the PAM-distal
               end of the protospacer when pam_side == "3" (from the PAM-proximal
               end when pam_side == "5")
ts_cut         same coordinate on the target strand (equal to nts_cut = blunt)
seed_len       PAM-proximal nt that dominate specificity
"""
from __future__ import annotations

from typing import Any, Dict, List

# --------------------------------------------------------------------------- #
# DNA-targeting nucleases (class 2, type II and type V)
# --------------------------------------------------------------------------- #
NUCLEASES: List[Dict[str, Any]] = [
    {
        "id": "SpCas9", "name": "SpCas9 (wild type)", "family": "Cas9",
        "species": "Streptococcus pyogenes SF370",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.62,
        "notes": "The canonical genome-editing nuclease: blunt DSB 3 bp upstream of the NGG PAM.",
        "ref": "Jinek 2012 Science 337:816",
    },
    {
        "id": "SaCas9", "name": "SaCas9", "family": "Cas9",
        "species": "Staphylococcus aureus",
        "system": "Type II-A", "class": 2, "size_aa": 1053,
        "pam": "NNGRRT", "pam_side": "3", "spacer_len": 21,
        "cut_geometry": "blunt", "nts_cut": 18, "ts_cut": 18,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.55,
        "notes": "Small enough for a single AAV vector; PAM NNGRRT (R = A/G).",
        "ref": "Ran 2015 Nature 520:186",
    },
    {
        "id": "Nme2Cas9", "name": "Nme2Cas9", "family": "Cas9",
        "species": "Neisseria meningitidis",
        "system": "Type II-C", "class": 2, "size_aa": 1082,
        "pam": "NNNNCC", "pam_side": "3", "spacer_len": 22,
        "cut_geometry": "blunt", "nts_cut": 19, "ts_cut": 19,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.48,
        "notes": "Compact high-fidelity ortholog; pyrimidine-rich N4CC PAM.",
        "ref": "Edraki 2019 Mol Cell 73:714",
    },
    {
        "id": "Nme1Cas9", "name": "Nme1Cas9", "family": "Cas9",
        "species": "Neisseria meningitidis",
        "system": "Type II-C", "class": 2, "size_aa": 1082,
        "pam": "NNNNGATT", "pam_side": "3", "spacer_len": 24,
        "cut_geometry": "blunt", "nts_cut": 21, "ts_cut": 21,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.45,
        "notes": "Long 8 nt PAM; used in AAV and for allele-specific targeting.",
        "ref": "Hou 2013 PNAS 110:15644",
    },
    {
        "id": "CjCas9", "name": "CjCas9", "family": "Cas9",
        "species": "Campylobacter jejuni",
        "system": "Type II-C", "class": 2, "size_aa": 984,
        "pam": "NNNNRYAC", "pam_side": "3", "spacer_len": 22,
        "cut_geometry": "blunt", "nts_cut": 19, "ts_cut": 19,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.40,
        "notes": "Smallest commonly used Cas9; strict long PAM keeps it very specific.",
        "ref": "Kim 2017 Nat Commun 8:14500",
    },
    {
        "id": "GeoCas9", "name": "GeoCas9", "family": "Cas9",
        "species": "Geobacillus stearothermophilus",
        "system": "Type II-C", "class": 2, "size_aa": 1087,
        "pam": "NNNNCRAA", "pam_side": "3", "spacer_len": 22,
        "cut_geometry": "blunt", "nts_cut": 19, "ts_cut": 19,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "thermostable", "default_efficiency": 0.42,
        "notes": "Thermostable; survives harsh delivery conditions.",
        "ref": "Harrington 2017 Nat Commun 8:1424",
    },
    {
        "id": "St1Cas9", "name": "St1Cas9 (Sth1)", "family": "Cas9",
        "species": "Streptococcus thermophilus CRISPR1",
        "system": "Type II-A", "class": 2, "size_aa": 1121,
        "pam": "NNAGAAW", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.50,
        "notes": "Natural CRISPR1 effector of the dairy bacterium S. thermophilus (PAM NNAGAAW).",
        "ref": "Deveau 2008 J Bacteriol 190:1390",
    },
    {
        "id": "St3Cas9", "name": "St3Cas9 (Sth3)", "family": "Cas9",
        "species": "Streptococcus thermophilus CRISPR3",
        "system": "Type II-A", "class": 2, "size_aa": 1388,
        "pam": "NGGNG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.52,
        "notes": "CRISPR3 effector; PAM NGGNG, closely related to SpCas9.",
        "ref": "Horvath 2008 J Bacteriol 190:1401",
    },
    {
        "id": "FnCas9", "name": "FnCas9", "family": "Cas9",
        "species": "Francisella novicida",
        "system": "Type II-B", "class": 2, "size_aa": 1629,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.35,
        "notes": "Large type II-B effector; low tolerance of mismatches.",
        "ref": "Zetsche 2015 Cell 163:759",
    },
    {
        "id": "TdCas9", "name": "TdCas9", "family": "Cas9",
        "species": "Treponema denticola",
        "system": "Type II-C", "class": 2, "size_aa": 1425,
        "pam": "NAAAAN", "pam_side": "3", "spacer_len": 22,
        "cut_geometry": "blunt", "nts_cut": 19, "ts_cut": 19,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.38,
        "notes": "A/T-rich PAM, useful where NGG sites are unavailable.",
        "ref": "Esvelt 2013 Nat Methods 10:1116",
    },
    {
        "id": "SpCas9-NG", "name": "SpCas9-NG", "family": "Cas9 (engineered)",
        "species": "S. pyogenes Cas9, PAM-interacting domain engineered",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.50,
        "notes": "Broadened NG PAM; reduced activity at non-NGG sites.",
        "ref": "Nishimasu 2018 Science 361:1259",
    },
    {
        "id": "xCas9-3.7", "name": "xCas9 3.7", "family": "Cas9 (engineered)",
        "species": "S. pyogenes Cas9, evolved",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.45,
        "notes": "Phage-assisted evolution product; NG / GAA / GAT PAMs, higher fidelity than SpCas9.",
        "ref": "Hu 2018 Nature 556:57",
    },
    {
        "id": "SpG", "name": "SpG", "family": "Cas9 (engineered)",
        "species": "S. pyogenes Cas9, structure-guided",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGN", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.55,
        "notes": "NGN PAM; prefers NGG > NGA > NGT.",
        "ref": "Walton 2020 Science 368:290",
    },
    {
        "id": "SpRY", "name": "SpRY (near-PAMless)", "family": "Cas9 (engineered)",
        "species": "S. pyogenes Cas9, structure-guided",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NRN", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.45,
        "notes": "Nearly PAMless (NRN high efficiency, NYN reduced). Widest targeting range of the Cas9 toolbox.",
        "ref": "Walton 2020 Science 368:290",
    },
    {
        "id": "SpRYc", "name": "SpRYc (chimeric)", "family": "Cas9 (engineered)",
        "species": "SpCas9-SpCas9-NG chimera",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NNG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.48,
        "notes": "SpRY plus NNG activity; high fidelity at NNG PAMs.",
        "ref": "Walton 2020 Science 368:290",
    },
    {
        "id": "eSpCas9-1.1", "name": "eSpCas9(1.1)", "family": "Cas9 (high fidelity)",
        "species": "S. pyogenes Cas9, K848A/K1003A/R1060A",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.55,
        "notes": "Weakened non-target-strand contacts -> fewer off-target cuts.",
        "ref": "Slaymaker 2016 Science 351:84",
    },
    {
        "id": "SpCas9-HF1", "name": "SpCas9-HF1", "family": "Cas9 (high fidelity)",
        "species": "S. pyogenes Cas9, N497A/R661A/Q695A/Q926A",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.50,
        "notes": "Disrupted DNA phosphate-backbone contacts; requires near-perfect guides.",
        "ref": "Kleinstiver 2016 Nature 529:490",
    },
    {
        "id": "HypaCas9", "name": "HypaCas9", "family": "Cas9 (high fidelity)",
        "species": "S. pyogenes Cas9, N692A/M694A/Q695A/H698A",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.50,
        "notes": "Increased proofreading; can show a lower on-target rate as well.",
        "ref": "Chen 2017 Nature 550:407",
    },
    {
        "id": "evoCas9", "name": "evoCas9", "family": "Cas9 (high fidelity)",
        "species": "S. pyogenes Cas9, evolved",
        "system": "Type II-A", "class": 2, "size_aa": 1368,
        "pam": "NGG", "pam_side": "3", "spacer_len": 20,
        "cut_geometry": "blunt", "nts_cut": 17, "ts_cut": 17,
        "seed_len": 12, "tracrrna": True, "collateral": False,
        "temperature": "37 °C", "default_efficiency": 0.48,
        "notes": "Bacterial selection for specificity; 79-fold better discrimination in some assays.",
        "ref": "Casini 2018 Nat Biotechnol 36:265",
    },
    {
        "id": "FnCas12a", "name": "FnCas12a (Cpf1)", "family": "Cas12a",
        "species": "Francisella novicida U112",
        "system": "Type V-A", "class": 2, "size_aa": 1300,
        "pam": "TTN", "pam_side": "5", "spacer_len": 23,
        "cut_geometry": "staggered", "nts_cut": 18, "ts_cut": 23,
        "seed_len": 6, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.42,
        "notes": "5' T-rich PAM, staggered cut with a 5 nt 5' overhang; also usable for diagnostics.",
        "ref": "Zetsche 2015 Cell 163:759",
    },
    {
        "id": "AsCas12a", "name": "AsCas12a (Cpf1)", "family": "Cas12a",
        "species": "Acidaminococcus sp. BV3L6",
        "system": "Type V-A", "class": 2, "size_aa": 1307,
        "pam": "TTTV", "pam_side": "5", "spacer_len": 23,
        "cut_geometry": "staggered", "nts_cut": 18, "ts_cut": 23,
        "seed_len": 6, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.40,
        "notes": "Single crRNA, staggered cut 18/23 nt from the PAM, AT-rich PAM preference.",
        "ref": "Zetsche 2015 Cell 163:759",
    },
    {
        "id": "LbCas12a", "name": "LbCas12a", "family": "Cas12a",
        "species": "Lachnospiraceae bacterium ND2006",
        "system": "Type V-A", "class": 2, "size_aa": 1228,
        "pam": "TTTV", "pam_side": "5", "spacer_len": 23,
        "cut_geometry": "staggered", "nts_cut": 18, "ts_cut": 23,
        "seed_len": 6, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.45,
        "notes": "Highest activity Cas12a homolog in human cells; 5' stagger leaves the PAM intact, so it re-cuts.",
        "ref": "Zetsche 2015 Cell 163:759",
    },
    {
        "id": "AaCas12b", "name": "AaCas12b (C2c1)", "family": "Cas12b",
        "species": "Alicyclobacillus acidiphilus",
        "system": "Type V-B", "class": 2, "size_aa": 1108,
        "pam": "VTTV", "pam_side": "5", "spacer_len": 20,
        "cut_geometry": "staggered", "nts_cut": 17, "ts_cut": 23,
        "seed_len": 8, "tracrrna": True, "collateral": True,
        "temperature": "thermostable", "default_efficiency": 0.40,
        "notes": "Thermostable type V-B effector with a 5' T-rich PAM and staggered cut.",
        "ref": "Shmakov 2015 Mol Cell 60:385",
    },
    {
        "id": "Un1Cas12f1", "name": "Un1Cas12f1 (Cas14a)", "family": "Cas12f",
        "species": "uncultured archaeon",
        "system": "Type V-F", "class": 2, "size_aa": 529,
        "pam": "TTTR", "pam_side": "5", "spacer_len": 20,
        "cut_geometry": "staggered", "nts_cut": 16, "ts_cut": 20,
        "seed_len": 8, "tracrrna": True, "collateral": True,
        "temperature": "37 °C (engineered variants)", "default_efficiency": 0.22,
        "notes": "~529 aa miniature nuclease originally described as Cas14; needs a modified sgRNA scaffold for mammalian editing.",
        "ref": "Harrington 2018 Science 362:839",
    },
    {
        "id": "Cas12i", "name": "Cas12i", "family": "Cas12i",
        "species": "unspecified metagenome isolate",
        "system": "Type V-I", "class": 2, "size_aa": 1063,
        "pam": "TTN", "pam_side": "5", "spacer_len": 20,
        "cut_geometry": "staggered", "nts_cut": 17, "ts_cut": 21,
        "seed_len": 8, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.30,
        "notes": "Compact type V-I effector with a 5' TTN PAM; processes its own crRNA array.",
        "ref": "Yan 2019 Science 363:88",
    },
    {
        "id": "CasPhi", "name": "Cas12j (CasΦ)", "family": "Cas12j",
        "species": "big phage-encoded Cas (metagenomic)",
        "system": "Type V-J", "class": 2, "size_aa": 757,
        "pam": "TBN", "pam_side": "5", "spacer_len": 21,
        "cut_geometry": "staggered", "nts_cut": 17, "ts_cut": 22,
        "seed_len": 8, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.25,
        "notes": "Phage-encoded miniature effector; extremely compact, B = C/G/T.",
        "ref": "Pausch 2020 Science 369:333",
    },
]

# --------------------------------------------------------------------------- #
# RNA-targeting effectors (class 2, type VI)
# --------------------------------------------------------------------------- #
RNA_EFFECTORS: List[Dict[str, Any]] = [
    {
        "id": "LwaCas13a", "name": "LwaCas13a (C2c2)", "family": "Cas13a",
        "species": "Leptotrichia wadei",
        "system": "Type VI-A", "class": 2, "size_aa": 1389,
        "pam": "PFS: non-G 3' of target", "pam_side": "none", "spacer_len": 28,
        "cut_geometry": "RNA", "nts_cut": 0, "ts_cut": 0,
        "seed_len": 8, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.55,
        "notes": "Programmable RNA knockdown; collateral RNase activity powers SHERLOCK detection.",
        "ref": "Abudayyeh 2016 Science 353:aaf5573",
    },
    {
        "id": "PspCas13b", "name": "PspCas13b", "family": "Cas13b",
        "species": "Prevotella sp. P5-125",
        "system": "Type VI-B", "class": 2, "size_aa": 1181,
        "pam": "PFS: D (A/G/T), NAN/NNA", "pam_side": "none", "spacer_len": 30,
        "cut_geometry": "RNA", "nts_cut": 0, "ts_cut": 0,
        "seed_len": 8, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.60,
        "notes": "Most efficient Cas13 for mammalian RNA knockdown.",
        "ref": "Cox 2017 Science 358:1019",
    },
    {
        "id": "RfxCas13d", "name": "RfxCas13d (CasRx)", "family": "Cas13d",
        "species": "Ruminococcus flavefaciens XPD3002",
        "system": "Type VI-D", "class": 2, "size_aa": 967,
        "pam": "PFS: none strict", "pam_side": "none", "spacer_len": 22,
        "cut_geometry": "RNA", "nts_cut": 0, "ts_cut": 0,
        "seed_len": 8, "tracrrna": False, "collateral": True,
        "temperature": "37 °C", "default_efficiency": 0.65,
        "notes": "Compact, no PFS constraint, extremely efficient knockdown; isoform-specific targeting.",
        "ref": "Konermann 2018 Cell 173:665",
    },
]

# --------------------------------------------------------------------------- #
# Base editors
# --------------------------------------------------------------------------- #
BASE_EDITORS: List[Dict[str, Any]] = [
    {
        "id": "BE3", "name": "BE3 (rAPOBEC1-nCas9-UGI)", "family": "CBE",
        "conversion": "C>T (C:G -> T:A)", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (4, 8), "preferred_motif": "TC",
        "motif_order": ["TC", "CC", "AC", "GC"], "bystander": "moderate",
        "indel_rate": 0.05, "max_efficiency": 0.45, "purity": 0.85,
        "notes": "First widely used cytosine base editor; 3rd-gen architecture with UGI.",
        "ref": "Komor 2016 Nature 533:420",
    },
    {
        "id": "BE4max", "name": "BE4max", "family": "CBE",
        "conversion": "C>T (C:G -> T:A)", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (4, 8), "preferred_motif": "TC",
        "motif_order": ["TC", "CC", "AC", "GC"], "bystander": "moderate",
        "indel_rate": 0.04, "max_efficiency": 0.55, "purity": 0.90,
        "notes": "Codon-optimised, bipartite NLS, two UGI domains; the workhorse CBE.",
        "ref": "Koblan 2018 Nat Biotechnol 36:843",
    },
    {
        "id": "AncBE4max", "name": "AncBE4max", "family": "CBE",
        "conversion": "C>T", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (4, 8), "preferred_motif": "TC",
        "motif_order": ["TC", "CC", "AC", "GC"], "bystander": "moderate",
        "indel_rate": 0.04, "max_efficiency": 0.58, "purity": 0.92,
        "notes": "Ancestral-codon reconstruction of BE4max; better expression in some cell types.",
        "ref": "Koblan 2018 Nat Biotechnol 36:843",
    },
    {
        "id": "YE1-BE4", "name": "YE1-BE4", "family": "CBE (narrow)",
        "conversion": "C>T", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (5, 6), "preferred_motif": "TC",
        "motif_order": ["TC", "CC", "AC", "GC"], "bystander": "low",
        "indel_rate": 0.02, "max_efficiency": 0.25, "purity": 0.95,
        "notes": "Narrowed window (positions 5-6) to avoid bystander cytidines.",
        "ref": "Kim 2017 Nat Biotechnol 35:475",
    },
    {
        "id": "evoCDA1-BE4max", "name": "evoCDA1-BE4max", "family": "CBE (wide)",
        "conversion": "C>T", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (1, 13), "preferred_motif": "TC/GC",
        "motif_order": ["TC", "GC", "CC", "AC"], "bystander": "high",
        "indel_rate": 0.07, "max_efficiency": 0.50, "purity": 0.80,
        "notes": "Enlarged window (positions 1-13); can edit GC-context cytosines outside the classic window.",
        "ref": "Thuronyi 2019 Nat Biotechnol 37:1070",
    },
    {
        "id": "ABE7.10", "name": "ABE7.10 (TadA-TadA*-nCas9)", "family": "ABE",
        "conversion": "A>G (A:T -> G:C)", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (5, 7), "preferred_motif": "no strong preference",
        "motif_order": [], "bystander": "low",
        "indel_rate": 0.01, "max_efficiency": 0.40, "purity": 0.97,
        "notes": "First efficient adenine base editor; narrow window, minimal indels.",
        "ref": "Gaudelli 2017 Nature 551:464",
    },
    {
        "id": "ABEmax", "name": "ABEmax", "family": "ABE",
        "conversion": "A>G", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (4, 7), "preferred_motif": "TA/AA disfavoured",
        "motif_order": [], "bystander": "low",
        "indel_rate": 0.01, "max_efficiency": 0.50, "purity": 0.97,
        "notes": "Codon-optimised ABE7.10 with improved NLS.",
        "ref": "Koblan 2018 Nat Biotechnol 36:843",
    },
    {
        "id": "ABE8e", "name": "ABE8e (TadA-8e)", "family": "ABE",
        "conversion": "A>G", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (4, 8), "preferred_motif": "wide, disfavours A at -1",
        "motif_order": [], "bystander": "high",
        "indel_rate": 0.02, "max_efficiency": 0.70, "purity": 0.93,
        "notes": "590-fold faster deamination than ABE7.10; fast and processive, wider window -> more bystander edits.",
        "ref": "Richter 2020 Nat Biotechnol 38:883",
    },
    {
        "id": "ABE9", "name": "ABE9 (N108Q/L145T)", "family": "ABE (narrow)",
        "conversion": "A>G", "cas": "SpCas9 D10A nickase",
        "pam": "NGG", "window": (5, 6), "preferred_motif": "TCN context",
        "motif_order": [], "bystander": "very low",
        "indel_rate": 0.01, "max_efficiency": 0.40, "purity": 0.99,
        "notes": "Engineered narrow window of 1-2 nt, suppresses bystander adenines.",
        "ref": "Chen 2023 Nat Biotechnol 41:1295",
    },
    {
        "id": "ABE8e-SpRY", "name": "ABE8e-SpRY", "family": "ABE (PAM-flexible)",
        "conversion": "A>G", "cas": "SpRY D10A nickase",
        "pam": "NRN / NYN", "window": (4, 8), "preferred_motif": "wide",
        "motif_order": [], "bystander": "high",
        "indel_rate": 0.03, "max_efficiency": 0.60, "purity": 0.90,
        "notes": "Near-PAMless adenine editor; broaden scope but raise the off-target search space.",
        "ref": "Richter 2020 Nat Biotechnol 38:883; Walton 2020 Science 368:290",
    },
]

# --------------------------------------------------------------------------- #
# Prime editors
# --------------------------------------------------------------------------- #
PRIME_EDITORS: List[Dict[str, Any]] = [
    {
        "id": "PE1", "name": "PE1", "family": "PE",
        "cas": "SpCas9 H840A nickase + M-MLV RT",
        "pam": "NGG", "max_efficiency": 0.05, "indel_rate": 0.01,
        "notes": "First-generation prime editor; RT not optimised, low efficiency.",
        "ref": "Anzalone 2019 Nature 576:149",
    },
    {
        "id": "PE2", "name": "PE2", "family": "PE",
        "cas": "SpCas9 H840A + engineered M-MLV RT (5 mutations)",
        "pam": "NGG", "max_efficiency": 0.20, "indel_rate": 0.02,
        "notes": "Enhanced RT improves editing; single nick, no second guide.",
        "ref": "Anzalone 2019 Nature 576:149",
    },
    {
        "id": "PE3", "name": "PE3", "family": "PE",
        "cas": "PE2 + second nickase guide",
        "pam": "NGG", "max_efficiency": 0.35, "indel_rate": 0.08,
        "notes": "Extra nicking guide biases resolution to the edited strand; higher efficiency and higher indel risk.",
        "ref": "Anzalone 2019 Nature 576:149",
    },
    {
        "id": "PE3b", "name": "PE3b", "family": "PE",
        "cas": "PE2 + edited-strand-specific nickase guide",
        "pam": "NGG", "max_efficiency": 0.28, "indel_rate": 0.03,
        "notes": "Nick guide matches only the edited allele -> matched indel/edited ratio.",
        "ref": "Anzalone 2019 Nature 576:149",
    },
    {
        "id": "PEmax", "name": "PEmax (PE2max architecture)", "family": "PE",
        "cas": "PE2 with R221K/N394K RT + optimised NLS",
        "pam": "NGG", "max_efficiency": 0.45, "indel_rate": 0.05,
        "notes": "Improved architecture/nuclear import; current therapeutic-grade PE variant.",
        "ref": "Chen 2021 Cell 184:5635",
    },
]


ALL_DNA_EFFECTORS: List[Dict[str, Any]] = NUCLEASES

EFFECTOR_BY_ID: Dict[str, Dict[str, Any]] = {
    e["id"]: e for e in (NUCLEASES + RNA_EFFECTORS + BASE_EDITORS + PRIME_EDITORS)
}

#: Table of the PAM-flexible SpCas9 family, used by the PAM-coverage explorer.
PAM_COVERAGE_SETS: Dict[str, str] = {
    "SpCas9 (NGG)": "NGG",
    "SpG (NGN)": "NGN",
    "SpCas9-NG (NG)": "NG",
    "xCas9 3.7 (NG)": "NG",
    "SpRY (NRN)": "NRN",
    "SpRY (NYN)": "NYN",
    "SpRYc (NNG)": "NNG",
}

#: Verification timeline — mutations reverted by base or prime editing.
POINT_MUTATION_CLASSES: List[Dict[str, str]] = [
    {"change": "A>G / T>C", "editor": "ABE (ABE8e, ABE9)", "comment": "Purine swap; no DSB, no need for a donor template."},
    {"change": "C>T / G>A", "editor": "CBE (BE4max, YE1)", "comment": "Transition; the most common pathogenic SNV class."},
    {"change": "C>G / G>C", "editor": "CBE + UGI tuning (or PE)", "comment": "Transversion — low-purity CBE products."},
    {"change": "A>T / T>A", "editor": "Prime editor", "comment": "Not accessible to base editors."},
    {"change": "A>C / T>G", "editor": "Prime editor", "comment": "Not accessible to base editors."},
    {"change": "G>T / C>A", "editor": "Prime editor", "comment": "Not accessible to base editors."},
    {"change": "Insertion / deletion", "editor": "PE2/PE3, HDR", "comment": "PE is the only DSB-free route for small indels."},
]
