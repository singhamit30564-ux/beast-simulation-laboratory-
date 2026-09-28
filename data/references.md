# Data provenance and references

Everything the simulator shows is traceable to one of these sources. Sequence
data are shipped in `data/genomes/` and are real: NCBI GenBank/RefSeq records
pulled with E-utilities, and GRCh38 windows from the Ensembl REST API.

## Sequence data

| Record | Coordinates / accession | Source |
| --- | --- | --- |
| HBB locus (incl. sickle codon) | chr11:5,225,464-5,227,071 (GRCh38) | Ensembl REST `sequence/region/human/11:5225464..5227071:1` |
| BCL11A erythroid enhancer | chr2:60,495,000-60,495,500 (GRCh38) | Ensembl REST `sequence/region/human/2:60495000..60495500:1` |
| HBG1 promoter (HPFH region) | chr11:5,249,390-5,249,890 (GRCh38) | Ensembl REST `sequence/region/human/11:5249390..5249890:1` |
| PCSK9 5' region | chr1:55,039,450-55,039,950 (GRCh38) | Ensembl REST `sequence/region/human/1:55039450..55039950:1` |
| TTR exon 3 region | chr18:31,594,000-31,594,500 (GRCh38) | Ensembl REST `sequence/region/human/18:31594000..31594500:1` |
| CEP290 IVS26 region | chr12:88,086,000-88,086,500 (GRCh38) | Ensembl REST `sequence/region/human/12:88086000..88086500:1` |
| CCR5 locus region | chr3:46,370,550-46,371,050 (GRCh38) | Ensembl REST `sequence/region/human/3:46370550..46371050:1` |
| AAVS1 / PPP1R12C region | chr19:55,115,000-55,115,500 (GRCh38) | Ensembl REST `sequence/region/human/19:55115000..55115500:1` |
| *E. coli* K-12 CRISPR-1 array (leader-proximal) | U00096.3:2,877,800-2,878,849 | NCBI E-utilities `efetch` |
| *E. coli* K-12 lac-region window | U00096.3:365,030-365,529 | NCBI E-utilities `efetch` |
| phiX174 complete genome | NC_001422.1 (5,386 bp) | NCBI E-utilities `efetch` |

The *E. coli* CRISPR array (repeat `GGTTTATCCCCGCTGGCGCGGGGAACTC`, 9 spacers of
33-34 nt) is **parsed out of the chromosome window at run time** — the app reads
the array rather than hard-coding it.

## CRISPR-Cas biology

- Jinek M. et al. *A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity.* Science 337:816 (2012). — SpCas9, blunt cut 3 bp upstream of NGG.
- Deltcheva E. et al. *CRISPR RNA maturation by trans-encoded small RNA and host factor RNase III.* Nature 471:602 (2011). — tracrRNA.
- Gasunas G. et al. / Nishimasu H. et al. — PAM recognition by the Cas9 PI domain (R1333/R1335 read NGG in the major groove).
- Zetsche B. et al. *Cpf1 is a single RNA-guided endonuclease of a class 2 CRISPR-Cas system.* Cell 163:759 (2015). — Cas12a: TTTV PAM, staggered cut 18/23 nt, 5 nt 5' overhang.
- Shmakov S. et al. *Discovery and functional characterization of diverse class 2 CRISPR-Cas systems.* Mol Cell 60:385 (2015). — Cas12b/C2c1.
- Harrington E. et al. *Programmed DNA destruction by miniature CRISPR-Cas14 enzymes.* Science 362:839 (2018). — Cas12f/Cas14.
- Yan W. et al. Science 363:88 (2019) — Cas12i; Pausch P. et al. Science 369:333 (2020) — CasΦ/Cas12j.
- Abudayyeh O. et al. Science 353:aaf5573 (2016) — Cas13a; Cox D. et al. Science 358:1019 (2017) — Cas13b; Konermann S. et al. Cell 173:665 (2018) — Cas13d/CasRx.
- Kleinstiver B. et al. Nature 529:490 (2016) — SpCas9-HF1; Slaymaker I. et al. Science 351:84 (2016) — eSpCas9(1.1); Chen J. et al. Nature 550:407 (2017) — HypaCas9; Casini A. et al. Nat Biotechnol 36:265 (2018) — evoCas9.
- Hu J. et al. Nature 556:57 (2018) — xCas9; Nishimasu H. et al. Science 361:1259 (2018) — SpCas9-NG; Walton R. et al. Science 368:290 (2020) — SpG and SpRY.
- Ran F. et al. Nature 520:186 (2015) — SaCas9; Kim E. et al. Nat Commun 8:14500 (2017) — CjCas9; Edraki A. et al. Mol Cell 73:714 (2019) — Nme2Cas9; Harrington L. et al. Nat Commun 8:1424 (2017) — GeoCas9.
- Makarova K. et al. *Evolutionary classification of CRISPR-Cas systems: a burst of class 2 and derived variants.* Nat Rev Microbiol 18:67 (2020). — system classification used in the effector table.
- Hsu P. et al. Nat Biotechnol 31:827 (2013) — MIT/Hsu position weights used for the off-target score.

## Bacterial immunity and anti-phage defence

- Barrangou R. et al. *CRISPR provides acquired resistance against viruses in prokaryotes.* Science 315:1709 (2007). — the yogurt experiment; S. thermophilus CRISPR1/CRISPR3.
- Horvath P. et al. J Bacteriol 190:1401 (2008). — repeat sequences and spacer content of S. thermophilus CRISPR1/2/3.
- Deveau H. et al. J Bacteriol 190:1390 (2008). — PAMs NNAGAAW (CRISPR1) and NGGNG (CRISPR3).
- Brouns S. et al. Science 321:960 (2008) — Cascade and the 61 nt crRNA; Jore M. et al. Nat Struct Mol Biol 18:529 (2011) — Cascade structure.
- Datsenko K. et al. PNAS 109:E2579 (2012) — *E. coli* type I-E adaptation; primed acquisition.
- Brouns 2008 / Pougach 2010 — polarised spacer insertion at the leader end.
- Cady K. et al. J Bacteriol 193:3433 (2011) — P. aeruginosa type I-F interference; Bondy-Denomy J. et al. Nature 493:429 (2013) — anti-CRISPR discovery.
- Rauch B. et al. Cell 168:150 (2017) — AcrIIA2/AcrIIA4; Harrington L. et al. Cell 170:1224 (2017) — AcrIIC1/AcrIIC3; Watters K. et al. Science 362:236 (2018) — AcrVA1; Meeske A. et al. Nature 578:381 (2020) — AcrVIB1.
- Sanger F. et al. Nature 265:687 (1977) — phiX174, the first sequenced genome.

## Repair, base editing and prime editing

- Komor A. et al. Nature 533:420 (2016) — BE3; Koblan L. et al. Nat Biotechnol 36:843 (2018) — BE4max/ABEmax.
- Gaudelli N. et al. Nature 551:464 (2017) — ABE7.10; Richter M. et al. Nat Biotechnol 38:883 (2020) — ABE8e (~590-fold faster deamination).
- Chen L. et al. Nat Biotechnol 41:1295 (2023) — ABE9 narrow window.
- Thuronyi B. et al. Nat Biotechnol 37:1070 (2019) — evoCDA1 wide-window CBE.
- Anzalone A. et al. *Search-and-replace genome editing without double-strand breaks or donor DNA.* Nature 576:149 (2019). — prime editing, PBS/RT-template design, PE2/PE3.
- Chen P. et al. Cell 184:5635 (2021) — PEmax.
- Ceccaldi R. et al. Trends Cell Biol 26:52 (2016) — repair pathway choice at DSBs; Sfeir & Symington 2015 — MMEJ.

## Clinical programmes

- Frangoul H. et al. *CRISPR-Cas9 gene editing for sickle cell disease and beta-thalassemia.* NEJM 384:252 (2021) — CTX001/Casgevy.
- FDA labels: CASGEVY STN 125787 (8 Dec 2023, SCD) and STN 125788 (16 Jan 2024, TDT).
- Gillmore J. et al. *CRISPR-Cas9 in vivo gene editing for transthyretin amyloidosis.* NEJM 385:493 (2021) — NTLA-2001.
- Kathiresan S. et al. Circulation 147:1002 (2022) — VERVE-101 non-human primate data; heart-1 (NCT05398029) human data 2023.
- Maeder M. et al. Nat Med 25:229 (2019) — CEP290 IVS26 targeting; BRILLIANCE (NCT03872479).
- Lu Y. et al. Nature 579:262 (2020); Stadtmauer E. et al. Science 367:eaba7365 (2020) — multiplex-edited T cells.
- Xu L. et al. NEJM 381:1240 (2019) — CCR5-edited T cells in an HIV-1 patient.
- Platt O. et al. NEJM 330:1639 (1994) — HbF and sickle-cell severity.

## What this simulator does *not* claim

- The on-target score is a transparent heuristic built from GC content, base
  preferences, homopolymers and specificity — **not** Doench/Azimuth Rule Set 2
  or the CFD matrix, which are trained on measured data and are not reproduced here.
- Indel spectra are drawn from an empirical distribution whose shape matches
  published Cas9 outcomes; individual alleles are not predicted from repair
  biophysics, except for microhomology-driven deletions, which *are* computed
  from the actual flanking sequence.
- Population dynamics use textbook phage-bacteria parameters. They reproduce the
  arms race qualitatively (collapse, rebound, escape) rather than any specific
  published time-course.
- Clinical cohort projections are illustrations of dose-response logic, not
  trial predictions.
