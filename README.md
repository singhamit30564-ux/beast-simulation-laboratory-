# 🧬 BEAST Simulation Laboratory

An interactive **CRISPR laboratory** built on Streamlit: real sequences, real effector
parameters, real repair biology — wrapped in simulations you can turn the knobs on.

There is no "run" button that returns a prediction. Every bench is a model with the assumptions
printed next to it, so the interesting question ("what happens if I change this?") has an answer
you can see change.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

---

## The benches

| Bench | What you do |
| --- | --- |
| 🛰️ **Mission control** | Orientation: what the app is, which question each bench answers, and the vocabulary you need. |
| 🦠 **Immunity arena** | Pick a bacterium and an invading phage. Read its CRISPR array (E. coli K-12's array is parsed live out of `data/genomes/ecoli_k12.fasta`), let Cas1–Cas2 write new spacers, then run the phage × bacteria arms race — and mutate the target to see an escape mutant get out of jail. |
| 🔬 **Cas explorer** | The whole effector toolbox: nucleases (**SpCas9, SaCas9, Nme1/2Cas9, CjCas9, GeoCas9, St1/3Cas9, FnCas9, TdCas9, SpCas9-NG, xCas9-3.7, SpG, SpRY, SpRYc, eSpCas9-1.1, SpCas9-HF1, HypaCas9, evoCas9, Fn/As/LbCas12a, AaCas12b, Un1Cas12f1, Cas12i, CasPhi**), RNA-targeting **Cas13a/b/d**, **10 base editors**, **5 prime editors** — with PAM coverage measured on real sequence. |
| 🎯 **Guide design** | Scan real IUPAC PAMs, rank guides, look at the off-target table, build the guide RNA, and design an HDR donor or a pegRNA. The on-target score is a transparent heuristic and the UI says so — it is *not* Doench/Azimuth Rule Set 2. |
| ✂️ **Genome surgery** | The full experiment: target, effector, cell type, delivery, dose, cell-cycle context → allele fates, indel spectrum, sequence-level product, simulated gel, collateral-damage radar. Four modalities: nuclease, HDR, base editing, prime editing. |
| 🩹 **Repair lab** | Where the outcome is actually decided: c-NHEJ vs MMEJ vs HDR vs SSA, the microhomologies that really exist around your cut, how the indel spectrum changes between cell types, and which HDR levers move the number. |
| 💊 **Therapeutics** | Five indications, seven clinical programmes (Casgevy, NTLA-2001, VERVE-101, EDIT-101, …), efficiency → endpoint curves, and a patient-cohort projection with between-patient variability. |
| 🗄️ **Database** | The reference layer: shipped genomes with accessions, twelve bacteria and their CRISPR systems, phages and anti-CRISPR proteins, the searchable effector catalogue, and the provenance file. |

## The simulation core (`beastlab/`)

```
beastlab/
  sequtils.py       dependency-free sequence toolkit (IUPAC, revcomp, Tm, ORFs, microhomology…)
  crispr.py         PAM scanning, guide scoring, off-target search, donor & pegRNA design
  repair.py         DSB repair pathway competition, indel spectra, base/prime editing, safety
  immunity.py       CRISPR arrays, adaptation, interference, phage × bacteria dynamics
  theme.py          the visual system (palette, cards, sequence renderers)
  viz.py            Plotly figures and SVG biology
  data/             cas_enzymes.py · microbes.py · targets.py · clinical.py
  pages/            one module per bench (pure Streamlit, no logic hidden in the UI)
```

Everything is pure Python + numpy/pandas/plotly — no Biopython, no compiled dependencies, so the
models are readable in place.

## Where the data comes from

* **Human loci** — real GRCh38 windows fetched from the Ensembl REST API (HBB sickle allele,
  BCL11A +58 enhancer, HBG1 HPFH promoter, PCSK9, TTR, CEP290 IVS26, AAVS1, CCR5).
* **Bacterial and phage genomes** — NCBI RefSeq: *E. coli* K-12 MG1655 U00096.3 (the type I-E
  CRISPR-1 array plus a coding-dense window) and phiX174 NC_001422.1, fetched with E-utilities.
* **Repeats, PAMs, editors, anti-CRISPRs, trials** — from the primary literature, FDA labels and
  trial press releases; every entry carries its reference. The full list, with the exact fetch
  commands, is in [`data/references.md`](data/references.md).

Type-level sequences (a system's canonical repeat rather than the strain's own array) are labelled
`type-level` in the data and in the UI. Nothing here pretends to be a strain it is not.

## What this simulator does **not** claim

* The on-target score is a heuristic, not a trained model; it is only useful for *ranking* guides.
* Indel sizes are drawn from an empirical spectrum, not from a repair biophysics simulation.
* Off-target search is exact matching with a mismatch budget inside the supplied sequence, not a
  whole-genome Cas-OFFinder/Bowtie search.
* The clinical layer maps published biomarker relationships — it is a projection, not a prognosis.

See the "what this simulator does not claim" section of `data/references.md`, and read the
assumption lists that accompany every model output.

## License

MIT — see [LICENSE](LICENSE).
