# 🧬 BEAST Simulation Laboratory

An interactive **CRISPR laboratory** built on Streamlit: real sequences, real effector
parameters, real repair biology — wrapped in simulations you can turn the knobs on.

Every bench is a model with the assumptions
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

The core uses Python + numpy/pandas/plotly; optional button-gated educational ML/RL uses CPU torch. PDF export uses ReportLab. All model source is readable in place.

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

* The original on-target ranking score is a heuristic. The separate MLP is trained only on simulated labels, not experimental data.
* Indel sizes are drawn from an empirical spectrum, not from a repair biophysics simulation.
* Off-target search is exact matching with a mismatch budget inside the supplied sequence, not a
  whole-genome Cas-OFFinder/Bowtie search.
* The clinical layer maps published biomarker relationships — it is a projection, not a prognosis.

See the "what this simulator does not claim" section of `data/references.md`, and read the
assumption lists that accompany every model output.

## License

MIT — see [LICENSE](LICENSE).

## Priority Feature Pack v2

All eight benches retain their navigation. New tools live here:

- **Therapeutics → Sickle cell case study**: HBB disease context → a labeled synthetic
  BCL11A enhancer teaching fragment (or your fragment/guide) → exact SpCas9/PAM check →
  chosen NHEJ deletion → toy HbF curve → classroom verdict. HBB is **not** the enhancer
  target. The demo is not a clinical guide. “Functional cure predicted” is explicitly
  an educational threshold rule, never evidence of a cure.
- **Genome surgery**: six 2D Plotly/SVG molecular cartoon frames, Play/Pause/Replay.
  Lite mode replaces animation with six numbered static frame panels.
- **Guide design**: strict single-record FASTA/raw DNA (10,000 bp maximum), RNA/DNA
  20-base custom spacer; user sequences are never transmitted to Ensembl. Invalid
  symbols are rejected, not silently removed. Results invalidate when inputs change.
- **Export center**: PDF, CSV and JSON for custom-tool results, ranked guides and the
  case study. The SHA-256 covers the canonical JSON `payload` (sorted keys, compact
  separators, ASCII encoding), **not** the surrounding file. CSV stores that complete
  payload in one JSON field; PDF prints it. Re-hash that payload to verify. This is
  integrity checking, not a digital signature or proof of biological correctness.
- **Scan vs GRCh38**: fetches one selected **bounded HBB or BCL11A reference window**
  from Ensembl REST. This is **not a whole-genome scan**. Searches both strands with
  immediately adjacent NGG, ≤3 substitutions, including exact matches; reports forward
  1-based coordinates and the PAM-proximal 12-nt seed. No bulges, chromatin, variants,
  or measured cleavage risks. Offline failures visibly switch to bundled phiX174,
  which says **NOT human**; a dedicated offline button also works without internet.
  phiX174 is scanned linearly (circular-origin spanning sites are not evaluated).
  Public API responses are RAM-cached for 24 hours, maximum four entries, timeout 8 s.
- **Dr. Titan**: three field notes per bench, SHA-256/seed-42 stable initial choice,
  and a session-stable Next button. Unicode scientist avatar needs no image download.
- **Sidebar Lite mode**: session-persisted toggle, 50-row display cap, supplementary
  Plotly charts omitted. Full calculated results remain in export payloads. Clear
  session inputs & results removes user data and optional in-memory retraining.

### Educational ML / RL

**Simulated training data — illustrative model, not clinical-grade.**

`beastlab/ml/efficiency.py` supplies a CPU-only torch MLP (24→32→16→1). The 24
features are GC fraction, 16 overlapping dinucleotide fractions, four positional
weighted quarter summaries, and a three-category PAM one-hot. Labels come from an
explicit GC/position/poly-T/PAM formula plus Gaussian noise (SD 0.09), clipped into
[0.01, 0.99]; **none are experimental measurements**. Seed 42, 5,000 guides, 4,000
train / 1,000 held out. The displayed “confidence band” is the held-out 90th percentile
absolute residual, clipped to [0,100]; it has **no biological coverage guarantee**.

The Guide design training panel shows loss history or optionally retrains in RAM.
Torch and weights load only on model-button clicks; shared shipped inference models
use `st.cache_resource`, are eval-only and CPU-only with one compute thread. Session
retraining does not mutate the shared resource or the RL reward model.

`beastlab/rl/guide_env.py` has Gymnasium-style `reset(seed, options)` and five-value
`step(action)` (no Gym dependency). Observation: 20×4 one-hot; actions 0–79 =
position×4 + base. Reward: change in MLP efficiency percentage points. Stop at 10
steps or >90% efficiency. DQN: 80→64→80, replay 1,000, target network, epsilon decay,
20 greedy imitation trajectories / 60 warm-start epochs, 120 training episodes.
Every fifth episode starts at a heuristic guide. UI reports the actual DQN/random/
greedy result on the same starting guide, seed 42 and step budget. Greedy uses more
score evaluations; this is **not an equal-compute or population benchmark**. Mutations
may destroy target binding. This is a score-optimization sandbox, not usable guide design.

Small shipped state dictionaries and metadata live in `models/efficiency_mlp/` and
`models/guide_dqn/`. Regenerate (in this order) with:

```bash
python -m beastlab.ml.efficiency
python -m beastlab.rl.dqn
```

No user input is used in training or written to checkpoints. Requirements select a
CPU torch wheel from the public PyTorch CPU index; no GPU, paid API, keys, sklearn,
or RL framework is needed. Model architecture and bounded training use little RAM;
Streamlit/torch themselves still have nonzero baseline overhead.

### Validation

```bash
pip install -r requirements-dev.txt
python -m py_compile app.py $(find beastlab -name '*.py')
python -m pytest -q
```

Tests cover all eight benches in normal/Lite mode using Streamlit AppTest, model
shapes/determinism/range, environment reset/actions/reward, exact PAM adjacency and
reverse coordinates, fallback, case verdicts, animation structure, input validation,
export hashes, stale results and model/export UI actions. Tests do not require a live API.

**Zero Data Retention — computed in RAM** means no application-level user-input disk
storage or guide submission to third parties. Inputs/results remain in the active
Streamlit session until cleared or expired; public reference caches are shared. This
is not a claim about hosting-provider/network logs. Downloaded exports are under the
user's control. Do not paste identifiable patient information.
