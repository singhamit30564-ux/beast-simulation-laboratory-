# beast-simulation-laboratory-
<p align="center">
  <img src="https://img.shields.io/badge/CRISPR-SIMULATOR-0a0e17?style=for-the-badge&logo=apachespark&logoColor=d4af37" />
  <img src="https://img.shields.io/badge/Python-3.10+-blue?style=for-the-badge&logo=python" />
</p>

<h1 align="center">🧬 Beast Simulation Laboratory</h1>
<h3 align="center">Digital Genome Surgery Simulator</h3>

---

## Overview

Beast Simulation Laboratory is a Python-based toolkit for simulating CRISPR-based genome editing experiments. It provides in-silico models for Cas9, Cas13, Base Editors, and Prime Editors with off-target analysis and repair pathway prediction.

## Features

### CRISPR-Cas9 Simulator
- gRNA design and on-target scoring
- PAM scanning (SpCas9, SpRY, xCas9, Sniper-Cas9)
- Double-strand break simulation
- HDR vs NHEJ repair pathway prediction
- Indel outcome profiling

### CRISPR-Cas13 Simulator
- ssRNA targeting and binding simulation
- Collateral cleavage modeling (SHERLOCK/DETECTR)
- Transcriptome off-target scanning
- Cas13a/b/d/x variant support

### Base Editor Simulator
- Cytosine Base Editor (CBE): C→T conversion
- Adenine Base Editor (ABE): A→G conversion
- Byproduct and indel analysis
- Editing window prediction

### Prime Editor Simulator
- PegRNA design with PBS and RT template
- All 12 types of point mutations
- Insertion and deletion simulation
- PE2/PE3/PEmax variant support

### Off-Target Analysis
- Whole genome in-silico scanning
- CFD (Cutting Frequency Determination) scoring
- Aggregate off-target risk assessment
- Genome-wide heatmap visualization

## Repository Structure

beast-simulation-laboratory/ ├── src/ │ ├── crispr/ │ │ ├── cas9_sim.py │ │ ├── cas13_sim.py │ │ ├── base_editor.py │ │ └── prime_editor.py │ ├── genome/ │ │ ├── pam_scanner.py │ │ ├── off_target.py │ │ ├── repair_pathway.py │ │ └── variant_library.py │ ├── visual/ │ │ ├── dna_render.py │ │ ├── cut_site.py │ │ └── repair_anim.py │ └── utils/ │ ├── sequence.py │ ├── scoring.py │ └── formatter.py ├── experiments/ │ ├── sickle_cell_cure.py │ ├── beta_thalassemia.py │ ├── muscular_dystrophy.py │ ├── malaria_resistance.py │ ├── gfp_reporter.py │ └── cancer_immunotherapy.py ├── notebooks/ │ ├── crispr_workflow.ipynb │ ├── base_editing_guide.ipynb │ └── prime_editing_advanced.ipynb ├── data/ │ ├── pam_library.json │ ├── repair_templates.json │ └── scoring_matrices/ ├── requirements.txt ├── LICENSE └── README.md

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/beast-simulation-laboratory.git
cd beast-simulation-laboratory
pip install -r requirements.txt

## Quick Start
from src.crispr.cas9_sim import Cas9Simulator
from src.genome.pam_scanner import PAMScanner

cas9 = Cas9Simulator(variant="SpCas9")
target_seq = "ATCGATCGATCGATCGNNGRRT"
pam = PAMScanner.find(target_seq, variant="SpCas9")

result = cas9.cut(
    target=target_seq,
    pam=pam,
    repair_pathway="NHEJ",
    cell_type="HEK293"
)

print(f"Cut Efficiency: {result.efficiency}%")
print(f"Indel Profile: {result.indels}")
print(f"Off-targets: {len(result.off_targets)}")

## Supported CRISPR Variants
Variant	PAM	Use Case	
SpCas9	NGG	Standard editing	
SpCas9-NG	NG	Relaxed PAM	
xCas9	NG, GAA, GAT	Broad PAM	
SpRY	NRN/NYN	Near PAM-less	
Sniper-Cas9	NGG	High fidelity	
Cas13a	Protospacer flanking site	RNA editing	
Cas13b	DR sequence	RNA editing	
Cas13d	Protospacer	Compact RNA editing	

## Output Formats
JSON
 
CSV
 
VCF (Variant Call Format)
 
HTML Report
 
MP4/GIF Animation

## Contributions 
Fork the repository
2. 
Create a feature branch
3. 
Commit your changes
4. 
Push to the branch
5. 
Open a Pull Request

## License 
MIT License 

## Acknowledgments
 Jennifer Doudna and Emmanuelle Charpentier for CRISPR-Cas9
 
David Liu for Base Editing and Prime Editing