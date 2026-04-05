# 🧬 GPCR-Nanodisc Integration Predictor

A bioinformatics web application that evaluates the structural compatibility of G-Protein Coupled Receptors (GPCRs) with lipid nanodiscs.

Combines **ESMFold** structure prediction, **Radius of Gyration (Rg)** analysis for scaffold sizing, and **Kyte-Doolittle hydrophobic belt** profiling to help researchers optimise nanodisc reconstitution conditions.

---

## Key Features

| Feature | Description |
|---|---|
| **Multi-Input Support** | Paste a sequence (ESMFold), upload a `.pdb` file, or fetch by PDB ID |
| **Scaffold Recommendation** | Auto-selects MSP1D1 / MSP1E3D1 / MSP2N2 based on estimated diameter |
| **Hydrophobic Mismatch Score** | Quantifies belt–lipid compatibility (0–1 scale) for DMPC, DPPC, POPC, DOPC |
| **Interactive 3D Viewer** | Cartoon representation with semi-transparent membrane belt overlay |
| **PDF Export** | Downloadable report with all metrics and the hydrophobicity plot |

---

## Scientific Background

Reconstituting GPCRs into nanodiscs requires precise matching of:

1. **Physical dimensions** — the protein must fit within the MSP scaffold belt.
2. **Hydrophobic match** — the lipid acyl-chain length should align with the GPCR's transmembrane hydrophobic surface.

This tool quantifies both aspects and recommends an appropriate scaffold protein.

---

## Getting Started

### Prerequisites

- Python 3.9+

### Installation

```bash
git clone https://github.com/TSUBAKI0531/gpcr-nanodisc-predictor.git
cd gpcr-nanodisc-predictor
pip install -r requirements.txt
```

### Usage

```bash
streamlit run app.py
```

---

## Tech Stack

- **Frontend:** Streamlit
- **Bioinformatics:** Biopython, ESMFold API
- **Visualisation:** Py3Dmol / Stmol, Matplotlib
- **Reporting:** FPDF
- **Analysis:** NumPy

---

## Author

**TSUBAKI0531**

Ph.D. in Agriculture · Specialising in Protein Engineering & Antibody Drug Discovery
