# 🧬 GPCR-Nanodisc Integration Predictor

**GPCRをナノディスクに再構成する際の構造適合性を定量評価するバイオインフォマティクスWebアプリケーション**

A bioinformatics web application that quantitatively evaluates the structural compatibility of G-Protein Coupled Receptors (GPCRs) with lipid nanodiscs — from sequence input to scaffold recommendation.

---

## 🎯 Motivation

Reconstituting GPCRs into lipid nanodiscs is essential for structural and functional studies of membrane proteins in a near-native environment. However, successful reconstitution depends on two critical parameters that researchers typically optimize through trial-and-error:

1. **Scaffold sizing** — The GPCR must physically fit within the MSP (Membrane Scaffold Protein) belt.
2. **Hydrophobic matching** — The lipid acyl-chain length must align with the protein's transmembrane hydrophobic surface.

Mismatch in either dimension causes protein aggregation, instability, or loss of activity, wasting expensive reagents and time.

**This tool replaces empirical guesswork with computational pre-screening**, allowing researchers to evaluate GPCR–nanodisc compatibility *before* entering the wet lab.

---

## 📐 Analysis Pipeline

```
┌─────────────────────┐
│  Input               │
│  ├─ Amino acid seq   │──→ ESMFold API ──→ 3D Structure (PDB)
│  ├─ PDB file upload  │──→──────────────┐
│  └─ PDB ID fetch     │──→ RCSB PDB ──→─┤
└─────────────────────┘                   ▼
                                ┌──────────────────┐
                                │  Structure Parse  │
                                │  (Biopython)      │
                                │  Cα coordinates   │
                                │  + residue types  │
                                └────────┬─────────┘
                         ┌───────────────┼───────────────┐
                         ▼               ▼               ▼
                ┌──────────────┐ ┌──────────────┐ ┌──────────────┐
                │  Rg Calc     │ │  KD Profile  │ │  3D Render   │
                │  → Est. Dia. │ │  → Belt Det. │ │  + Belt Box  │
                │  → MSP Rec.  │ │  → Score     │ │  (py3Dmol)   │
                └──────────────┘ └──────────────┘ └──────────────┘
                         │               │               │
                         ▼               ▼               ▼
                ┌────────────────────────────────────────────────┐
                │  Output: Metrics / 3D View / Plot / PDF Report │
                └────────────────────────────────────────────────┘
```

---

## 🌟 Key Features

| Feature | Description |
|---|---|
| **Multi-Input Support** | Paste a sequence (→ ESMFold), upload a `.pdb` file, or fetch by PDB ID from RCSB |
| **Radius of Gyration (Rg)** | Computed from Cα coordinates; used to estimate protein diameter and select MSP size |
| **Hydrophobic Belt Detection** | Z-axis Kyte-Doolittle profile with sliding-window smoothing; auto-detects transmembrane belt boundaries |
| **Compatibility Score** | Quantifies belt width vs. lipid bilayer thickness mismatch (0–1 scale) |
| **MSP Recommendation** | Auto-selects MSP1D1 / MSP1E3D1 / MSP2N2 based on estimated disc diameter |
| **Interactive 3D Viewer** | Cartoon representation with semi-transparent membrane boundary overlay (py3Dmol) |
| **PDF Export** | Downloadable report containing all metrics, interpretation, and hydrophobicity plot |

---

## 🧪 Scientific Basis

### Radius of Gyration → Scaffold Selection

$$R_g = \sqrt{ \frac{1}{N} \sum_{i=1}^{N} (\mathbf{r}_i - \bar{\mathbf{r}})^2 }$$

The estimated disc diameter is derived as $d \approx R_g \times 2.5 / 10$ (nm), then matched to the smallest compatible MSP:

| MSP | Max Disc Diameter |
|---|---|
| MSP1D1 | ≤ 9.7 nm |
| MSP1E3D1 | ≤ 12.1 nm |
| MSP2N2 | ≤ 16.5 nm |

> Ref: Denisov, I.G. & Sligar, S.G. *Chem. Rev.* **117**, 4669–4713 (2017)

### Hydrophobic Belt → Lipid Compatibility

The Kyte-Doolittle hydrophobicity scale is applied per residue along the Z-axis. Residues with KD score > 1.0 are classified as hydrophobic; the Z-span of the hydrophobic cluster defines the belt width.

Compatibility score:

$$S = \max\!\left(0,\ 1 - \frac{|W_{\text{belt}} - T_{\text{lipid}}|}{T_{\text{lipid}}}\right)$$

where $W_{\text{belt}}$ is the detected belt width (Å) and $T_{\text{lipid}}$ is the target bilayer thickness (Å).

| Lipid | Hydrophobic Thickness (Å) |
|---|---|
| DMPC (C14:0) | 23.0 |
| DPPC (C16:0) | 28.5 |
| DOPC (C18:1) | 29.0 |
| POPC (C16:0/C18:1) | 30.0 |

> Ref: Kučerka, N. et al. *Biochim. Biophys. Acta* **1261**, 1512–1520 (2011)

---

## 🚀 Getting Started

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

Then open `http://localhost:8501` in your browser.

**Quick test:** Enter PDB ID `2RH1` (β2 adrenergic receptor) → Select POPC → View results.

---

## 🗂️ Project Structure

```
gpcr-nanodisc-predictor/
├── app.py               # Main application (Streamlit + analysis logic)
├── requirements.txt     # Python dependencies
├── README.md
└── .gitignore
```

---

## 🛠️ Tech Stack

| Category | Tools |
|---|---|
| Frontend | Streamlit |
| Structure Prediction | ESMFold API |
| Bioinformatics | Biopython (PDB parsing) |
| Visualization | py3Dmol / stmol, Matplotlib |
| Reporting | fpdf2 |
| Numerical Analysis | NumPy |

---

## ⚠️ Limitations & Future Work

- **Orientation dependency:** The hydrophobic belt detection assumes the transmembrane axis is roughly aligned with the Z-axis. Pre-aligned structures (e.g., from OPM database) yield more accurate results.
- **Single-chain analysis:** Multi-chain complexes (e.g., GPCR–G protein) are parsed as a single entity; per-chain analysis is not yet supported.
- **Lipid diversity:** Currently limited to four symmetric phospholipids. Future versions could incorporate asymmetric bilayer compositions and cholesterol effects.
- **Validation:** Systematic benchmarking against experimentally validated GPCR–nanodisc pairs is planned.

---

## 🧑‍🔬 Author

**TSUBAKI0531**

Ph.D. in Agriculture — Protein Engineering & Antibody Drug Discovery

This tool was developed as part of a self-directed effort to bridge wet-lab membrane protein research with computational pre-screening, combining domain expertise in GPCR reconstitution with bioinformatics implementation.

---

## 📄 License

MIT
