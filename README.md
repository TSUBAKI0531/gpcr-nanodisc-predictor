# 🧬 GPCR-Nanodisc Integration Predictor

**GPCR-Nanodisc Integration Predictor** is a bioinformatics web application designed to evaluate the stability of G-Protein Coupled Receptors (GPCRs) when reconstituted into lipid nanodiscs.

By integrating **ESMFold** for structure prediction, **Radius of Gyration ($R_g$)** analysis for scaffold sizing, and **Hydrophobic Belt** analysis for lipid compatibility, this tool aids researchers in optimizing experimental conditions for membrane protein studies.

---

## 🌟 Key Features

### 1. Multi-Input Support
* **PDB ID Fetch:** Directly download structures from the RCSB PDB (e.g., input `2RH1`).
* **File Upload:** Support for custom `.pdb` files (e.g., AlphaFold predictions).
* **AI Prediction:** Generate 3D structures from amino acid sequences using the **ESMFold API**.

### 2. Smart Scaffold Recommendation
* Automatically calculates the protein's **Radius of Gyration ($R_g$)**.
* Suggests the optimal Membrane Scaffold Protein (MSP) size:
    * **MSP1D1** (~9.7 nm)
    * **MSP1E3D1** (~12.1 nm)
    * **MSP2N2** (~16.5 nm)

### 3. Hydrophobic Mismatch Analysis
* Visualizes the transmembrane region using the **Kyte-Doolittle** hydrophobicity scale.
* Scores the compatibility between the protein's hydrophobic belt and the target lipid bilayer (e.g., DMPC, POPC).

### 4. Interactive Visualization & Reporting
* **3D Viewer:** Real-time visualization of the protein and predicted membrane boundary using `py3Dmol`.
* **PDF Export:** Generate and download a comprehensive analysis report containing metrics and hydrophobicity plots for lab notebooks.

---

## 🧪 Scientific Background

Reconstituting GPCRs into nanodiscs requires precise matching of the lipid bilayer thickness and the scaffold diameter to the target protein. Mismatches can lead to protein aggregation or instability.

This tool quantifies:
1.  **Physical Dimensions:** Ensuring the GPCR fits within the MSP belt.
2.  **Hydrophobic Match:** Ensuring the lipid acyl chain length matches the GPCR's hydrophobic surface.

---

## 🚀 Getting Started

### Prerequisites
* Python 3.9+

### Installation
```bash
git clone [https://github.com/your-username/gpcr-nanodisc-predictor.git](https://github.com/your-username/gpcr-nanodisc-predictor.git)
cd gpcr-nanodisc-predictor
pip install -r requirements.txt
Usage
Bash
streamlit run app.py
🛠️ Tech Stack
Frontend: Streamlit

Bioinformatics: Biopython, ESMFold API

Visualization: Py3Dmol, Stmol, Matplotlib

Reporting: FPDF

Analysis: NumPy

🧑‍🔬 Author
[Your Name / ユーザー名]

Ph.D. in Agriculture

Specializing in Protein Engineering & Antibody Drug Discovery