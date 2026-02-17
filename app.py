import streamlit as st
import requests
import io
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBList
from stmol import showmol
import py3Dmol
import os
import tempfile
from fpdf import FPDF # 追加

# --- データベース設定 ---
LIPID_DB = {
    "DMPC (C14:0)": 23.0,
    "DPPC (C16:0)": 28.5,
    "POPC (Standard)": 30.0,
    "DOPC (Unsaturated)": 29.0
}

KD_SCALE = {
    'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
    'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
    'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
    'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
}

# --- 関数群 ---
def predict_structure(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        res = requests.post(url, data=sequence, timeout=180) 
        return res.text if res.status_code == 200 else None
    except:
        return None

def fetch_pdb_by_id(pdb_id):
    try:
        pdbl = PDBList()
        file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        with open(file_path, 'r') as f:
            pdb_data = f.read()
        os.remove(file_path)
        return pdb_data
    except:
        return None

def analyze_protein(pdb_string):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", io.StringIO(pdb_string))
    
    coords, hydro = [], []
    for model in structure:
        for residue in model.get_residues():
            res_name = residue.get_resname()
            if res_name in KD_SCALE and 'CA' in residue:
                coords.append(residue['CA'].get_coord())
                hydro.append(KD_SCALE[res_name])
    
    coords = np.array(coords)
    if len(coords) == 0: return 0, 0, ([], [], 0, 0), coords

    rg = np.sqrt(np.mean(np.sum((coords - np.mean(coords, axis=0))**2, axis=1)))
    
    z_coords = coords[:, 2]
    sorted_idx = np.argsort(z_coords)
    window = 10
    if len(z_coords) > window:
        z_smooth = np.convolve(z_coords[sorted_idx], np.ones(window)/window, mode='valid')
        h_smooth = np.convolve(np.array(hydro)[sorted_idx], np.ones(window)/window, mode='valid')
    else:
        z_smooth, h_smooth = z_coords, np.array(hydro)
    
    belt_z = z_smooth[h_smooth > 1.0]
    b_min, b_max = (np.min(belt_z), np.max(belt_z)) if len(belt_z) > 0 else (0, 0)
    
    return rg, b_max - b_min, (z_smooth, h_smooth, b_min, b_max), coords

# --- PDF生成関数 (新規追加) ---
def create_pdf(rg, belt_width, score, msp_name, est_diam, fig_plot, lipid_name):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", 'B', 16)
    pdf.cell(0, 10, "GPCR-Nanodisc Analysis Report", ln=True, align='C')
    pdf.ln(10)

    pdf.set_font("Arial", size=12)
    pdf.cell(0, 10, f"Target Lipid: {lipid_name}", ln=True)
    pdf.ln(5)

    # Metrics Table
    pdf.set_font("Arial", 'B', 12)
    pdf.cell(60, 10, "Metric", 1)
    pdf.cell(60, 10, "Value", 1)
    pdf.ln()
    
    pdf.set_font("Arial", size=12)
    metrics = [
        ("Radius of Gyration (Rg)", f"{rg:.2f} A"),
        ("Hydrophobic Belt Width", f"{belt_width:.2f} A"),
        ("Compatibility Score", f"{score:.2f}"),
        ("Estimated Diameter", f"{est_diam:.2f} nm")
    ]
    for metric, value in metrics:
        pdf.cell(60, 10, metric, 1)
        pdf.cell(60, 10, value, 1)
        pdf.ln()
    
    pdf.ln(10)
    pdf.set_font("Arial", 'B', 14)
    pdf.cell(0, 10, f"Recommendation: {msp_name}", ln=True, align='L')
    pdf.ln(5)

    # Save plot to temp file and add to PDF
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmpfile:
        fig_plot.savefig(tmpfile.name, dpi=150)
        pdf.image(tmpfile.name, x=10, y=None, w=190)
    
    # Return PDF bytes
    return pdf.output(dest='S').encode('latin-1')

# --- Streamlit UI ---
st.set_page_config(page_title="GPCR-Nanodisc Predictor", layout="wide")
st.title("🧬 GPCR-Nanodisc Integration Predictor")

st.sidebar.header("1. Target Selection")
input_method = st.sidebar.radio("Input Method", ["Paste Sequence (ESMFold)", "Upload PDB File", "Enter PDB ID"])

pdb_data = None
pdb_name = "Custom Protein" # PDF用

if input_method == "Paste Sequence (ESMFold)":
    seq_input = st.sidebar.text_area("Amino Acid Sequence", height=150)
    if st.sidebar.button("Predict Structure"):
        if seq_input:
            with st.spinner("Predicting..."):
                pdb_data = predict_structure(seq_input.strip())
                if not pdb_data: st.error("Failed.")
elif input_method == "Upload PDB File":
    uploaded_file = st.sidebar.file_uploader("Upload .pdb file", type="pdb")
    if uploaded_file:
        pdb_data = uploaded_file.getvalue().decode("utf-8")
        pdb_name = uploaded_file.name
elif input_method == "Enter PDB ID":
    pdb_id_input = st.sidebar.text_input("PDB ID (e.g., 2RH1)", max_chars=4)
    if st.sidebar.button("Fetch PDB"):
        with st.spinner(f"Fetching {pdb_id_input}..."):
            pdb_data = fetch_pdb_by_id(pdb_id_input)
            pdb_name = f"PDB ID: {pdb_id_input}"
            if not pdb_data: st.error("Failed.")

st.sidebar.header("2. Lipid Settings")
lipid_choice = st.sidebar.selectbox("Select Target Lipid", list(LIPID_DB.keys()))
target_thickness = LIPID_DB[lipid_choice]

if pdb_data:
    rg, b_width, plot_data, coords = analyze_protein(pdb_data)
    
    if rg == 0:
        st.error("Could not parse structure.")
    else:
        score = max(0, 1.0 - abs(b_width - target_thickness) / target_thickness)

        st.subheader("Analysis Results")
        c1, c2, c3 = st.columns(3)
        c1.metric("Rg", f"{rg:.2f} Å")
        c2.metric("Belt Width", f"{b_width:.2f} Å")
        c3.metric("Score", f"{score:.2f}")

        est_diam = rg * 2.5 / 10
        msp_name = "MSP1D1" if est_diam < 9.7 else "MSP1E3D1" if est_diam < 12.1 else "MSP2N2"
        st.info(f"Recommended Scaffold: **{msp_name}**")

        # 3D View
        st.subheader("🧊 3D Visualization")
        view = py3Dmol.view(width=700, height=400)
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        center_z = (plot_data[2] + plot_data[3]) / 2
        view.addShape({'type': 'Box', 'center': {'x': 0, 'y': 0, 'z': center_z},
                       'dimensions': {'w': 60, 'h': 60, 'd': b_width}, 'color': 'yellow', 'opacity': 0.3})
        view.zoomTo()
        showmol(view, height=400, width=800)

        # Plot
        st.subheader("📊 Hydrophobic Profile")
        fig, ax = plt.subplots(figsize=(10, 4)) # グラフサイズ調整
        ax.plot(plot_data[0], plot_data[1], color='blue', label='Hydrophobicity')
        ax.axvspan(plot_data[2], plot_data[3], color='yellow', alpha=0.3, label='Belt Area')
        ax.axhline(y=1.0, linestyle='--', color='gray', alpha=0.5)
        ax.set_title(f"Profile vs {lipid_choice}")
        ax.set_xlabel("Z-coordinate (Å)")
        ax.legend()
        st.pyplot(fig)

        # --- PDFダウンロードボタン (新規追加) ---
        st.divider()
        st.subheader("📥 Export Report")
        
        # PDF生成
        pdf_bytes = create_pdf(rg, b_width, score, msp_name, est_diam, fig, lipid_choice)
        
        st.download_button(
            label="Download Analysis Report (PDF)",
            data=pdf_bytes,
            file_name=f"analysis_report_{msp_name}.pdf",
            mime="application/pdf"
        )