import streamlit as st
import requests
import io
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from stmol import showmol
import py3Dmol

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

# --- ロジック関数 ---
def predict_structure(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        res = requests.post(url, data=sequence, timeout=120)
        return res.text if res.status_code == 200 else None
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
    # 1. 回転半径 (Rg) 計算
    rg = np.sqrt(np.mean(np.sum((coords - np.mean(coords, axis=0))**2, axis=1)))
    
    # 2. 疎水性ベルト解析 (Z軸方向)
    z_coords = coords[:, 2]
    sorted_idx = np.argsort(z_coords)
    z_smooth = np.convolve(z_coords[sorted_idx], np.ones(10)/10, mode='valid')
    h_smooth = np.convolve(np.array(hydro)[sorted_idx], np.ones(10)/10, mode='valid')
    
    belt_z = z_smooth[h_smooth > 1.0] # 閾値1.0以上を疎水性領域とみなす
    b_min, b_max = (np.min(belt_z), np.max(belt_z)) if len(belt_z) > 0 else (0, 0)
    
    return rg, b_max - b_min, (z_smooth, h_smooth, b_min, b_max), coords

# --- Streamlit UI ---
st.set_page_config(page_title="GPCR-Nanodisc Predictor", layout="wide")
st.title("🧬 GPCR-Nanodisc Integration Predictor")

st.sidebar.header("Analysis Settings")
lipid_choice = st.sidebar.selectbox("Select Target Lipid", list(LIPID_DB.keys()))
target_thickness = LIPID_DB[lipid_choice]
seq_input = st.sidebar.text_area("Amino Acid Sequence", height=200)

if st.sidebar.button("Run Prediction"):
    if seq_input:
        with st.spinner("Processing..."):
            pdb_data = predict_structure(seq_input)
            if pdb_data:
                rg, b_width, plot_data, coords = analyze_protein(pdb_data)
                score = max(0, 1.0 - abs(b_width - target_thickness) / target_thickness)

                # メトリクス表示
                c1, c2, c3 = st.columns(3)
                c1.metric("Radius of Gyration ($R_g$)", f"{rg:.2f} Å")
                c2.metric("Hydrophobic Belt", f"{b_width:.2f} Å")
                c3.metric("Compatibility Score", f"{score:.2f}")

                # MSP推奨
                est_diam = rg * 2.5 / 10 # 近似計算
                msp_name = "MSP1D1" if est_diam < 9.7 else "MSP1E3D1" if est_diam < 12.1 else "MSP2N2"
                st.info(f"Recommended Scaffold: **{msp_name}** (Estimated Diameter: {est_diam:.2f} nm)")

                # 3D可視化
                st.subheader("🧊 3D Structure & Membrane Belt")
                view = py3Dmol.view(width=700, height=400)
                view.addModel(pdb_data, 'pdb')
                view.setStyle({'cartoon': {'color': 'spectrum'}})
                view.addShape({'type': 'Box', 'center': {'x': 0, 'y': 0, 'z': (plot_data[2]+plot_data[3])/2},
                               'dimensions': {'w': 60, 'h': 60, 'd': b_width}, 'color': 'yellow', 'opacity': 0.3})
                view.zoomTo()
                showmol(view, height=400, width=800)

                # プロファイル表示
                st.subheader("📊 Hydrophobic Profile")
                fig, ax = plt.subplots(figsize=(10, 3))
                ax.plot(plot_data[0], plot_data[1], color='blue', label='Hydrophobicity')
                ax.axvspan(plot_data[2], plot_data[3], color='yellow', alpha=0.3, label='Belt Area')
                ax.set_xlabel("Z-coordinate (Å)")
                ax.legend()
                st.pyplot(fig)
            else:
                st.error("Prediction failed. Try a shorter sequence.")