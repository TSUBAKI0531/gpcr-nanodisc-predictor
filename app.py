import streamlit as st
import requests
import io
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBList
from stmol import showmol
import py3Dmol
import os

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
        # タイムアウトを少し延長
        res = requests.post(url, data=sequence, timeout=180) 
        return res.text if res.status_code == 200 else None
    except:
        return None

def fetch_pdb_by_id(pdb_id):
    try:
        pdbl = PDBList()
        # 一時ディレクトリにダウンロード
        file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')
        with open(file_path, 'r') as f:
            pdb_data = f.read()
        os.remove(file_path) # クリーンアップ
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
    # スムージング処理
    window = 10
    if len(z_coords) > window:
        z_smooth = np.convolve(z_coords[sorted_idx], np.ones(window)/window, mode='valid')
        h_smooth = np.convolve(np.array(hydro)[sorted_idx], np.ones(window)/window, mode='valid')
    else:
        z_smooth, h_smooth = z_coords, np.array(hydro)
    
    belt_z = z_smooth[h_smooth > 1.0]
    b_min, b_max = (np.min(belt_z), np.max(belt_z)) if len(belt_z) > 0 else (0, 0)
    
    return rg, b_max - b_min, (z_smooth, h_smooth, b_min, b_max), coords

# --- Streamlit UI ---
st.set_page_config(page_title="GPCR-Nanodisc Predictor", layout="wide")
st.title("🧬 GPCR-Nanodisc Integration Predictor")

# サイドバー設定
st.sidebar.header("1. Target Selection")
input_method = st.sidebar.radio("Input Method", ["Paste Sequence (ESMFold)", "Upload PDB File", "Enter PDB ID"])

pdb_data = None

# 入力メソッドごとの処理
if input_method == "Paste Sequence (ESMFold)":
    seq_input = st.sidebar.text_area("Amino Acid Sequence", height=150)
    if st.sidebar.button("Predict Structure"):
        if seq_input:
            with st.spinner("Predicting structure via ESMFold API (may take time)..."):
                pdb_data = predict_structure(seq_input.strip()) # 空白除去
                if not pdb_data:
                    st.error("ESMFold API failed. Please try 'Enter PDB ID' with '2RH1' instead.")

elif input_method == "Upload PDB File":
    uploaded_file = st.sidebar.file_uploader("Upload .pdb file", type="pdb")
    if uploaded_file:
        pdb_data = uploaded_file.getvalue().decode("utf-8")

elif input_method == "Enter PDB ID":
    pdb_id_input = st.sidebar.text_input("PDB ID (e.g., 2RH1)", max_chars=4)
    if st.sidebar.button("Fetch PDB"):
        with st.spinner(f"Fetching {pdb_id_input}..."):
            pdb_data = fetch_pdb_by_id(pdb_id_input)
            if not pdb_data:
                st.error("Failed to fetch PDB ID.")

st.sidebar.header("2. Lipid Settings")
lipid_choice = st.sidebar.selectbox("Select Target Lipid", list(LIPID_DB.keys()))
target_thickness = LIPID_DB[lipid_choice]

# 解析実行と表示
if pdb_data:
    rg, b_width, plot_data, coords = analyze_protein(pdb_data)
    
    if rg == 0:
        st.error("Could not parse structure data.")
    else:
        score = max(0, 1.0 - abs(b_width - target_thickness) / target_thickness)

        # メトリクス
        st.subheader("Analysis Results")
        c1, c2, c3 = st.columns(3)
        c1.metric("Radius of Gyration ($R_g$)", f"{rg:.2f} Å")
        c2.metric("Hydrophobic Belt", f"{b_width:.2f} Å")
        c3.metric("Compatibility Score", f"{score:.2f}")

        # 推奨MSP
        est_diam = rg * 2.5 / 10
        msp_name = "MSP1D1" if est_diam < 9.7 else "MSP1E3D1" if est_diam < 12.1 else "MSP2N2"
        st.info(f"Recommended Scaffold: **{msp_name}** (Estimated Diameter: {est_diam:.2f} nm)")

        # 3D可視化
        st.subheader("🧊 3D Visualization")
        view = py3Dmol.view(width=700, height=400)
        view.addModel(pdb_data, 'pdb')
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        # 膜領域ボックス
        center_z = (plot_data[2] + plot_data[3]) / 2
        view.addShape({'type': 'Box', 'center': {'x': 0, 'y': 0, 'z': center_z},
                       'dimensions': {'w': 60, 'h': 60, 'd': b_width}, 'color': 'yellow', 'opacity': 0.3})
        view.zoomTo()
        showmol(view, height=400, width=800)

        # グラフ
        st.subheader("📊 Hydrophobic Profile")
        fig, ax = plt.subplots(figsize=(10, 3))
        ax.plot(plot_data[0], plot_data[1], color='blue', label='Hydrophobicity')
        ax.axvspan(plot_data[2], plot_data[3], color='yellow', alpha=0.3, label='Belt Area')
        ax.axhline(y=1.0, linestyle='--', color='gray', alpha=0.5)
        ax.set_xlabel("Z-coordinate (Å)")
        ax.set_ylabel("Score")
        ax.legend()
        st.pyplot(fig)

elif input_method == "Paste Sequence (ESMFold)" and not pdb_data:
    st.info("👈 Enter a sequence and click 'Predict Structure' to start.")