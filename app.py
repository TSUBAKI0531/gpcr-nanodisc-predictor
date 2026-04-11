"""
GPCR-Nanodisc Integration Predictor
====================================
Evaluates GPCR stability in lipid nanodiscs by combining ESMFold structure
prediction, Radius-of-Gyration analysis, and hydrophobic-belt profiling.

Author: TSUBAKI0531
"""

from __future__ import annotations

import io
import os
import tempfile
from dataclasses import dataclass
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import requests
import streamlit as st
from Bio.PDB import PDBParser
from fpdf import FPDF
import streamlit.components.v1 as components

# =============================================================================
# Constants & Configuration
# =============================================================================

LIPID_DB: dict[str, float] = {
    "DMPC (C14:0)": 23.0,
    "DPPC (C16:0)": 28.5,
    "POPC (Standard)": 30.0,
    "DOPC (Unsaturated)": 29.0,
}

# Kyte-Doolittle hydrophobicity scale
KD_SCALE: dict[str, float] = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

# MSP scaffold database: (max diameter nm, full name)
MSP_CATALOG: list[tuple[float, str]] = [
    (9.7, "MSP1D1"),
    (12.1, "MSP1E3D1"),
    (16.5, "MSP2N2"),
]

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ESMFOLD_TIMEOUT_SEC = 180
SMOOTHING_WINDOW_DEFAULT = 10
HYDROPHOBIC_THRESHOLD = 1.0  # KD score above which residue is "hydrophobic"
RG_TO_DIAMETER_FACTOR = 2.5  # empirical factor: diameter ~ Rg * 2.5 / 10


# =============================================================================
# Data classes
# =============================================================================

@dataclass
class ProfileData:
    """Smoothed Z-coordinate and hydrophobicity arrays plus belt boundaries."""
    z_smooth: np.ndarray
    h_smooth: np.ndarray
    belt_min: float
    belt_max: float


@dataclass
class AnalysisResult:
    """All outputs of a single protein analysis."""
    rg: float
    belt_width: float
    profile: ProfileData
    coords: np.ndarray
    valid: bool = True


# =============================================================================
# Structure acquisition
# =============================================================================

def predict_structure(sequence: str) -> Optional[str]:
    """Predict 3D structure from amino-acid sequence via ESMFold API."""
    try:
        resp = requests.post(
            ESMFOLD_API_URL,
            data=sequence,
            timeout=ESMFOLD_TIMEOUT_SEC,
        )
        resp.raise_for_status()
        return resp.text
    except requests.RequestException as exc:
        st.error(f"ESMFold API error: {exc}")
        return None


def fetch_pdb_by_id(pdb_id: str) -> Optional[str]:
    """Download a PDB structure from RCSB via REST API."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return resp.text
    except requests.RequestException as exc:
        st.error(f"PDB fetch failed for '{pdb_id}': {exc}")
        return None


# =============================================================================
# Analysis logic
# =============================================================================

def _compute_rg(coords: np.ndarray) -> float:
    """Radius of gyration from CA coordinates."""
    centroid = coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1))))


def _smooth(arr: np.ndarray, window: int) -> np.ndarray:
    """Simple moving-average smoothing."""
    if len(arr) <= window:
        return arr
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="valid")


def analyze_protein(
    pdb_string: str,
    window: int = SMOOTHING_WINDOW_DEFAULT,
) -> AnalysisResult:
    """
    Parse a PDB string and compute Rg, hydrophobic belt width, and profile.

    Parameters
    ----------
    pdb_string : str
        PDB-format text.
    window : int
        Smoothing window size for the Z/hydrophobicity profile.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", io.StringIO(pdb_string))
    except Exception:
        empty = ProfileData(np.array([]), np.array([]), 0.0, 0.0)
        return AnalysisResult(0.0, 0.0, empty, np.array([]), valid=False)

    coords: list[np.ndarray] = []
    hydro: list[float] = []

    for model in structure:
        for residue in model.get_residues():
            res_name = residue.get_resname()
            if res_name in KD_SCALE and "CA" in residue:
                coords.append(residue["CA"].get_coord())
                hydro.append(KD_SCALE[res_name])

    if not coords:
        empty = ProfileData(np.array([]), np.array([]), 0.0, 0.0)
        return AnalysisResult(0.0, 0.0, empty, np.array([]), valid=False)

    coords_arr = np.array(coords)
    hydro_arr = np.array(hydro)
    rg = _compute_rg(coords_arr)

    # Sort along Z-axis and smooth
    z_coords = coords_arr[:, 2]
    sorted_idx = np.argsort(z_coords)
    z_smooth = _smooth(z_coords[sorted_idx], window)
    h_smooth = _smooth(hydro_arr[sorted_idx], window)

    # Identify hydrophobic belt boundaries
    belt_z = z_smooth[h_smooth > HYDROPHOBIC_THRESHOLD]
    if len(belt_z) > 0:
        b_min, b_max = float(np.min(belt_z)), float(np.max(belt_z))
    else:
        b_min, b_max = 0.0, 0.0

    profile = ProfileData(z_smooth, h_smooth, b_min, b_max)
    return AnalysisResult(rg, b_max - b_min, profile, coords_arr)


# =============================================================================
# Scaffold recommendation
# =============================================================================

def recommend_msp(estimated_diameter_nm: float) -> str:
    """Choose the smallest MSP scaffold that accommodates the estimated diameter."""
    for max_diam, name in MSP_CATALOG:
        if estimated_diameter_nm < max_diam:
            return name
    return MSP_CATALOG[-1][1]  # largest available


def compute_compatibility_score(belt_width: float, target_thickness: float) -> float:
    """
    Score [0, 1] representing how well the belt matches the target lipid thickness.
    1.0 = perfect match, 0.0 = completely mismatched.
    """
    if target_thickness == 0:
        return 0.0
    return max(0.0, 1.0 - abs(belt_width - target_thickness) / target_thickness)


# =============================================================================
# PDF report
# =============================================================================

def create_pdf(
    rg: float,
    belt_width: float,
    score: float,
    msp_name: str,
    est_diam: float,
    fig: plt.Figure,
    lipid_name: str,
    protein_name: str = "Custom Protein",
) -> bytes:
    """Generate a PDF analysis report and return raw bytes."""
    pdf = FPDF()
    pdf.add_page()

    # --- Title ---
    pdf.set_font("Arial", "B", 16)
    pdf.cell(0, 10, "GPCR-Nanodisc Analysis Report", ln=True, align="C")
    pdf.ln(5)
    pdf.set_font("Arial", size=10)
    pdf.cell(0, 8, f"Protein: {protein_name}  |  Lipid: {lipid_name}", ln=True, align="C")
    pdf.ln(8)

    # --- Metrics table ---
    pdf.set_font("Arial", "B", 11)
    col_w = 80
    pdf.cell(col_w, 9, "Metric", 1, 0, "C")
    pdf.cell(col_w, 9, "Value", 1, 1, "C")
    pdf.set_font("Arial", size=11)

    angstrom = "\u00c5"
    metrics = [
        ("Radius of Gyration (Rg)", f"{rg:.2f} {angstrom}"),
        ("Hydrophobic Belt Width", f"{belt_width:.2f} {angstrom}"),
        ("Compatibility Score", f"{score:.2f}"),
        ("Estimated Diameter", f"{est_diam:.2f} nm"),
        ("Recommended Scaffold", msp_name),
    ]
    for label, value in metrics:
        pdf.cell(col_w, 9, label, 1)
        pdf.cell(col_w, 9, value, 1)
        pdf.ln()

    pdf.ln(8)

    # --- Plot image ---
    tmp_path: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as tmp:
            tmp_path = tmp.name
            fig.savefig(tmp_path, dpi=150, bbox_inches="tight")
        pdf.image(tmp_path, x=10, w=190)
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

    return pdf.output(dest="S").encode("latin-1")


# =============================================================================
# Visualization helpers
# =============================================================================

def render_3d_view(pdb_data: str, profile: ProfileData, belt_width: float) -> None:
    """Render interactive 3D structure with membrane belt overlay."""
    center_z = (profile.belt_min + profile.belt_max) / 2.0

    pdb_escaped = (
        pdb_data
        .replace("\\", "\\\\")
        .replace("`", "\\`")
        .replace("${", "\\${")
    )

    html = f"""<!DOCTYPE html>
<html>
<head><style>
    body {{ margin:0; padding:0; overflow:hidden; }}
    #viewer3d {{ width:100%; height:450px; position:relative; }}
</style></head>
<body>
<div id="viewer3d"></div>
<script>
function initViewer() {{
    var viewer = $3Dmol.createViewer("viewer3d", {{backgroundColor: "white"}});
    viewer.addModel(`{pdb_escaped}`, "pdb");
    viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});
    viewer.addShape({{
        type: "box",
        center: {{x: 0, y: 0, z: {center_z:.2f}}},
        dimensions: {{w: 60, h: 60, d: {belt_width:.2f}}},
        color: "yellow",
        opacity: 0.25
    }});
    viewer.zoomTo();
    viewer.render();
}}

// Poll until 3Dmol is loaded
var script = document.createElement("script");
script.src = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js";
script.onload = initViewer;
document.head.appendChild(script);
</script>
</body>
</html>"""
    components.html(html, height=480, width=800, scrolling=False)


def plot_hydrophobic_profile(
    profile: ProfileData,
    lipid_name: str,
) -> plt.Figure:
    """Create a hydrophobicity-vs-Z-coordinate plot and return the Figure."""
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(
        profile.z_smooth, profile.h_smooth,
        color="#1f77b4", linewidth=1.5, label="Hydrophobicity",
    )
    ax.axvspan(
        profile.belt_min, profile.belt_max,
        color="#ffd700", alpha=0.25, label="Hydrophobic Belt",
    )
    ax.axhline(
        y=HYDROPHOBIC_THRESHOLD, linestyle="--", color="gray",
        alpha=0.5, label=f"Threshold ({HYDROPHOBIC_THRESHOLD})",
    )
    ax.set_title(f"Hydrophobic Profile \u2014 {lipid_name}", fontsize=13)
    ax.set_xlabel("Z-coordinate (\u00c5)")
    ax.set_ylabel("KD Hydrophobicity")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    return fig


# =============================================================================
# Streamlit UI
# =============================================================================

def _init_session_state() -> None:
    """Initialise session-state keys once."""
    if "pdb_data" not in st.session_state:
        st.session_state.pdb_data = None
    if "pdb_name" not in st.session_state:
        st.session_state.pdb_name = "Custom Protein"


def _sidebar_input() -> None:
    """Sidebar: structure input & lipid selection."""
    st.sidebar.header("1. Target Selection")
    input_method = st.sidebar.radio(
        "Input Method",
        ["Paste Sequence (ESMFold)", "Upload PDB File", "Enter PDB ID"],
    )

    if input_method == "Paste Sequence (ESMFold)":
        seq = st.sidebar.text_area(
            "Amino Acid Sequence", height=150,
            help="One-letter amino acid codes, no FASTA header.",
        )
        if st.sidebar.button("Predict Structure"):
            if not seq.strip():
                st.sidebar.warning("Please enter a sequence.")
            else:
                with st.spinner("Running ESMFold prediction..."):
                    result = predict_structure(seq.strip())
                    if result:
                        st.session_state.pdb_data = result
                        st.session_state.pdb_name = "ESMFold Prediction"

    elif input_method == "Upload PDB File":
        uploaded = st.sidebar.file_uploader("Upload .pdb file", type=["pdb"])
        if uploaded:
            st.session_state.pdb_data = uploaded.getvalue().decode("utf-8")
            st.session_state.pdb_name = uploaded.name

    elif input_method == "Enter PDB ID":
        pdb_id = st.sidebar.text_input("PDB ID (e.g. 2RH1)", max_chars=4)
        if st.sidebar.button("Fetch PDB"):
            if len(pdb_id) != 4:
                st.sidebar.warning("PDB ID must be exactly 4 characters.")
            else:
                with st.spinner(f"Fetching {pdb_id.upper()}..."):
                    result = fetch_pdb_by_id(pdb_id.upper())
                    if result:
                        st.session_state.pdb_data = result
                        st.session_state.pdb_name = f"PDB: {pdb_id.upper()}"

    st.sidebar.divider()
    st.sidebar.header("2. Lipid Settings")


def main() -> None:
    st.set_page_config(page_title="GPCR-Nanodisc Predictor", layout="wide")
    st.title("\U0001f9ec GPCR-Nanodisc Integration Predictor")
    st.caption(
        "Evaluate GPCR\u2013nanodisc compatibility via Rg sizing, "
        "hydrophobic belt profiling, and scaffold recommendation."
    )

    _init_session_state()
    _sidebar_input()

    lipid_choice = st.sidebar.selectbox("Target Lipid", list(LIPID_DB.keys()))
    target_thickness = LIPID_DB[lipid_choice]

    st.sidebar.divider()
    with st.sidebar.expander("Advanced Settings"):
        smooth_window = st.slider(
            "Smoothing Window",
            min_value=3, max_value=30, value=SMOOTHING_WINDOW_DEFAULT,
            help="Moving-average window for the hydrophobicity profile.",
        )

    # ----- Main panel -----
    pdb_data = st.session_state.pdb_data
    if pdb_data is None:
        st.info("Select a protein structure from the sidebar to begin analysis.")
        return

    result = analyze_protein(pdb_data, window=smooth_window)

    if not result.valid:
        st.error(
            "Could not parse any residues from the structure. "
            "Please check the PDB data."
        )
        return

    # -- Metrics --
    score = compute_compatibility_score(result.belt_width, target_thickness)
    est_diam = result.rg * RG_TO_DIAMETER_FACTOR / 10.0
    msp_name = recommend_msp(est_diam)

    st.subheader("Analysis Results")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Rg", f"{result.rg:.2f} \u00c5")
    c2.metric("Belt Width", f"{result.belt_width:.2f} \u00c5")
    c3.metric("Score", f"{score:.2f}")
    c4.metric("Est. Diameter", f"{est_diam:.1f} nm")

    st.info(f"**Recommended Scaffold:** {msp_name}")

    # -- 3D view --
    st.subheader("\U0001f9ca 3D Visualization")
    render_3d_view(pdb_data, result.profile, result.belt_width)

    # -- Profile plot --
    st.subheader("\U0001f4ca Hydrophobic Profile")
    fig = plot_hydrophobic_profile(result.profile, lipid_choice)
    st.pyplot(fig)

    # -- PDF export --
    st.divider()
    st.subheader("\U0001f4e5 Export Report")
    pdf_bytes = create_pdf(
        rg=result.rg,
        belt_width=result.belt_width,
        score=score,
        msp_name=msp_name,
        est_diam=est_diam,
        fig=fig,
        lipid_name=lipid_choice,
        protein_name=st.session_state.pdb_name,
    )
    st.download_button(
        label="Download Analysis Report (PDF)",
        data=pdf_bytes,
        file_name=f"nanodisc_report_{msp_name}.pdf",
        mime="application/pdf",
    )


if __name__ == "__main__":
    main()
