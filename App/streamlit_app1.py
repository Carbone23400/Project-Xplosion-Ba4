import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from coordchem.parser import parse_formula, FormulaParseError
from coordchem.geometry import geometry_report
from coordchem.complex import Complex
from data.database.ir_ra_bands import IRBandDB


def render_3d_view(complex_obj: Complex, width: int = 500, height: int = 400):
    """Render a 3D view of ``complex_obj`` inside the current Streamlit page."""
    import streamlit.components.v1 as components

    try:
        html = complex_obj.draw_3d_html(width=width, height=height)
    except Exception as exc:
        st.warning(f"Could not build a 3D view: {exc}")
        return
    components.html(html, width=width, height=height + 20)


# ---------------------------------------------------------------------------
# Spectrum helpers
# ---------------------------------------------------------------------------

INTENSITY_SCALE = {
    "very strong":  1.00,
    "strong":       0.75,
    "strong broad": 0.70,
    "medium":       0.50,
    "medium broad": 0.45,
    "weak":         0.25,
}


def gaussian(x, center, height, sigma):
    """Return a Gaussian curve centered at `center` with given height and width."""
    return height * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def build_spectrum(bands, sigma=20.0, wn_range=(100, 4000)):
    """
    Combine all bands into a single Gaussian-broadened spectrum.
    Returns (x, y) arrays normalized to [0, 1].
    """
    x = np.linspace(wn_range[0], wn_range[1], 3000)
    y = np.zeros_like(x)
    for band in bands:
        height = INTENSITY_SCALE.get(band.intensity.lower().strip(), 0.5)
        y += gaussian(x, band.center, height, sigma)
    if y.max() > 0:
        y /= y.max()
    return x, y


def plot_spectrum(bands, title, sigma):
    """Build a Plotly figure from a list of BandRecord objects."""
    x, y = build_spectrum(bands, sigma=sigma)
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x, y=y,
        mode="lines",
        fill="tozeroy",
        line=dict(color="#2563eb", width=2),
        fillcolor="rgba(37,99,235,0.08)",
    ))
    for band in bands:
        fig.add_vline(x=band.center, line=dict(color="rgba(100,116,139,0.3)", width=1, dash="dot"))
    fig.update_layout(
        title=title,
        xaxis_title="Wavenumber (cm⁻¹)",
        yaxis_title="Relative Intensity",
        xaxis=dict(autorange="reversed", showgrid=True, gridcolor="#f1f5f9"),
        yaxis=dict(range=[0, 1.08], showgrid=True, gridcolor="#f1f5f9"),
        plot_bgcolor="white",
        height=380,
        margin=dict(l=50, r=20, t=50, b=45),
    )
    return fig


def bands_to_df(bands):
    """Convert a list of BandRecord objects to a display DataFrame."""
    rows = []
    for b in sorted(bands, key=lambda b: b.center):
        rows.append({
            "Ligand":        b.ligand,
            "Mode":          b.coordination,
            "Center (cm⁻¹)": f"{b.center:.0f}",
            "Range (cm⁻¹)":  f"{b.wn_min:.0f}–{b.wn_max:.0f}",
            "Intensity":     b.intensity,
            "Assignment":    b.assignment,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------

st.set_page_config(page_title="CoordChem", page_icon="🔬", layout="centered")

st.title("🔬 CoordChem — IR/Raman Spectrum Predictor")
st.caption("Educational tool · ±20–50 cm⁻¹ accuracy · data from Nakamoto 6th ed.")

formula_input = st.text_input("Enter a coordination complex formula", value="[Fe(CN)6]4-")
spectrum_type = st.radio("Spectrum type", ["IR", "Raman"], horizontal=True)
sigma = st.slider("Peak width σ (cm⁻¹)", min_value=5, max_value=60, value=20, step=5)

if not st.button("Analyze", type="primary"):
    st.stop()

# Parse formula
try:
    parsed = parse_formula(formula_input.strip())
    report = geometry_report(parsed)
except FormulaParseError as e:
    st.error(f"Could not parse formula: {e}")
    st.stop()

for w in parsed.warnings:
    st.warning(w)

# Complex info
st.subheader("Complex Information")
ox = parsed.oxidation_state
st.write(f"**Metal:** {parsed.metal}  |  "
         f"**Oxidation state:** {'+' + str(ox) if ox and ox > 0 else ox}  |  "
         f"**Coordination number:** {parsed.coordination_number}  |  "
         f"**d-count:** {report['d_count'] if report['d_count'] is not None else '—'}  |  "
         f"**Geometry:** {report['geometry']}")

st.markdown("**Ligands:**")
for lig, count in parsed.ligands.items():
    name   = parsed.ligand_names.get(lig, lig)
    charge = parsed.ligand_charges.get(lig, 0)
    dent   = parsed.ligand_denticity.get(lig, 1)
    donor  = parsed.donor_atoms.get(lig, "?")
    st.write(f"- {count}× **{lig}** ({name}) — charge {charge:+d}, denticity {dent}, donor atom: {donor}")

# 3D view
st.subheader("3D Structure")
render_3d_view(Complex(parsed))

# Fetch bands
db = IRBandDB()
all_bands = []
for lig in parsed.ligands:
    all_bands += db.get_bands(lig, spectrum_type=spectrum_type, metal=parsed.metal)
db.close()

# Spectrum
st.subheader(f"Predicted {spectrum_type} Spectrum")
if all_bands:
    st.plotly_chart(plot_spectrum(all_bands, f"{spectrum_type} — {formula_input}", sigma), use_container_width=True)
    st.subheader("Band Assignments")
    st.dataframe(bands_to_df(all_bands), use_container_width=True, hide_index=True)
else:
    st.info(f"No {spectrum_type} data available for the ligands in this complex.")

st.divider()
st.caption("Data: Nakamoto, K. *Infrared and Raman Spectra of Inorganic and Coordination Compounds*, 6th ed., Wiley, 2009.")
