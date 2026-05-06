import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "data"))

import numpy as np
import pandas as pd
import streamlit as st

from coordchem.parser import parse_formula, FormulaParseError
from coordchem.geometry import geometry_report
from coordchem.spectra.predictor import predict_spectrum
from coordchem.spectra.renderer import plot_spectrum, build_spectrum
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

def bands_to_df(result):
    """Convert a PredictionResult to a display DataFrame."""
    rows = []
    for b in sorted(result.bands, key=lambda b: b.center):
        rows.append({
            "Ligand":         b.ligand,
            "Mode":           b.coordination,
            "Center (cm⁻¹)": f"{b.center:.0f}",
            "Range (cm⁻¹)":  f"{b.wn_min:.0f}–{b.wn_max:.0f}",
            "Intensity":      b.intensity,
            "Assignment":     b.assignment,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------

st.set_page_config(page_title="CoordAnalyst", page_icon="⚛️", layout="centered")
st.title("⚛️ CoordAnalyst — Coordination Complex Spectra Predictor")
st.caption("Educational tool · ±20–50 cm⁻¹ accuracy · data from Nakamoto 6th ed.")

input_mode = st.radio("Input mode", ["Formula", "Name"], horizontal=True)

if input_mode == "Formula":
    user_input = st.text_input("Enter a coordination complex formula", value="[Fe(CN)6]4-")
else:
    user_input = st.text_input("Enter a complex name", value="hexacyanoferrate(II)")

spectrum_type = st.radio("Spectrum type", ["IR", "Raman"], horizontal=True)
sigma         = st.slider("Peak width σ (cm⁻¹)", min_value=5, max_value=60, value=20, step=5)

if not st.button("Analyze", type="primary"):
    st.stop()

# Parse formula
from coordchem.name import parse_name   

try:
    if input_mode == "Formula":
        parsed = parse_formula(user_input.strip())
    else:
        parsed = parse_name(user_input.strip())
    report = geometry_report(parsed)
except FormulaParseError as e:
    st.error(f"Could not parse formula: {e}")
    st.stop()
except Exception as e:
    st.error(f"Could not resolve name: {e}")
    st.stop()

for w in parsed.warnings:
    st.warning(w)

# Complex info
st.subheader("Complex Information")
ox = parsed.oxidation_state
st.write(
    f"**Metal:** {parsed.metal}  |  "
    f"**Oxidation state:** {'+' + str(ox) if ox and ox > 0 else ox}  |  "
    f"**Coordination number:** {parsed.coordination_number}  |  "
    f"**d-count:** {report['d_count'] if report['d_count'] is not None else '—'}  |  "
    f"**Geometry:** {report['geometry']}"
)

st.markdown("**Ligands:**")
for lig, count in parsed.ligands.items():
    name   = parsed.ligand_names.get(lig, lig)
    charge = parsed.ligand_charges.get(lig, 0)
    dent   = parsed.ligand_denticity.get(lig, 1)
    donor  = parsed.donor_atoms.get(lig, "?")
    st.write(
        f"- {count}× **{lig}** ({name}) — "
        f"charge {charge:+d}, denticity {dent}, donor atom: {donor}"
    )

# Run predictor
result    = predict_spectrum(parsed, spectrum_type=spectrum_type)
all_bands = result.bands

for w in result.warnings:
    st.warning(w)

# Spectrum
st.subheader(f"Predicted {spectrum_type} Spectrum")
if all_bands:
    st.plotly_chart(
        plot_spectrum(result.bands, result.intensities, f"{spectrum_type} — {user_input}", sigma),
        use_container_width=True
    )
    st.subheader("Band Assignments")
    st.dataframe(bands_to_df(result), use_container_width=True, hide_index=True)
else:
    st.info(f"No {spectrum_type} data available for the ligands in this complex.")

st.divider()
st.caption("Data: Nakamoto, K. *Infrared and Raman Spectra of Inorganic and Coordination Compounds*, 6th ed., Wiley, 2009.")