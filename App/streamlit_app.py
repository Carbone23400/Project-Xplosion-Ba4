import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "data"))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from coordchem.parser import parse_formula, FormulaParseError
from coordchem.geometry import geometry_report
from coordchem.complex import Complex
from coordchem.spectra.predictor import predict_spectrum
from coordchem.spectra.renderer import plot_spectrum 

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

def get_iupac_name(parsed) -> str:
    """
    Build a basic IUPAC name from a ParsedComplex object.
    e.g. [Fe(CN)6]4- → hexacyanidoferrate(II)
    """
    # Number prefixes
    PREFIXES = {
        1: "mono", 2: "di",    3: "tri",   4: "tetra",
        5: "penta", 6: "hexa", 7: "hepta", 8: "octa",
        9: "nona",  10: "deca"
    }

    # Build ligand part — alphabetical order by ligand name (IUPAC rule)
    ligand_parts = []
    for lig, count in sorted(
        parsed.ligands.items(),
        key=lambda x: parsed.ligand_names.get(x[0], x[0])
    ):
        name   = parsed.ligand_names.get(lig, lig)
        prefix = PREFIXES.get(count, str(count))
        # Use bis/tris/tetrakis for complex ligand names containing a number
        if any(p in name for p in ("di", "tri", "tetra", "bis", "tris")):
            complex_prefixes = {2: "bis", 3: "tris", 4: "tetrakis",
                                5: "pentakis", 6: "hexakis"}
            prefix = complex_prefixes.get(count, str(count))
            ligand_parts.append(f"{prefix}({name})")
        else:
            ligand_parts.append(f"{prefix}{name}")

    ligand_str = "".join(ligand_parts)

    # Metal name — lowercase
    METAL_NAMES = {
        "Fe": "ferrate", "Cu": "cuprate", "Co": "cobaltate",
        "Ni": "nickelate", "Cr": "chromate", "Mn": "manganate",
        "Pt": "platinate", "Pd": "palladate", "Au": "aurate",
        "Ag": "argentate", "Zn": "zincate",  "Mo": "molybdate",
        "W":  "tungstate",  "Ru": "ruthenate","Os": "osmate",
        "Rh": "rhodate",    "Ir": "iridate",  "V":  "vanadate",
        "Ti": "titanate",   "Cd": "cadmate",
    }

    ox = parsed.oxidation_state
    ox_str     = f"({_roman(ox)})" if ox is not None else ""
    metal_name = METAL_NAMES.get(parsed.metal, parsed.metal.lower() + "ate")

    # If complex charge is positive or zero, use metal name directly (cation)
    if parsed.complex_charge >= 0:
        metal_name = parsed.metal.lower()

    return f"{ligand_str}{metal_name}{ox_str}"


def _roman(n: int) -> str:
    """Convert an integer to a Roman numeral string."""
    if n is None:
        return ""
    if n == 0:
        return "0"
    vals  = [10,"X", 9,"IX", 5,"V", 4,"IV", 1,"I"]
    vals  = [(10,"X"),(9,"IX"),(5,"V"),(4,"IV"),(1,"I")]
    neg   = n < 0
    n     = abs(n)
    result = ""
    for value, numeral in vals:
        while n >= value:
            result += numeral
            n      -= value
    return f"−{result}" if neg else result

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

def render_2d_view(complex_obj: Complex, size: int = 520):
    """Render a 2D SVG depiction of ``complex_obj`` inside the current Streamlit page."""
    import streamlit.components.v1 as components

    try:
        svg = complex_obj.draw_2d_svg(size=size, title="")
    except Exception as exc:
        st.warning(f"Could not build a 2D drawing: {exc}")
        return
    components.html(svg, width=size, height=size + 20)

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
# App layout
# ---------------------------------------------------------------------------

st.set_page_config(page_title="CoordAnalyst", page_icon="⚛️", layout="wide")
st.title("CoordAnalyst : Coordination Complex Spectra Predictor")
st.caption("±20–50 cm⁻¹ accuracy · data from Nakamoto 6th ed.")

with st.sidebar:
    st.header("Complex Input")

    input_mode = st.radio("Input mode", ["Formula", "Name"], horizontal=True)

    if input_mode == "Formula":
        user_input = st.text_input("Enter a coordination complex formula", value="[Fe(CN)6]4-")
    else:
        user_input = st.text_input("Enter a coordination complex name", value="hexacyanoferrate(II)")

    analyze = st.button("Analyze", type="primary", use_container_width=True)

    st.divider()

    spectrum_type = st.radio("Spectrum type", ["IR", "Raman", "Both"], horizontal=True)
    sigma         = st.slider("Peak width σ (cm⁻¹)", min_value=5, max_value=60, value=20, step=5)

    st.divider()
    st.caption(
        "**Supported ligands:** CN CO NH₃ H₂O Cl Br I OH NO₂ ONO SCN NCS "
        "NO N₃ en ox acac py dmso PPh₃"
    )

if not user_input.strip():
    st.info("Enter a formula or name in the sidebar and click **Analyze**.")
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
    label = "formula" if input_mode == "Formula" else "name"
    st.error(f"Could not resolve {label}: {e}")
    st.stop()

for w in parsed.warnings:
    st.warning(w)

# Complex info
if input_mode == "Formula":
    st.subheader(f"Formula : {user_input}")
    st.subheader(f"Name: {get_iupac_name(parsed)}")
else:
    st.subheader(f"Formula: {parsed.raw_formula}")
    st.subheader(f"Name : {user_input}")

st.subheader("Complex Information")
c1, c2, c3, c4, c5 = st.columns(5)
ox = parsed.oxidation_state
c1.metric("Metal",           parsed.metal)
c2.metric("Oxidation State", f"+{ox}" if ox and ox > 0 else str(ox))
c3.metric("Coord. Number",   parsed.coordination_number)
c4.metric("d-count",         report["d_count"] if report["d_count"] is not None else "—")
c5.metric("Geometry",        report["geometry"])

with st.expander("Ligand details", expanded=True):
    lig_rows = []
    for lig, count in parsed.ligands.items():
        charge = parsed.ligand_charges.get(lig, 0)
        lig_rows.append({
            "Formula":    lig,
            "Count":      str(count),
            "Name":       parsed.ligand_names.get(lig, lig),
            "Charge":     f"+{charge}" if charge > 0 else str(charge),
            "Denticity":  str(parsed.ligand_denticity.get(lig, 1)),
            "Donor atom": parsed.donor_atoms.get(lig, "?"),
        })
    st.dataframe(pd.DataFrame(lig_rows), use_container_width=True, hide_index=True)

complex_obj = Complex(parsed)
with st.expander("2D Diagram", expanded=True):
    render_2d_view(complex_obj)
with st.expander("3D Structure", expanded=True):
    render_3d_view(complex_obj)

st.divider()

# Run predictor(s)
want_ir    = spectrum_type in ("IR", "Both")
want_raman = spectrum_type in ("Raman", "Both")

ir_result    = predict_spectrum(parsed, spectrum_type="IR",    apply_corrections=True) if want_ir    else None
raman_result = predict_spectrum(parsed, spectrum_type="Raman", apply_corrections=True) if want_raman else None

for r in (ir_result, raman_result):
    if r is not None:
        for w in r.warnings:
            st.warning(w)


def _plot_block(result, label):
    if result.bands:
        st.plotly_chart(
            plot_spectrum(result.bands, result.intensities, f"{label} — {user_input}", sigma),
            use_container_width=True,
        )
    else:
        st.info(f"No {label} data available for this complex.")


# Spectrum plots
st.subheader("Predicted Spectra")
if spectrum_type == "Both":
    col_ir, col_ra = st.columns(2)
    with col_ir:
        _plot_block(ir_result, "IR")
    with col_ra:
        _plot_block(raman_result, "Raman")
else:
    _plot_block(ir_result if want_ir else raman_result, spectrum_type)

# Band assignments
st.subheader("Band Assignments")
if spectrum_type == "Both":
    tab_ir, tab_ra = st.tabs(["IR Bands", "Raman Bands"])
    with tab_ir:
        if ir_result and ir_result.bands:
            st.dataframe(bands_to_df(ir_result), use_container_width=True, hide_index=True)
        else:
            st.info("No IR band data available.")
    with tab_ra:
        if raman_result and raman_result.bands:
            st.dataframe(bands_to_df(raman_result), use_container_width=True, hide_index=True)
        else:
            st.info("No Raman band data available.")
else:
    active = ir_result if want_ir else raman_result
    if active and active.bands:
        st.dataframe(bands_to_df(active), use_container_width=True, hide_index=True)
    else:
        st.info(f"No {spectrum_type} band data available.")

st.divider()
st.caption("Data: Nakamoto, K. *Infrared and Raman Spectra of Inorganic and Coordination Compounds*, 6th ed., Wiley, 2009.")