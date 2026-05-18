import sys
import os
from dataclasses import dataclass

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "data"))

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from coordchem.parser import KNOWN_LIGANDS, parse_formula, FormulaParseError
from coordchem.geometry import geometry_report
from coordchem.complex import Complex
from coordchem.spectra.predictor import predict_spectrum
from coordchem.spectra.renderer import plot_spectrum 

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


@dataclass(frozen=True)
class GeometryChoice:
    """User-facing geometry label paired with the renderer-compatible value."""
    label: str
    value: str


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


def render_centered_html(html: str, height: int):
    """Render HTML centered inside a full-width Streamlit component."""
    import streamlit.components.v1 as components

    components.html(
        f"""
        <div style="width:100%; display:flex; justify-content:center; align-items:center;">
            <div style="flex:0 0 auto; max-width:100%;">
                {html}
            </div>
        </div>
        """,
        height=height,
        scrolling=False,
    )


def render_2d_view(
    complex_obj: Complex,
    size: int = 520,
    geometry_override: str | None = None,
):
    """Render a 2D SVG depiction of ``complex_obj`` inside the current Streamlit page."""
    from coordchem.viz.diagram_2d import diagram_2d_svg

    try:
        svg = diagram_2d_svg(
            complex_obj.parsed,
            size=size,
            title="",
            geometry_override=geometry_override,
        )
    except Exception as exc:
        st.warning(f"Could not build a 2D drawing: {exc}")
        return
    render_centered_html(svg, height=size + 20)

def render_3d_view(
    complex_obj: Complex,
    width: int = 500,
    height: int = 400,
    geometry_override: str | None = None,
):
    """Render a 3D view of ``complex_obj`` inside the current Streamlit page."""
    from coordchem.viz.molecule3D import complex_3d_html

    try:
        html = complex_3d_html(
            complex_obj.parsed,
            width=width,
            height=height,
            geometry=geometry_override,
        )
    except Exception as exc:
        st.warning(f"Could not build a 3D view: {exc}")
        return
    render_centered_html(html, height=height + 20)


def supported_ligands_text() -> str:
    """Return ligands in the same order as the parser's KNOWN_LIGANDS table."""
    return " ".join(KNOWN_LIGANDS)


def geometry_choices(predicted_geometry: str) -> list[GeometryChoice]:
    """Return drawable geometry choices from a possibly ambiguous prediction."""
    choices = [
        GeometryChoice(
            label=part.strip(),
            value=normalized_geometry_choice(part.strip()),
        )
        for part in predicted_geometry.split(" or ")
        if part.strip()
    ]
    if not choices:
        choices = [
            GeometryChoice(
                label=predicted_geometry,
                value=normalized_geometry_choice(predicted_geometry),
            )
        ]
    return list(dict.fromkeys(choices))


def normalized_geometry_choice(geometry: str) -> str:
    """Normalize explanatory geometry labels to values used by renderers/rules."""
    return geometry


def geometry_label(geometry: str | GeometryChoice) -> str:
    """Human-readable label for geometry selection controls."""
    if isinstance(geometry, GeometryChoice):
        geometry = geometry.label
    return geometry[:1].upper() + geometry[1:]


def dmso_donor_choices(parsed) -> list[str]:
    """Return DMSO donor choices when HSAB leaves the binding mode ambiguous."""
    donor_info = parsed.donor_atoms.get("dmso")
    if donor_info is None:
        return []

    donors = [
        donor.strip()
        for donor in str(donor_info).split("/")
        if donor.strip()
    ]
    return donors if len(donors) > 1 else []


def dmso_donor_label(donor: str) -> str:
    """Display label for DMSO linkage choices."""
    labels = {
        "S": "S-bound DMSO",
        "O": "O-bound DMSO",
    }
    return labels.get(donor, donor)


def dmso_assignment_note(parsed) -> str | None:
    """Return a concise professional note for non-ambiguous DMSO assignments."""
    if "dmso" not in parsed.ligands:
        return None

    donor = parsed.donor_atoms.get("dmso")
    if donor == "O":
        return "DMSO linkage: O-bound, estimated from HSAB metal hardness."
    if donor == "S":
        return "DMSO linkage: S-bound, estimated from HSAB metal softness."
    return None


def _plot_block(result, label, sigma, user_input, invert=False):
    if result.bands:
        st.plotly_chart(
            plot_spectrum(result.bands, result.intensities, f"{label} — {user_input}", sigma, invert=invert),
            use_container_width=True,
        )
    else:
        st.info(f"No {label} data available for this complex.")

# ---------------------------------------------------------------------------
# App layout
# ---------------------------------------------------------------------------

st.set_page_config(page_title="CoordAnalyst", page_icon="⚛️", layout="wide")
st.markdown(
    """
    <style>
    div[data-testid="stMetricValue"] {
        white-space: normal;
        overflow-wrap: break-word;
        word-break: normal;
        line-height: 1.15;
    }
    div[data-testid="stMetricValue"] > div {
        white-space: normal;
        overflow-wrap: break-word;
        word-break: normal;
    }
    </style>
    """,
    unsafe_allow_html=True,
)
st.title("CoordAnalyst : Coordination Complex Spectra Predictor")
st.caption("±20–50 cm⁻¹ accuracy · data from Nakamoto 6th ed.")

with st.sidebar:
    st.header("Complex Input")
    input_mode = st.radio("Input mode", ["Formula", "Name"], horizontal=True)
    if input_mode == "Formula":
        user_input = st.text_input("Enter a coordination complex formula", value="[Fe(CN)6]4-")
    else:
        user_input = st.text_input("Enter a coordination complex name", value="hexacyanoferrate(II)")
    analyze       = st.button("Analyze", type="primary", use_container_width=True)
    st.divider()
    spectrum_type = st.radio("Spectrum type", ["IR", "Raman", "Both"], index=2, horizontal=True)
    sigma         = st.slider("Peak width σ (cm⁻¹)", min_value=5, max_value=60, value=20, step=5)
    if spectrum_type in ("IR", "Both"):
        invert    = st.checkbox("Transmittance style (downward peaks)", value=False, key="invert_ir")
    else:
        invert    = False


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

predicted_geometry = str(report["geometry"])
available_geometries = geometry_choices(predicted_geometry)
selected_geometry_choice = available_geometries[0]

if len(available_geometries) > 1:
    st.info(
        "Several conformations are plausible for this complex. "
        "Choose the one to use for the drawings and spectrum rules."
    )
    selected_geometry_choice = st.radio(
        "Conformation",
        available_geometries,
        format_func=geometry_label,
        horizontal=True,
    )

selected_geometry = selected_geometry_choice.value
selected_geometry_label = selected_geometry_choice.label
parsed.geometry = selected_geometry

available_dmso_donors = dmso_donor_choices(parsed)
selected_dmso_donor = None

if available_dmso_donors:
    st.info(
        "DMSO can bind through sulfur or oxygen for this metal. "
        "Choose the linkage to use for the drawings and spectra."
    )
    selected_dmso_donor = st.radio(
        "DMSO linkage",
        available_dmso_donors,
        format_func=dmso_donor_label,
        horizontal=True,
    )
    parsed.donor_atoms["dmso"] = selected_dmso_donor

for w in parsed.warnings:
    if w.startswith("DMSO "):
        continue
    st.warning(w)

dmso_note = dmso_assignment_note(parsed)
if dmso_note and not available_dmso_donors:
    st.info(dmso_note)

# Complex identity
st.subheader("Complex Identity")
name_col, formula_col = st.columns(2)
name_col.markdown(f"**Name:** {parsed.iupac_name}")
formula_col.markdown(f"**Formula:** {parsed.raw_formula}")

st.subheader("Complex Information")
c1, c2, c3, c4, c5 = st.columns([1, 1, 1, 1, 1.7])
ox = parsed.oxidation_state
c1.metric("Metal",           parsed.metal)
c2.metric("Oxidation State", f"+{ox}" if ox and ox > 0 else str(ox))
c3.metric("Coord. Number",   parsed.coordination_number)
c4.metric("d-count",         report["d_count"] if report["d_count"] is not None else "—")
c5.metric("Geometry",        geometry_label(selected_geometry_label))
if selected_geometry_label != predicted_geometry:
    st.caption(f"Predicted geometry: {predicted_geometry}")

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
    render_2d_view(complex_obj, geometry_override=selected_geometry)
with st.expander("3D Structure", expanded=True):
    render_3d_view(complex_obj, geometry_override=selected_geometry)

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


# Spectrum plots
st.subheader("Predicted Spectra")
if spectrum_type == "Both":
    col_ir, col_ra = st.columns(2)
    with col_ir:
        _plot_block(ir_result,    "IR",    sigma, user_input, invert=invert)
    with col_ra:
        _plot_block(raman_result, "Raman", sigma, user_input, invert=False)
elif want_ir:
    _plot_block(ir_result,    "IR",    sigma, user_input, invert=invert)
else:
    _plot_block(raman_result, "Raman", sigma, user_input, invert=False)

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
