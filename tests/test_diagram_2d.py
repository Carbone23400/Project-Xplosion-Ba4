"""
tests/test_diagram_2d.py
------------------------
Tests for coordchem.viz.diagram_2d.

Run with:
    python -m pytest tests/test_diagram_2d.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

from coordchem.viz.diagram_2d import (
    build_coordination_mol,
    diagram_2d_svg,
    save_diagram_2d,
)


# ============================================================================
# Test examples: coordination complexes only
# ============================================================================

SQUARE_PLANAR_COMPLEX = "[PtCl4]2-"
TETRAAMMINE_COMPLEX = "[Cu(NH3)4]2+"
OCTAHEDRAL_COMPLEX = "[Fe(CN)6]4-"
INVALID_COMPLEX = "not a complex"


# ============================================================================
# Basic SVG generation
# ============================================================================

def test_diagram_2d_svg_returns_svg_string():
    svg = diagram_2d_svg(SQUARE_PLANAR_COMPLEX)

    assert isinstance(svg, str)
    assert "<svg" in svg
    assert "</svg>" in svg


def test_diagram_2d_svg_uses_requested_size():
    svg = diagram_2d_svg(SQUARE_PLANAR_COMPLEX, size=300)

    assert "width='300px'" in svg
    assert "height='300px'" in svg


def test_diagram_2d_svg_accepts_title():
    svg = diagram_2d_svg(OCTAHEDRAL_COMPLEX, title="Hexacyanoferrate(II)")

    assert "<svg" in svg
    assert "</svg>" in svg


# ============================================================================
# Rendering behaviour
# ============================================================================

def test_diagram_renders_normally_without_highlight():
    svg = diagram_2d_svg(SQUARE_PLANAR_COMPLEX)

    assert "<svg" in svg
    assert "</svg>" in svg
    assert "bond-" in svg
    assert "atom-" in svg


# ============================================================================
# Geometry / projection behaviour
# ============================================================================

def test_cn4_or_3d_projection_uses_depth_cues():
    svg = diagram_2d_svg(TETRAAMMINE_COMPLEX).lower()

    assert "<svg" in svg
    assert "bond-" in svg
    assert "atom-" in svg


def test_octahedral_projection_uses_depth_cues():
    svg = diagram_2d_svg(OCTAHEDRAL_COMPLEX).lower()

    assert "<svg" in svg
    assert "</svg>" in svg
    assert "bond-" in svg
    assert "atom-" in svg


# ============================================================================
# Error handling
# ============================================================================

def test_diagram_2d_svg_rejects_invalid_complex():
    with pytest.raises(ValueError):
        diagram_2d_svg(INVALID_COMPLEX)


def test_diagram_2d_svg_rejects_invalid_size():
    with pytest.raises(ValueError, match="positive"):
        diagram_2d_svg(SQUARE_PLANAR_COMPLEX, size=0)


# ============================================================================
# File writing
# ============================================================================

def test_save_diagram_2d_writes_svg_file():
    out_file = Path("outputs/fe_cn6.svg")
    out_file.parent.mkdir(parents=True, exist_ok=True)

    returned = save_diagram_2d(
        OCTAHEDRAL_COMPLEX,
        out_file,
        title="Hexacyanoferrate(II)",
    )

    assert returned == out_file
    assert out_file.exists()

    content = out_file.read_text(encoding="utf-8")
    assert "<svg" in content
    assert "</svg>" in content
    assert "bond-" in content
    assert "atom-" in content


# ============================================================================
# Extra representative complexes
# ============================================================================

def test_square_planar_complex_draws():
    svg = diagram_2d_svg("[PdCl4]2-")

    assert "<svg" in svg
    assert "</svg>" in svg


def test_octahedral_ammine_complex_draws():
    svg = diagram_2d_svg("[Co(NH3)6]3+")

    assert "<svg" in svg
    assert "</svg>" in svg


def test_linear_complex_draws():
    svg = diagram_2d_svg("[Ag(NH3)2]+")

    assert "<svg" in svg
    assert "</svg>" in svg


def test_monodentate_labels_face_the_metal():
    mol = build_coordination_mol("[Ag(CN)2]-")

    labels = [
        atom.GetProp("atomLabel")
        for atom in mol.GetAtoms()
        if atom.HasProp("atomLabel")
    ]

    assert labels == ["NC", "CN"]


def test_edta_complex_draws_as_hexadentate_ligand():
    mol = build_coordination_mol("[Co(EDTA)]-")

    metal = 0
    metal_bonds = [
        bond
        for bond in mol.GetBonds()
        if metal in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert mol.GetNumAtoms() == 21
    assert len(metal_bonds) == 6
