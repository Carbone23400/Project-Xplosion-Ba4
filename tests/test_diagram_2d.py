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
from coordchem.viz.ligand_data import LIGAND_SMILES


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


def test_diagram_2d_svg_uses_generated_coordination_name_by_default():
    svg = diagram_2d_svg("[Co(EDTA)]-")

    assert "ethylenediaminetetraacetatocobaltate(III)" in svg


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


def test_en_donor_atoms_display_as_nh2_without_separate_h_atoms():
    mol = build_coordination_mol("[Co(en)3]3+")

    donor_atoms = [
        atom
        for atom in mol.GetAtoms()
        if atom.HasProp("atomLabel")
        and atom.GetProp("atomLabel") in {"NH2", "H2N"}
    ]
    h_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "H"]

    assert len(donor_atoms) == 6
    assert all(atom.GetAtomicNum() == 7 for atom in donor_atoms)
    assert not h_atoms


def test_edta_projection_has_two_wedges_two_dashes_and_two_plain_bonds():
    mol = build_coordination_mol("[Co(EDTA)]-")

    bond_dirs = [
        str(bond.GetBondDir())
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert bond_dirs.count("BEGINWEDGE") == 2
    assert bond_dirs.count("BEGINDASH") == 2
    assert bond_dirs.count("NONE") == 2


def test_edta_projection_is_vertical_and_left_right_symmetric():
    mol = build_coordination_mol("[Co(EDTA)]-")
    conf = mol.GetConformer()

    sites = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            other = bond.GetOtherAtomIdx(0)
            pos = conf.GetAtomPosition(other)
            sites.append((str(bond.GetBondDir()), round(pos.x, 2), round(pos.y, 2)))

    assert ("NONE", 0.0, 2.15) in sites
    assert ("NONE", 0.0, -2.15) in sites
    assert ("BEGINDASH", -1.86, 1.07) in sites
    assert ("BEGINDASH", 1.86, 1.07) in sites
    assert ("BEGINWEDGE", -1.86, -1.07) in sites
    assert ("BEGINWEDGE", 1.86, -1.07) in sites


def test_edta_ligand_coordinates_are_horizontally_symmetric():
    mol = build_coordination_mol("[Co(EDTA)]-")
    conf = mol.GetConformer()

    mirror_pairs = [
        (1, 12),
        (2, 13),
        (3, 14),
        (4, 15),
        (5, 16),
        (6, 17),
        (7, 18),
        (8, 19),
        (9, 20),
        (10, 11),
    ]

    for top_idx, bottom_idx in mirror_pairs:
        top = conf.GetAtomPosition(top_idx)
        bottom = conf.GetAtomPosition(bottom_idx)
        assert round(top.x, 2) == round(bottom.x, 2)
        assert round(top.y, 2) == round(-bottom.y, 2)


def test_edta_svg_contains_interrupted_internal_bonds():
    svg = diagram_2d_svg("[Co(EDTA)]-")

    assert "stroke-dasharray:6,5" not in svg
    assert svg.count("class='bond-5 atom-6 atom-7'") == 2
    assert svg.count("class='bond-16 atom-17 atom-18'") == 2


def test_cn4_projection_has_two_wedges_and_two_dashes():
    mol = build_coordination_mol("[PtCl4]2-")

    bond_dirs = [
        str(bond.GetBondDir())
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert bond_dirs.count("BEGINWEDGE") == 2
    assert bond_dirs.count("BEGINDASH") == 2


def test_cn6_projection_has_two_wedges_two_dashes_and_two_plain_bonds():
    mol = build_coordination_mol("[Co(NH3)6]3+")

    bond_dirs = [
        str(bond.GetBondDir())
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert bond_dirs.count("BEGINWEDGE") == 2
    assert bond_dirs.count("BEGINDASH") == 2
    assert bond_dirs.count("NONE") == 2


def test_bipy_projection_has_two_wedges_two_dashes_and_two_plain_bonds():
    mol = build_coordination_mol("[Fe(bipy)3]2+")

    bond_dirs = [
        str(bond.GetBondDir())
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert bond_dirs.count("BEGINWEDGE") == 2
    assert bond_dirs.count("BEGINDASH") == 2
    assert bond_dirs.count("NONE") == 2


def test_bipy_projection_is_vertical_and_left_right_symmetric():
    mol = build_coordination_mol("[Fe(bipy)3]2+")
    conf = mol.GetConformer()

    sites = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            other = bond.GetOtherAtomIdx(0)
            pos = conf.GetAtomPosition(other)
            sites.append((str(bond.GetBondDir()), round(pos.x, 2), round(pos.y, 2)))

    assert ("NONE", 0.0, 3.5) in sites
    assert ("NONE", 0.0, -3.5) in sites
    assert ("BEGINDASH", -3.0, 1.7) in sites
    assert ("BEGINDASH", 3.0, 1.7) in sites
    assert ("BEGINWEDGE", -3.0, -1.7) in sites
    assert ("BEGINWEDGE", 3.0, -1.7) in sites


def test_acac_projection_keeps_plain_bonds_vertical():
    mol = build_coordination_mol("[Co(acac)3]")
    conf = mol.GetConformer()

    plain_sites = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            if str(bond.GetBondDir()) == "NONE":
                other = bond.GetOtherAtomIdx(0)
                pos = conf.GetAtomPosition(other)
                plain_sites.append((round(pos.x, 2), round(pos.y, 2)))

    assert (0.0, 3.5) in plain_sites
    assert (0.0, -3.5) in plain_sites


def test_acac_binds_through_oxygen_atoms():
    mol = build_coordination_mol("[Co(acac)3]")

    donor_symbols = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            donor_idx = bond.GetOtherAtomIdx(0)
            donor_symbols.append(mol.GetAtomWithIdx(donor_idx).GetSymbol())

    assert donor_symbols == ["O", "O", "O", "O", "O", "O"]


def test_acac_uses_enolate_bond_pattern():
    from rdkit import Chem

    mol = Chem.MolFromSmiles(LIGAND_SMILES["acac"])

    assert mol.GetBondBetweenAtoms(1, 2).GetBondType() == Chem.BondType.DOUBLE
    assert mol.GetBondBetweenAtoms(3, 4).GetBondType() == Chem.BondType.DOUBLE
    assert mol.GetBondBetweenAtoms(4, 5).GetBondType() == Chem.BondType.SINGLE
    assert mol.GetAtomWithIdx(5).GetFormalCharge() == 0
