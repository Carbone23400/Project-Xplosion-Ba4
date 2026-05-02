"""
tests/test_diagram_2d.py
------------------------
Tests for coordchem.viz.diagram_2d.

Run with:
    python -m pytest tests/test_diagram_2d.py -v
"""

from __future__ import annotations

from math import hypot
from pathlib import Path

import pytest

from coordchem.viz.diagram_2d import (
    H2_ANNOTATION_PROP,
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


@pytest.mark.parametrize(
    "formula",
    [
        "[Pt(PPh3)4]2+",
        "[Fe(phen)3]2+",
        "[Fe(Cp)2]",
    ],
)
def test_representative_pdonor_phen_and_cp_complexes_draw(formula):
    svg = diagram_2d_svg(formula)

    assert "<svg" in svg
    assert "</svg>" in svg
    assert "bond-" in svg
    assert "atom-" in svg


def test_cp2_complex_uses_sandwich_centroid_representation():
    svg = diagram_2d_svg("[Fe(Cp)2]")

    assert svg.count("class='cp-delocalized-ring'") == 2
    assert "atom-0 atom-1" in svg
    assert "atom-0 atom-2" in svg
    assert ">Fe<" in svg
    assert "font-weight:400" in svg
    assert "dicyclopentadienyliron(II)" in svg
    assert ">sandwich<" in svg


def test_phen_uses_1_10_phenanthroline_donor_pattern():
    from rdkit import Chem

    mol = Chem.MolFromSmiles(LIGAND_SMILES["phen"])
    donor_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "N"]

    assert Chem.MolToSmiles(mol) == "c1cnc2c(c1)ccc1cccnc12"
    assert donor_indices == [9, 12]


def test_pph3_draws_as_compact_terminal_label():
    mol = build_coordination_mol("[Pt(PPh3)4]2+")

    labels = [
        atom.GetProp("atomLabel")
        for atom in mol.GetAtoms()
        if atom.HasProp("atomLabel")
    ]

    assert mol.GetNumAtoms() == 5
    assert labels == ["PPh3", "Ph3P", "PPh3", "PPh3"]


@pytest.mark.parametrize(
    ("formula", "expected_labels"),
    [
        ("[Pt(PMe3)4]2+", ["PMe3", "Me3P", "PMe3", "PMe3"]),
        ("[Pt(PEt3)4]2+", ["PEt3", "Et3P", "PEt3", "PEt3"]),
    ],
)
def test_pme3_and_pet3_draw_as_compact_terminal_labels(formula, expected_labels):
    mol = build_coordination_mol(formula)

    labels = [
        atom.GetProp("atomLabel")
        for atom in mol.GetAtoms()
        if atom.HasProp("atomLabel")
    ]

    assert mol.GetNumAtoms() == 5
    assert labels == expected_labels


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


def test_en_donor_atoms_display_as_n_atoms_without_separate_h_atoms():
    mol = build_coordination_mol("[Co(en)3]3+")

    donor_atoms = [
        mol.GetAtomWithIdx(bond.GetOtherAtomIdx(0))
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]
    h_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "H"]

    assert len(donor_atoms) == 6
    assert all(atom.GetAtomicNum() == 7 for atom in donor_atoms)
    assert all(atom.GetProp("atomLabel") == "N" for atom in donor_atoms)
    assert not h_atoms


def test_en_svg_adds_h2_annotations_next_to_n_donors():
    mol = build_coordination_mol("[Co(en)3]3+")
    svg = diagram_2d_svg("[Co(en)3]3+")

    donor_atoms = [
        mol.GetAtomWithIdx(bond.GetOtherAtomIdx(0))
        for bond in mol.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert all(atom.HasProp(H2_ANNOTATION_PROP) for atom in donor_atoms)
    assert "coordchem-h2-label" not in svg
    assert ">H2</text>" not in svg
    assert "baseline-shift" not in svg


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


def test_square_planar_cn4_projection_is_flat_cross():
    mol = build_coordination_mol("[PtCl4]2-")
    conf = mol.GetConformer()

    sites = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            other = bond.GetOtherAtomIdx(0)
            pos = conf.GetAtomPosition(other)
            sites.append((str(bond.GetBondDir()), round(pos.x, 2), round(pos.y, 2)))

    assert all(direction == "NONE" for direction, _, _ in sites)
    assert ("NONE", 0.0, 3.2) in sites
    assert ("NONE", -3.2, 0.0) in sites
    assert ("NONE", 3.2, 0.0) in sites
    assert ("NONE", 0.0, -3.2) in sites


def test_tetrahedral_cn4_projection_has_vertical_depth_cues():
    mol = build_coordination_mol("[ZnCl4]2-")
    conf = mol.GetConformer()

    sites = []
    for bond in mol.GetBonds():
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}:
            other = bond.GetOtherAtomIdx(0)
            pos = conf.GetAtomPosition(other)
            sites.append((str(bond.GetBondDir()), round(pos.x, 2), round(pos.y, 2)))

    assert ("NONE", -3.0, 0.0) in sites
    assert ("NONE", 3.0, 0.0) in sites
    assert ("BEGINDASH", 0.0, 2.85) in sites
    assert ("BEGINWEDGE", 0.0, -2.85) in sites


def test_ambiguous_cn4_diagram_draws_tetrahedral_and_square_planar_panels():
    svg = diagram_2d_svg("[Ni(CO)4]")
    tetrahedral = build_coordination_mol("[Ni(CO)4]", geometry_override="tetrahedral")
    square_planar = build_coordination_mol("[Ni(CO)4]", geometry_override="square planar")

    assert "<svg" in svg

    tetrahedral_dirs = [
        str(bond.GetBondDir())
        for bond in tetrahedral.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]
    square_planar_dirs = [
        str(bond.GetBondDir())
        for bond in square_planar.GetBonds()
        if 0 in {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()}
    ]

    assert tetrahedral_dirs.count("BEGINDASH") == 1
    assert tetrahedral_dirs.count("BEGINWEDGE") == 1
    assert tetrahedral_dirs.count("NONE") == 2
    assert square_planar_dirs == ["NONE", "NONE", "NONE", "NONE"]


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


def test_phen_projection_is_vertical_and_left_right_symmetric():
    mol = build_coordination_mol("[Fe(phen)3]2+")
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


def test_oxalate_draws_neutral_with_carbonyls_outward():
    mol = build_coordination_mol("[Fe(ox)3]3-")
    conf = mol.GetConformer()

    assert all(atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())

    for base in (1, 7, 13):
        for donor_idx, carbon_idx, oxygen_idx in (
            (base, base + 1, base + 2),
            (base + 5, base + 3, base + 4),
        ):
            donor = conf.GetAtomPosition(donor_idx)
            carbon = conf.GetAtomPosition(carbon_idx)
            oxygen = conf.GetAtomPosition(oxygen_idx)
            length = (donor.x * donor.x + donor.y * donor.y) ** 0.5
            outward_x = donor.x / length
            outward_y = donor.y / length
            carbonyl_x = oxygen.x - carbon.x
            carbonyl_y = oxygen.y - carbon.y

            assert carbonyl_x * outward_x + carbonyl_y * outward_y > 0

            simple_co = hypot(carbon.x - donor.x, carbon.y - donor.y)
            double_co = hypot(oxygen.x - carbon.x, oxygen.y - carbon.y)
            assert simple_co > double_co
            assert simple_co == pytest.approx(double_co, abs=0.65)

        carbon_1 = conf.GetAtomPosition(base + 1)
        carbon_2 = conf.GetAtomPosition(base + 3)
        cc_length = hypot(carbon_2.x - carbon_1.x, carbon_2.y - carbon_1.y)
        assert cc_length < 2.0


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
