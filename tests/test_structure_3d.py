"""
tests/test_structure_3d.py
--------------------------
Tests for ``coordchem.viz.molecule3D``.

Run with:
    python -m pytest tests/test_structure_3d.py -v
"""

import pytest

pytest.importorskip("rdkit")

from rdkit import Chem  # noqa: E402

from coordchem.complex import Complex  # noqa: E402
from coordchem.parser import parse_formula  # noqa: E402
from coordchem.viz.molecule3D import (  # noqa: E402
    build_complex_3d,
    build_ligand_3d,
    find_donor_atom,
    geometry_positions,
    octahedral_positions,
    parse_complex_input,
)


class TestGeometryPositions:
    def test_octahedral_returns_six_axes(self):
        pos = octahedral_positions(distance=2.0)
        assert len(pos) == 6
        # Sites are 2 Å away from origin along the axes
        for x, y, z in pos:
            assert pytest_approx(x ** 2 + y ** 2 + z ** 2) == 4.0

    def test_geometry_positions_octahedral(self):
        pos = geometry_positions("octahedral", 6)
        assert len(pos) == 6

    def test_geometry_positions_tetrahedral(self):
        pos = geometry_positions("tetrahedral", 4)
        assert len(pos) == 4

    def test_geometry_positions_unknown_falls_back(self):
        # Unknown geometry but n=6 should still give 6 sites
        pos = geometry_positions("unknown geometry", 6)
        assert len(pos) == 6

    def test_geometry_positions_square_antiprismatic(self):
        pos = geometry_positions("square antiprismatic or dodecahedral", 8)

        assert len(pos) == 8
        assert all(any(abs(component) > 0 for component in site) for site in pos)


class TestBuildLigand:
    def test_build_water(self):
        mol = build_ligand_3d("O")
        assert mol.GetNumConformers() == 1
        # H2O has 1 O + 2 H after AddHs
        assert mol.GetNumAtoms() == 3

    def test_build_invalid_smiles(self):
        with pytest.raises(ValueError):
            build_ligand_3d("not_a_smiles_string!!")

    def test_find_donor_in_ammonia(self):
        mol = build_ligand_3d("N")
        idx = find_donor_atom(mol, "N")
        assert mol.GetAtomWithIdx(idx).GetSymbol() == "N"

    def test_find_donor_with_override(self):
        mol = build_ligand_3d("NCCN")  # ethylenediamine
        idx = find_donor_atom(mol, "N", override=0)
        assert idx == 0
        assert mol.GetAtomWithIdx(idx).GetSymbol() == "N"

    def test_find_donor_missing_raises(self):
        mol = build_ligand_3d("N")
        with pytest.raises(ValueError):
            find_donor_atom(mol, "P")


class TestBuildComplex:
    def test_parse_complex_input_accepts_methyl_formula(self):
        parsed = parse_complex_input("[Ti(CH3)4]")

        assert parsed.metal == "Ti"
        assert parsed.ligands == {"CH3": 4}

    def test_hexacyanoferrate_builds(self):
        parsed = parse_formula("[Fe(CN)6]4-")
        mol = build_complex_3d(parsed)

        assert isinstance(mol, Chem.Mol)
        assert mol.GetNumConformers() == 1

        # Metal at index 0 should be Fe at the origin
        assert mol.GetAtomWithIdx(0).GetSymbol() == "Fe"
        conf = mol.GetConformer()
        origin = conf.GetAtomPosition(0)
        assert pytest_approx(origin.x ** 2 + origin.y ** 2 + origin.z ** 2) == 0.0

        # 6 dative bonds from the metal
        metal_atom = mol.GetAtomWithIdx(0)
        dative_bonds = [
            b for b in metal_atom.GetBonds()
            if b.GetBondType() == Chem.BondType.DATIVE
        ]
        assert len(dative_bonds) == 6

    def test_complex_class_draw_3d_html(self):
        py3Dmol = pytest.importorskip("py3Dmol")  # noqa: F841

        c = Complex.from_formula("[Fe(CN)6]4-")
        html = c.draw_3d_html(width=300, height=300)

        assert isinstance(html, str)
        assert len(html) > 0

    def test_complex_class_build_3d(self):
        c = Complex.from_formula("[Fe(CN)6]4-")
        mol = c.build_3d()

        assert mol.GetNumConformers() >= 1
        assert mol.GetAtomWithIdx(0).GetSymbol() == "Fe"

    def test_tetrahedral_complex_has_four_sites(self):
        parsed = parse_formula("[Zn(NH3)4]2+")
        mol = build_complex_3d(parsed)

        metal_atom = mol.GetAtomWithIdx(0)
        dative_bonds = [
            b for b in metal_atom.GetBonds()
            if b.GetBondType() == Chem.BondType.DATIVE
        ]
        assert len(dative_bonds) == 4

    @pytest.mark.parametrize(
        ("formula", "expected_donor_symbols"),
        [
            ("[TaF8]3-", ["F"] * 8),
            ("[Zn(ox)4]6-", ["O"] * 8),
            ("[Zr(ox)2F4]4-", ["O"] * 4 + ["F"] * 4),
        ],
    )
    def test_cn8_3d_handles_monodentate_bidentate_and_mixed_ligands(
        self,
        formula,
        expected_donor_symbols,
    ):
        parsed = parse_formula(formula)
        mol = build_complex_3d(parsed)

        assert parsed.coordination_number == 8
        assert _metal_donor_symbols(mol) == expected_donor_symbols

    def test_methyl_complex_builds_four_ch3_ligands(self):
        parsed = parse_formula("[Ti(CH3)4]")
        mol = build_complex_3d(parsed)

        donor_symbols = _metal_donor_symbols(mol)
        methyl_carbons = [
            atom
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "C"
        ]

        assert parsed.ligands == {"CH3": 4}
        assert donor_symbols == ["C"] * 4
        assert len(methyl_carbons) == 4
        assert all(
            sum(
                1 for neighbor in atom.GetNeighbors()
                if neighbor.GetSymbol() == "H"
            ) == 3
            for atom in methyl_carbons
        )

    def test_dmso_hard_metal_uses_oxygen_donor(self):
        parsed = parse_formula("[Fe(dmso)6]3+")
        mol = build_complex_3d(parsed)

        donor_symbols = _metal_donor_symbols(mol)
        assert donor_symbols == ["O"] * 6

    def test_dmso_soft_metal_uses_sulfur_donor(self):
        parsed = parse_formula("[Pt(dmso)4]2+")
        mol = build_complex_3d(parsed)

        donor_symbols = _metal_donor_symbols(mol)
        assert donor_symbols == ["S"] * 4

    def test_dmso_3d_sulfur_donor_has_no_hydrogen_neighbor(self):
        parsed = parse_formula("[Pt(dmso)4]2+")
        mol = build_complex_3d(parsed)

        sulfur_donors = [
            mol.GetAtomWithIdx(bond.GetOtherAtomIdx(0))
            for bond in mol.GetAtomWithIdx(0).GetBonds()
            if mol.GetAtomWithIdx(bond.GetOtherAtomIdx(0)).GetSymbol() == "S"
        ]

        assert len(sulfur_donors) == 4
        assert all(
            all(neighbor.GetSymbol() != "H" for neighbor in atom.GetNeighbors())
            for atom in sulfur_donors
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _metal_donor_symbols(mol):
    metal_atom = mol.GetAtomWithIdx(0)
    return [
        mol.GetAtomWithIdx(bond.GetOtherAtomIdx(0)).GetSymbol()
        for bond in metal_atom.GetBonds()
        if bond.GetBondType() == Chem.BondType.DATIVE
    ]


def pytest_approx(value, rel=1e-6, abs_=1e-6):
    return pytest.approx(value, rel=rel, abs=abs_)
