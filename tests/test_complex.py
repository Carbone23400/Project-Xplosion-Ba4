"""
tests/test_complex.py
---------------------
Tests for coordchem.complex

Run with:
    python -m pytest tests/test_complex.py -v
"""

from coordchem.complex import Complex
from coordchem.parser import ParsedComplex


class TestComplexFromFormula:
    def test_from_formula_fe_cyanide(self):
        c = Complex.from_formula("[Fe(CN)6]4-")

        assert isinstance(c.parsed, ParsedComplex)
        assert c.metal == "Fe"
        assert c.ligands == {"CN": 6}
        assert c.oxidation_state == 2
        assert c.coordination_number == 6
        assert c.geometry == "octahedral"
        assert c.d_count == 6
        assert c.complex_charge == -4

    def test_from_formula_copper_ammine(self):
        c = Complex.from_formula("[Cu(NH3)4]2+")

        assert c.metal == "Cu"
        assert c.ligands == {"NH3": 4}
        assert c.oxidation_state == 2
        assert c.coordination_number == 4
        assert c.geometry == "distorted square planar or tetrahedral"
        assert c.d_count == 9
        assert c.complex_charge == 2

    def test_from_formula_unbracketed(self):
        c = Complex.from_input("PtCl2(NH3)2")

        assert c.metal == "Pt"
        assert c.ligands == {"Cl": 2, "NH3": 2}
        assert c.coordination_number == 4


class TestComplexFromName:
    def test_from_name_tetraamminecopper(self):
        c = Complex.from_name("tetraamminecopper(II)")

        assert isinstance(c.parsed, ParsedComplex)
        assert c.metal == "Cu"
        assert c.ligands == {"NH3": 4}
        assert c.oxidation_state == 2
        assert c.coordination_number == 4
        assert c.complex_charge == 2

    def test_from_input_name_fallback(self):
        c = Complex.from_input("hexaamminecobalt(III)")

        assert c.metal == "Co"
        assert c.ligands == {"NH3": 6}
        assert c.oxidation_state == 3
        assert c.coordination_number == 6
        assert c.geometry == "octahedral"


class TestComplexProperties:
    def test_report_keys(self):
        c = Complex.from_formula("[Fe(CN)6]4-")
        report = c.report()

        assert report["metal"] == "Fe"
        assert report["oxidation_state"] == 2
        assert report["coordination_number"] == 6
        assert report["geometry"] == "octahedral"
        assert report["d_count"] == 6
        assert report["ligands"] == {"CN": 6}

    def test_str_contains_main_information(self):
        c = Complex.from_formula("[Fe(CN)6]4-")
        text = str(c)

        assert "Fe" in text
        assert "CN" in text
        assert "Coordination no." in text
        assert "Oxidation state" in text

    def test_ligand_metadata(self):
        c = Complex.from_formula("[Co(en)3]3+")

        assert c.ligand_denticity["en"] == 2
        assert c.donor_atoms["en"] == "N"
        assert c.ligand_names["en"] == "ethylenediamine"
        assert c.coordination_number == 6
        assert c.geometry == "octahedral"


class TestComplexInputValidation:
    def test_from_input_rejects_empty_string(self):
        try:
            Complex.from_input("   ")
            assert False, "Expected ValueError for empty input"
        except ValueError:
            assert True

    def test_from_input_rejects_non_string(self):
        try:
            Complex.from_input(123)  # type: ignore[arg-type]
            assert False, "Expected TypeError for non-string input"
        except TypeError:
            assert True