"""
tests/test_name2.py
-------------------
Test suite for coordchem.name2.

Run with: python -m pytest tests/test_name2.py -v
"""

import pytest

from coordchem.name2 import (
    extract_complex_charge_from_name,
    ligand_data,
    metal_data,
    parse_name,
)
from src.coordchem.parser import ParsedComplex


# ===========================================================================
# 1. Oxidation state from Roman numerals
# ===========================================================================

class TestOxidationStateFromName:

    def test_extract_iron_ii(self):
        assert extract_complex_charge_from_name("hexacyanidoferrate(II)") == 2

    def test_extract_cobalt_iii(self):
        assert extract_complex_charge_from_name("hexaamminecobalt(III)") == 3

    def test_extract_none_without_roman_numeral(self):
        assert extract_complex_charge_from_name("chromium carbonyl") is None


# ===========================================================================
# 2. Metal name parsing
# ===========================================================================

class TestMetalData:

    def test_cation_metal_name(self):
        assert metal_data("tetraamminecopper(II)") == "Cu"

    def test_anion_metal_name(self):
        assert metal_data("hexacyanidoferrate(II)") == "Fe"

    def test_platinum_is_not_detected_as_tin(self):
        assert metal_data("diamminedichloridoplatinum(II)") == "Pt"

    def test_unknown_metal_raises(self):
        with pytest.raises(ValueError):
            metal_data("tetraammineunknownium(II)")


# ===========================================================================
# 3. Ligand parsing from names
# ===========================================================================

class TestLigandData:

    def test_tetraammine(self):
        assert ligand_data("tetraamminecopper(II)") == {"NH3": 4}

    def test_hexacyanido(self):
        assert ligand_data("hexacyanidoferrate(II)") == {"CN": 6}

    def test_mixed_ligands(self):
        assert ligand_data("diamminedichloridoplatinum(II)") == {"Cl": 2, "NH3": 2}

    def test_name_normalization(self):
        assert ligand_data("di ammine-di chlorido platinum(II)") == {"Cl": 2, "NH3": 2}


# ===========================================================================
# 4. Full name parsing
# ===========================================================================

class TestParseName:

    def test_parse_name_returns_parsed_complex(self):
        result = parse_name("tetraamminecopper(II)")
        assert isinstance(result, ParsedComplex)

    def test_parse_tetraamminecopper_ii(self):
        result = parse_name("tetraamminecopper(II)")
        assert result.metal == "Cu"
        assert result.ligands == {"NH3": 4}
        assert result.oxidation_state == 2
        assert result.complex_charge == 2
        assert result.coordination_number == 4

    def test_parse_diamminedichloridoplatinum_ii(self):
        result = parse_name("diamminedichloridoplatinum(II)")
        assert result.metal == "Pt"
        assert result.ligands == {"Cl": 2, "NH3": 2}
        assert result.oxidation_state == 2
        assert result.complex_charge == 0
        assert result.coordination_number == 4

    def test_parse_hexacyanidoferrate_ii(self):
        result = parse_name("hexacyanidoferrate(II)")
        assert result.metal == "Fe"
        assert result.ligands == {"CN": 6}
        assert result.oxidation_state == 2
        assert result.complex_charge == -4
        assert result.coordination_number == 6
