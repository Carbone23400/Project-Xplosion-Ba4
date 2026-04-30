"""
tests/test_name.py
-------------------
Test suite for name parsing (coordchem.name)

Run with: python -m pytest tests/test_name.py -v
"""

import pytest

from src.coordchem.name2 import parse_name, ligand_data, metal_data, extract_complex_charge_from_name


# ===========================================================================
# 1. Metal extraction
# ===========================================================================

class TestMetalExtraction:

    def test_iron(self):
        assert metal_data("hexacyanoferrate") == "Fe"

    def test_copper(self):
        assert metal_data("tetraamminecopper") == "Cu"

    def test_platinum(self):
        assert metal_data("diamminedichloroplatinum") == "Pt"


# ===========================================================================
# 2. Oxidation state extraction
# ===========================================================================

class TestOxidationState:

    def test_fe_ii(self):
        assert extract_complex_charge_from_name("hexacyanoferrate(II)") == 2

    def test_fe_iii(self):
        assert extract_complex_charge_from_name("hexacyanoferrate(III)") == 3

    def test_no_oxidation(self):
        assert extract_complex_charge_from_name("tetraamminecopper") is None


# ===========================================================================
# 3. Ligand parsing
# ===========================================================================

class TestLigandParsing:

    def test_single_ligand(self):
        lig = ligand_data("ammine")
        assert lig == {"NH3": 1}

    def test_multiple_ligands(self):
        lig = ligand_data("tetraammine")
        assert lig == {"NH3": 4}

    def test_mixed_ligands(self):
        lig = ligand_data("diamminedichloro")
        assert lig.get("NH3", 0) == 2
        assert lig.get("Cl", 0) == 2


# ===========================================================================
# 4. Full parsing
# ===========================================================================

class TestParseName:

    def test_simple_complex(self):
        parsed = parse_name("tetraamminecopper(II)")
        assert parsed.metal == "Cu"
        assert parsed.oxidation_state == 2
        assert parsed.ligands.get("NH3", 0) == 4

    def test_hexacyanoferrate_ii(self):
        parsed = parse_name("hexacyanoferrate(II)")
        assert parsed.metal == "Fe"
        assert parsed.oxidation_state == 2
        assert parsed.ligands.get("CN", 0) == 6

    def test_hexacyanoferrate_iii(self):
        parsed = parse_name("hexacyanoferrate(III)")
        assert parsed.oxidation_state == 3

    def test_platinum_complex(self):
        parsed = parse_name("diamminedichloroplatinum(II)")
        assert parsed.metal == "Pt"
        assert parsed.ligands.get("NH3", 0) == 2
        assert parsed.ligands.get("Cl", 0) == 2


# ===========================================================================
# 5. Charge consistency
# ===========================================================================

class TestChargeConsistency:

    def test_charge_calculation(self):
        parsed = parse_name("hexacyanoferrate(II)")
        assert parsed.complex_charge is not None


# ===========================================================================
# 6. Error handling
# ===========================================================================

class TestErrors:

    def test_invalid_input(self):
        with pytest.raises(Exception):
            parse_name(123)

    def test_unknown_name(self):
        parsed = parse_name("blabla")
        raise ValueError