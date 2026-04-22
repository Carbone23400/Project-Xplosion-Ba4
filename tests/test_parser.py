"""
tests/test_parser.py
--------------------
Test suite for Project-Xplosion-Ba4.parser

Run with:  pytest tests/test_parser.py -v
"""

import pytest
from coordchem.parser import parse_formula, FormulaParseError


# ===========================================================================
# Helper
# ===========================================================================

def parse(formula):
    """Shorthand."""
    return parse_formula(formula)


# ===========================================================================
# 1. Metal extraction
# ===========================================================================

class TestMetalExtraction:

    def test_two_letter_metal(self):
        r = parse("[Fe(CN)6]4-")
        assert r.metal == "Fe"

    def test_one_letter_metal(self):
        r = parse("[V(H2O)6]3+")
        assert r.metal == "V"

    def test_platinum(self):
        r = parse("[PtCl2(NH3)2]")
        assert r.metal == "Pt"

    def test_cobalt(self):
        r = parse("[Co(en)3]3+")
        assert r.metal == "Co"

    def test_copper(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.metal == "Cu"

    def test_chromium(self):
        r = parse("[Cr(H2O)6]3+")
        assert r.metal == "Cr"

    def test_nickel(self):
        r = parse("[NiCl4]2-")
        assert r.metal == "Ni"


# ===========================================================================
# 2. Complex charge parsing
# ===========================================================================

class TestChargeExtraction:

    def test_negative_charge_4(self):
        r = parse("[Fe(CN)6]4-")
        assert r.complex_charge == -4

    def test_positive_charge_2(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.complex_charge == 2

    def test_positive_charge_3(self):
        r = parse("[Co(en)3]3+")
        assert r.complex_charge == 3

    def test_neutral_complex(self):
        r = parse("[PtCl2(NH3)2]")
        assert r.complex_charge == 0

    def test_charge_1_minus(self):
        r = parse("[AuCl4]-")
        assert r.complex_charge == -1

    def test_charge_1_plus(self):
        r = parse("[Ag(NH3)2]+")
        assert r.complex_charge == 1


# ===========================================================================
# 3. Ligand parsing — monodentate
# ===========================================================================

class TestMonodentateLigands:

    def test_cyanide_6(self):
        r = parse("[Fe(CN)6]4-")
        assert r.ligands == {"CN": 6}

    def test_ammonia_4(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.ligands == {"NH3": 4}

    def test_ammonia_6(self):
        r = parse("[Co(NH3)6]3+")
        assert r.ligands == {"NH3": 6}

    def test_water_6(self):
        r = parse("[Fe(H2O)6]3+")
        assert r.ligands == {"H2O": 6}

    def test_mixed_chloro_ammine(self):
        r = parse("[PtCl2(NH3)2]")
        assert r.ligands.get("Cl") == 2
        assert r.ligands.get("NH3") == 2

    def test_chloride_4(self):
        r = parse("[NiCl4]2-")
        assert r.ligands == {"Cl": 4}

    def test_mixed_water_chloro(self):
        r = parse("[CrCl2(H2O)4]+")
        assert r.ligands.get("Cl") == 2
        assert r.ligands.get("H2O") == 4

    def test_carbonyl(self):
        r = parse("[Cr(CO)6]")
        assert r.ligands == {"CO": 6}
        assert r.complex_charge == 0


# ===========================================================================
# 4. Ligand parsing — polydentate
# ===========================================================================

class TestPolydentateLigands:

    def test_ethylenediamine(self):
        r = parse("[Co(en)3]3+")
        assert r.ligands == {"en": 3}

    def test_en_denticity(self):
        r = parse("[Co(en)3]3+")
        assert r.ligand_denticity["en"] == 2

    def test_en_coordination_number(self):
        r = parse("[Co(en)3]3+")
        assert r.coordination_number == 6   # 3 × bidentate

    def test_mixed_en_chloro(self):
        r = parse("[CoCl2(en)2]+")
        assert r.ligands.get("Cl") == 2
        assert r.ligands.get("en") == 2

    def test_mixed_en_cn(self):
        """Coordination number of [Co(en)2(CN)2]+ should be 6."""
        r = parse("[Co(en)2(CN)2]+")
        assert r.coordination_number == 6

    def test_oxalate(self):
        r = parse("[Fe(ox)3]3-")
        assert r.ligands == {"ox": 3}
        assert r.ligand_denticity["ox"] == 2


# ===========================================================================
# 5. Oxidation state calculation
# ===========================================================================

class TestOxidationState:

    def test_fe_hexacyano(self):
        """[Fe(CN)6]4-  → Fe(II),  6×CN(-1) = -6,  charge=-4 → OS = -4-(-6) = +2"""
        r = parse("[Fe(CN)6]4-")
        assert r.oxidation_state == 2

    def test_fe_hexacyano_III(self):
        """[Fe(CN)6]3-  → Fe(III)"""
        r = parse("[Fe(CN)6]3-")
        assert r.oxidation_state == 3

    def test_cu_tetraammine(self):
        """[Cu(NH3)4]2+  → Cu(II),  NH3 neutral"""
        r = parse("[Cu(NH3)4]2+")
        assert r.oxidation_state == 2

    def test_co_hexaammine(self):
        """[Co(NH3)6]3+  → Co(III)"""
        r = parse("[Co(NH3)6]3+")
        assert r.oxidation_state == 3

    def test_pt_cisplatin(self):
        """[PtCl2(NH3)2]  → Pt(II),  2×Cl(-1) + 2×NH3(0) = -2,  charge=0 → OS=+2"""
        r = parse("[PtCl2(NH3)2]")
        assert r.oxidation_state == 2

    def test_cr_carbonyl(self):
        """[Cr(CO)6]  → Cr(0),  6×CO neutral"""
        r = parse("[Cr(CO)6]")
        assert r.oxidation_state == 0

    def test_au_tetrachloro(self):
        """[AuCl4]-  → Au(III),  4×Cl(-1) = -4,  charge=-1 → OS=+3"""
        r = parse("[AuCl4]-")
        assert r.oxidation_state == 3

    def test_co_en(self):
        """[Co(en)3]3+  → Co(III),  en neutral"""
        r = parse("[Co(en)3]3+")
        assert r.oxidation_state == 3

    def test_fe_hexaaqua(self):
        """[Fe(H2O)6]3+  → Fe(III)"""
        r = parse("[Fe(H2O)6]3+")
        assert r.oxidation_state == 3


# ===========================================================================
# 6. Coordination number
# ===========================================================================

class TestCoordinationNumber:

    def test_octahedral_mono(self):
        assert parse("[Fe(CN)6]4-").coordination_number == 6

    def test_square_planar(self):
        assert parse("[PtCl2(NH3)2]").coordination_number == 4

    def test_tetrahedral(self):
        assert parse("[NiCl4]2-").coordination_number == 4

    def test_linear(self):
        assert parse("[Ag(NH3)2]+").coordination_number == 2

    def test_bidentate_3(self):
        assert parse("[Co(en)3]3+").coordination_number == 6

    def test_bidentate_mixed(self):
        assert parse("[CoCl2(en)2]+").coordination_number == 6


# ===========================================================================
# 7. Counter ions
# ===========================================================================

class TestCounterIons:

    def test_potassium_counter_ion(self):
        r = parse("K4[Fe(CN)6]")
        assert r.counter_ions.get("K") == 4

    def test_sodium_counter_ion(self):
        r = parse("Na2[PtCl4]")
        assert r.counter_ions.get("Na") == 2

    def test_no_counter_ion(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.counter_ions == {}


# ===========================================================================
# 8. Ligand metadata
# ===========================================================================

class TestLigandMetadata:

    def test_cyanide_charge(self):
        r = parse("[Fe(CN)6]4-")
        assert r.ligand_charges["CN"] == -1

    def test_ammonia_neutral(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.ligand_charges["NH3"] == 0

    def test_chloride_charge(self):
        r = parse("[PtCl2(NH3)2]")
        assert r.ligand_charges["Cl"] == -1

    def test_en_name(self):
        r = parse("[Co(en)3]3+")
        assert "ethylenediamine" in r.ligand_names["en"]

    def test_cyanide_donor_atom(self):
        r = parse("[Fe(CN)6]4-")
        assert r.donor_atoms["CN"] == "C"

    def test_ammonia_donor_atom(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.donor_atoms["NH3"] == "N"

    def test_water_donor_atom(self):
        r = parse("[Fe(H2O)6]3+")
        assert r.donor_atoms["H2O"] == "O"


# ===========================================================================
# 9. Edge cases and robustness
# ===========================================================================

class TestEdgeCases:

    def test_whitespace_stripped(self):
        r = parse("  [Fe(CN)6]4-  ")
        assert r.metal == "Fe"

    def test_neutral_no_charge_suffix(self):
        r = parse("[Cr(CO)6]")
        assert r.complex_charge == 0
        assert r.oxidation_state == 0

    def test_unknown_ligand_gives_warning(self):
        r = parse("[Fe(XYZ)3]3+")
        assert any("XYZ" in w for w in r.warnings)

    def test_single_ligand_no_count(self):
        """[AuCl4]- — count inferred as 4 from (Cl)4 or Cl4."""
        r = parse("[AuCl4]-")
        assert r.ligands.get("Cl") == 4

    def test_invalid_metal_raises(self):
        with pytest.raises(FormulaParseError):
            parse("[Xx(CN)6]4-")

    def test_str_representation(self):
        """__str__ should not raise."""
        r = parse("[Fe(CN)6]4-")
        output = str(r)
        assert "Fe" in output
        assert "CN" in output


# ===========================================================================
# 10. Real-world benchmark complexes
# ===========================================================================

class TestBenchmarkComplexes:
    """Validate against well-known coordination chemistry textbook examples."""

    def test_hexacyanoferrate_II(self):
        r = parse("[Fe(CN)6]4-")
        assert r.metal == "Fe"
        assert r.oxidation_state == 2
        assert r.coordination_number == 6
        assert r.ligands == {"CN": 6}

    def test_hexacyanoferrate_III(self):
        r = parse("[Fe(CN)6]3-")
        assert r.oxidation_state == 3

    def test_cisplatin(self):
        r = parse("[PtCl2(NH3)2]")
        assert r.metal == "Pt"
        assert r.oxidation_state == 2
        assert r.coordination_number == 4

    def test_tetraamminecopper(self):
        r = parse("[Cu(NH3)4]2+")
        assert r.metal == "Cu"
        assert r.oxidation_state == 2
        assert r.coordination_number == 4

    def test_tris_en_cobalt(self):
        r = parse("[Co(en)3]3+")
        assert r.metal == "Co"
        assert r.oxidation_state == 3
        assert r.coordination_number == 6

    def test_hexaaquairon_III(self):
        r = parse("[Fe(H2O)6]3+")
        assert r.metal == "Fe"
        assert r.oxidation_state == 3
        assert r.coordination_number == 6

    def test_chromium_carbonyl(self):
        r = parse("[Cr(CO)6]")
        assert r.metal == "Cr"
        assert r.oxidation_state == 0
        assert r.coordination_number == 6

    def test_tetrachloroaurate(self):
        r = parse("[AuCl4]-")
        assert r.metal == "Au"
        assert r.oxidation_state == 3
        assert r.coordination_number == 4

    def test_diamminesilver(self):
        r = parse("[Ag(NH3)2]+")
        assert r.metal == "Ag"
        assert r.oxidation_state == 1
        assert r.coordination_number == 2

    def test_prussian_blue_unit(self):
        """[Fe(CN)6]3- is the building block of Prussian Blue."""
        r = parse("[Fe(CN)6]3-")
        assert r.oxidation_state == 3
        assert r.coordination_number == 6
