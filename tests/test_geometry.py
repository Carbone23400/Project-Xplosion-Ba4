"""
tests/test_geometry.py
----------------------
Test suite for coordchem.geometry

Run with: python -m pytest tests/test_geometry.py -v
"""

import pytest

from src.coordchem.parser import parse_formula
from src.coordchem.geometry import get_geometry, predict_geometry, geometry_report, get_d_count

# ===========================================================================
# Helper
# ===========================================================================

def geom(formula: str) -> str:
    """Shorthand."""
    return get_geometry(formula)

# ===========================================================================
# 1. Basic geometry from formula
# ===========================================================================

class TestBasicGeometry:
    def test_linear_ag(self):
        assert geom("[Ag(NH3)2]+") == "linear"
    
    def test_octahedral_fe_cyanide(self):
        assert geom("[Fe(CN)6]4-") == "octahedral"
    
    def test_octahedral_co_en(self):
        assert geom("[Co(en)3]3+") == "octahedral"
    
    def test_octahedral_cr_carbonyl(self):
        assert geom("[Cr(CO)6]") == "octahedral"
    
    def test_square_planar_pt(self):
        assert geom("[Pt(NH3)2Cl2]") == "square planar"
    
    def test_square_planar_pd(self):
        assert geom("[PdCl4]2-") == "square planar"
    
    def test_ambiguous_ni(self):
        assert geom("[NiCl4]2-") == "tetrahedral or square planar"
    
    def test_ambiguous_cu(self):
        assert geom("[Cu(NH3)4]2+") == "distorted square planar or tetrahedral"

# ===========================================================================
# 2. Geometry from ParsedComplex object
# ===========================================================================

class TestParsedComplexInput:
    
    def test_predict_geometry_from_parsed(self):
        parsed = parse_formula("[Fe(CN)6]4-")
        assert predict_geometry(parsed) == "octahedral"
    
    def test_get_geometry_from_parsed(self):
        parsed = parse_formula("[Pt(NH3)2Cl2]")
        assert get_geometry(parsed) == "square planar"


# ===========================================================================
# 2b. Geometry from compound name
# ===========================================================================

class TestCompoundNameInput:

    def test_geometry_from_full_name(self):
        assert get_geometry("tetraamminecopper(II)") == "distorted square planar or tetrahedral"

    def test_d_count_from_full_name(self):
        assert get_d_count("diamminedichloridoplatinum(II)") == 8

    def test_report_from_full_name(self):
        report = geometry_report("hexacyanidoferrate(II)")
        assert report["metal"] == "Fe"
        assert report["coordination_number"] == 6
        assert report["geometry"] == "octahedral"

# ===========================================================================
# 3. Coordination number logic
# ===========================================================================

class TestCoordinationNumberRules:

    def test_cn2_linear(self):
        assert geom("[Ag(NH3)2]+") == "linear"
    
    def test_cn4_pt_square_planar(self):
        assert geom("[Pt(NH3)2Cl2]") == "square planar"
    
    def test_cn5(self):
        assert geom("[Fe(CO)5]") == "trigonal bipyramidal or square pyramidal"
    
    def test_cn6_octahedral(self):
        assert geom("[Fe(H2O)6]3+") == "octahedral"
    
    def test_cn7(self):
        assert geom("[Mo(CN)7]4-") == "pentagonal bipyramidal or capped octahedral"
    
    def test_cn8(self):
        assert geom("[Zn(ox)4]4-") == "square antiprismatic or dodecahedral"

# ===========================================================================
# 4. CN = 4 special cases 
# ===========================================================================

class TestCN4SpecialCases:

    def test_pt_ii_square_planar(self):
        assert geom("[PtCl4]2-") == "square planar"

    def test_pd_ii_square_planar(self):
        assert geom("[PdCl4]2-") == "square planar"

    def test_ni_ii_strong_field_square_planar(self):
        assert geom("[Ni(CN)4]2-") == "square planar"

    def test_ni_ii_weak_field_ambiguous(self):
        assert geom("[NiCl4]2-") == "tetrahedral or square planar"

    def test_cu_ii_distorted(self):
        assert geom("[CuCl4]2-") == "distorted square planar or tetrahedral"

    def test_zn_ii_tetrahedral(self):
        assert geom("[ZnCl4]2-") == "tetrahedral"

    def test_cd_ii_tetrahedral(self):
        assert geom("[CdCl4]2-") == "tetrahedral"


# ===========================================================================
# 5. Geometry report
# ===========================================================================

class TestGeometryReport:

    def test_report_is_dict(self):
        report = geometry_report("[Fe(CN)6]4-")
        assert isinstance(report, dict)

    def test_report_contains_geometry(self):
        report = geometry_report("[Fe(CN)6]4-")
        assert report["geometry"] == "octahedral"

    def test_report_contains_metal(self):
        report = geometry_report("[PtCl2(NH3)2]")
        assert report["metal"] == "Pt"

    def test_report_contains_oxidation_state(self):
        report = geometry_report("[PtCl2(NH3)2]")
        assert report["oxidation_state"] == 2

    def test_report_contains_coordination_number(self):
        report = geometry_report("[Co(en)3]3+")
        assert report["coordination_number"] == 6

    def test_report_contains_ligands(self):
        report = geometry_report("[Co(en)3]3+")
        assert report["ligands"] == {"en": 3}


# ===========================================================================
# 6. Error handling
# ===========================================================================

class TestErrorHandling:

    def test_invalid_input_type(self):
        with pytest.raises(TypeError):
            get_geometry(123)

    def test_invalid_formula_raises(self):
        with pytest.raises(Exception):
            get_geometry("[Xx(CN)6]4-")

# ===========================================================================
# 7. d-electron count
# ===========================================================================

class TestDCount:

    def test_fe_ii_d6(self):
        assert get_d_count("[Fe(CN)6]4-") == 6

    def test_fe_iii_d5(self):
        assert get_d_count("[Fe(CN)6]3-") == 5

    def test_cu_ii_d9(self):
        assert get_d_count("[Cu(NH3)4]2+") == 9

    def test_pt_ii_d8(self):
        assert get_d_count("[PtCl2(NH3)2]") == 8

    def test_cr_zero_d6(self):
        assert get_d_count("[Cr(CO)6]") == 6

    def test_d_count_from_parsed_complex(self):
        parsed = parse_formula("[Co(en)3]3+")
        assert get_d_count(parsed) == 6

    def test_report_contains_d_count(self):
        report = geometry_report("[PtCl2(NH3)2]")
        assert report["d_count"] == 8

# ===========================================================================
# 8. Real-world benchmark complexes
# ===========================================================================

class TestBenchmarkComplexes:

    def test_hexacyanoferrate_II(self):
        assert geom("[Fe(CN)6]4-") == "octahedral"

    def test_hexacyanoferrate_III(self):
        assert geom("[Fe(CN)6]3-") == "octahedral"

    def test_cisplatin(self):
        assert geom("[PtCl2(NH3)2]") == "square planar"

    def test_tetraamminecopper(self):
        assert geom("[Cu(NH3)4]2+") == "distorted square planar or tetrahedral"

    def test_tris_en_cobalt(self):
        assert geom("[Co(en)3]3+") == "octahedral"

    def test_hexaaquairon_III(self):
        assert geom("[Fe(H2O)6]3+") == "octahedral"

    def test_chromium_carbonyl(self):
        assert geom("[Cr(CO)6]") == "octahedral"

    def test_tetrachloroaurate(self):
        assert geom("[AuCl4]-") == "tetrahedral or square planar"

    def test_diamminesilver(self):
        assert geom("[Ag(NH3)2]+") == "linear"
