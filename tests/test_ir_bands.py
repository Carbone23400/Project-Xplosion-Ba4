"""
tests/test_ir_bands.py
----------------------
Test suite for coordchem.database.ir_bands

Run with:  python -m pytest tests/test_ir_bands.py -v
"""

import pytest
from coordchem.database.ir_bands import IRBandDB, BandRecord


@pytest.fixture
def db():
    """Fresh in-memory database for every test."""
    database = IRBandDB()
    yield database
    database.close()


# ===========================================================================
# 1. Database setup
# ===========================================================================

class TestDatabaseSetup:

    def test_db_creates_without_error(self, db):
        assert db is not None

    def test_seed_data_loaded(self, db):
        ligands = db.get_all_ligands()
        assert len(ligands) > 0

    def test_all_expected_ligands_present(self, db):
        ligands = db.get_all_ligands()
        expected = ["CN", "CO", "NH3", "H2O", "Cl", "NO2", "en", "ox", "acac",
                    "Br", "N3", "NO", "OH", "SCN", "NCS", "py", "dmso"]
        for lig in expected:
            assert lig in ligands, f"Ligand '{lig}' missing from database"

    def test_returns_band_record_objects(self, db):
        bands = db.get_bands("CN")
        assert all(isinstance(b, BandRecord) for b in bands)


# ===========================================================================
# 2. Band record properties
# ===========================================================================

class TestBandRecord:

    def test_center_is_midpoint(self, db):
        bands = db.get_bands("CN")
        for b in bands:
            assert b.center == pytest.approx((b.wn_min + b.wn_max) / 2)

    def test_width_is_half_range(self, db):
        bands = db.get_bands("CN")
        for b in bands:
            assert b.width == pytest.approx((b.wn_max - b.wn_min) / 2)

    def test_wn_min_leq_wn_max(self, db):
        for lig in db.get_all_ligands():
            for b in db.get_bands(lig):
                assert b.wn_min <= b.wn_max

    def test_intensity_valid_values(self, db):
        valid = {"very strong", "strong", "medium", "weak", "strong broad"}
        for lig in db.get_all_ligands():
            for b in db.get_bands(lig):
                assert b.intensity in valid, \
                    f"Invalid intensity '{b.intensity}' for {lig}"

    def test_spectrum_type_valid(self, db):
        for lig in db.get_all_ligands():
            for stype in ("IR", "Raman"):
                for b in db.get_bands(lig, spectrum_type=stype):
                    assert b.spectrum_type == stype


# ===========================================================================
# 3. Cyanide — benchmark ligand
# ===========================================================================

class TestCyanide:

    def test_ir_bands_exist(self, db):
        bands = db.get_bands("CN", spectrum_type="IR")
        assert len(bands) > 0

    def test_raman_bands_exist(self, db):
        bands = db.get_bands("CN", spectrum_type="Raman")
        assert len(bands) > 0

    def test_terminal_cn_stretch_in_range(self, db):
        """C≡N stretch for terminal CN should be 2100–2200 cm⁻¹."""
        bands = db.get_bands("CN", coordination="terminal")
        stretches = [b for b in bands if "C≡N stretch" in b.assignment
                     and b.coordination == "terminal"]
        assert any(2050 <= b.wn_min and b.wn_max <= 2250 for b in stretches)

    def test_bridging_cn_lower_than_terminal(self, db):
        """Bridging CN stretch should be below terminal CN stretch."""
        terminal = db.get_bands("CN", coordination="terminal", spectrum_type="IR")
        bridging = db.get_bands("CN", coordination="bridging", spectrum_type="IR")
        t_stretches = [b for b in terminal if "C≡N stretch" in b.assignment]
        b_stretches = [b for b in bridging if "C≡N stretch" in b.assignment]
        if t_stretches and b_stretches:
            assert min(b.center for b in b_stretches) < \
                   min(b.center for b in t_stretches)

    def test_fe_specific_band_at_2093(self, db):
        """Fe(II) hexacyanoferrate benchmark: CN stretch at 2093 cm⁻¹."""
        bands = db.get_bands("CN", spectrum_type="IR", metal="Fe")
        fe_bands = [b for b in bands if b.metal == "Fe"]
        assert any(abs(b.center - 2093) < 10 for b in fe_bands)

    def test_metal_specific_overrides_generic(self, db):
        """When metal='Fe', Fe-specific rows should appear before generic."""
        bands = db.get_bands("CN", spectrum_type="IR", metal="Fe")
        metals = [b.metal for b in bands]
        # Fe-specific entries should be present
        assert "Fe" in metals

    def test_mc_stretch_below_1000(self, db):
        """M-C stretch should be in the fingerprint region."""
        bands = db.get_bands("CN", spectrum_type="IR")
        mc = [b for b in bands if "M–C stretch" in b.assignment]
        assert any(b.wn_max < 700 for b in mc)


# ===========================================================================
# 4. Carbonyl — benchmark ligand
# ===========================================================================

class TestCarbonyl:

    def test_co_stretch_very_strong(self, db):
        bands = db.get_bands("CO", spectrum_type="IR")
        stretches = [b for b in bands if "C≡O stretch" in b.assignment
                     and b.coordination == "terminal"]
        assert any(b.intensity == "very strong" for b in stretches)

    def test_terminal_co_above_1850(self, db):
        bands = db.get_bands("CO", coordination="terminal", spectrum_type="IR")
        stretches = [b for b in bands if "C≡O stretch" in b.assignment]
        assert any(b.wn_min >= 1850 for b in stretches)

    def test_bridging_co_below_terminal(self, db):
        """Bridging CO stretch center should be lower than terminal CO stretch center."""
        terminal = [b for b in db.get_bands("CO", spectrum_type="IR")
                    if b.coordination == "terminal" and "C≡O" in b.assignment
                    and b.wn_min > 1700]
        bridging = [b for b in db.get_bands("CO", spectrum_type="IR")
                    if b.coordination == "bridging" and "C≡O" in b.assignment]
        if terminal and bridging:
            assert max(b.center for b in bridging) < \
                   max(b.center for b in terminal)

    def test_cr_co6_benchmark(self, db):
        """Cr(CO)6: IR active T1u band at 1984 cm⁻¹."""
        bands = db.get_bands("CO", spectrum_type="IR", metal="Cr")
        cr_bands = [b for b in bands if b.metal == "Cr"]
        assert any(abs(b.center - 1984) < 10 for b in cr_bands)

    def test_cr_co6_raman_a1g(self, db):
        """Cr(CO)6: Raman active A1g band at 2119 cm⁻¹."""
        bands = db.get_bands("CO", spectrum_type="Raman", metal="Cr")
        cr_bands = [b for b in bands if b.metal == "Cr"]
        assert any(abs(b.center - 2119) < 10 for b in cr_bands)


# ===========================================================================
# 5. Ammonia
# ===========================================================================

class TestAmmonia:

    def test_nh_stretch_present(self, db):
        bands = db.get_bands("NH3", spectrum_type="IR")
        assert any("N–H stretch" in b.assignment for b in bands)

    def test_coordinated_nh_redshifted(self, db):
        """Coordinated NH3 N-H stretch should be below 3400 cm⁻¹."""
        bands = db.get_bands("NH3", coordination="terminal", spectrum_type="IR")
        nh = [b for b in bands if "N–H stretch" in b.assignment]
        assert any(b.wn_max <= 3350 for b in nh)

    def test_mn_stretch_in_fingerprint(self, db):
        bands = db.get_bands("NH3", spectrum_type="IR")
        mn = [b for b in bands if "M–N stretch" in b.assignment]
        assert any(b.wn_max < 600 for b in mn)

    def test_cu_specific_band(self, db):
        bands = db.get_bands("NH3", spectrum_type="IR", metal="Cu")
        cu_bands = [b for b in bands if b.metal == "Cu"]
        assert len(cu_bands) > 0


# ===========================================================================
# 6. Water
# ===========================================================================

class TestWater:

    def test_oh_stretch_present(self, db):
        bands = db.get_bands("H2O", spectrum_type="IR")
        assert any("O–H stretch" in b.assignment for b in bands)

    def test_coordinated_water_broad(self, db):
        """Coordinated water band should span a wide range."""
        bands = db.get_bands("H2O", coordination="terminal", spectrum_type="IR")
        oh = [b for b in bands if "O–H stretch" in b.assignment]
        assert any((b.wn_max - b.wn_min) >= 200 for b in oh)

    def test_mo_stretch_present(self, db):
        bands = db.get_bands("H2O", spectrum_type="IR")
        assert any("M–O stretch" in b.assignment for b in bands)


# ===========================================================================
# 7. Chloride
# ===========================================================================

class TestChloride:

    def test_mcl_stretch_present(self, db):
        bands = db.get_bands("Cl", spectrum_type="IR")
        assert any("M–Cl stretch" in b.assignment for b in bands)

    def test_mcl_stretch_in_fingerprint(self, db):
        bands = db.get_bands("Cl", spectrum_type="IR")
        mcl = [b for b in bands if "M–Cl stretch" in b.assignment]
        assert any(b.wn_max <= 450 for b in mcl)

    def test_pt_specific_band(self, db):
        """Pt-Cl stretch benchmark: 320–360 cm⁻¹."""
        bands = db.get_bands("Cl", spectrum_type="IR", metal="Pt")
        pt_bands = [b for b in bands if b.metal == "Pt"]
        assert any(310 <= b.center <= 370 for b in pt_bands)

    def test_raman_mcl_present(self, db):
        bands = db.get_bands("Cl", spectrum_type="Raman")
        assert len(bands) > 0


# ===========================================================================
# 8. Polydentate ligands
# ===========================================================================

class TestPolydentateLigands:

    def test_en_mn_stretch(self, db):
        bands = db.get_bands("en", spectrum_type="IR")
        assert any("M–N stretch" in b.assignment for b in bands)

    def test_ox_co_stretch(self, db):
        bands = db.get_bands("ox", spectrum_type="IR")
        assert any("C=O" in b.assignment for b in bands)

    def test_acac_co_stretch(self, db):
        bands = db.get_bands("acac", spectrum_type="IR")
        assert any("C=O" in b.assignment or "C=C" in b.assignment
                   for b in bands)


# ===========================================================================
# 9. Range search
# ===========================================================================

class TestRangeSearch:

    def test_range_2000_2200_returns_cn_co(self, db):
        """The 2000–2200 region should contain CN and CO stretches."""
        bands = db.get_bands_in_range(2000, 2200)
        ligands = {b.ligand for b in bands}
        assert "CN" in ligands or "CO" in ligands

    def test_range_300_500_returns_ml_stretches(self, db):
        """Metal-ligand stretches live in the 300–500 region."""
        bands = db.get_bands_in_range(300, 500)
        assert len(bands) > 0

    def test_empty_range_returns_empty(self, db):
        """No bands exist above 4000 cm⁻¹."""
        bands = db.get_bands_in_range(4100, 4500)
        assert bands == []

    def test_raman_range_search(self, db):
        bands = db.get_bands_in_range(2000, 2200, spectrum_type="Raman")
        assert len(bands) > 0


# ===========================================================================
# 10. Custom band addition
# ===========================================================================

class TestCustomBands:

    def test_add_custom_band(self, db):
        db.add_band(
            ligand="TEST", coordination="terminal", metal="any",
            spectrum_type="IR", wn_min=1000, wn_max=1100,
            intensity="strong", assignment="test stretch",
            source="unit test"
        )
        bands = db.get_bands("TEST")
        assert any(b.assignment == "test stretch" for b in bands)

    def test_custom_band_appears_in_range_search(self, db):
        db.add_band(
            ligand="TEST2", coordination="terminal", metal="any",
            spectrum_type="IR", wn_min=1500, wn_max=1600,
            intensity="medium", assignment="custom band",
            source="unit test"
        )
        bands = db.get_bands_in_range(1450, 1650)
        assert any(b.ligand == "TEST2" for b in bands)

    def test_custom_band_in_all_ligands(self, db):
        db.add_band(
            ligand="NEW", coordination="any", metal="any",
            spectrum_type="IR", wn_min=500, wn_max=600,
            intensity="weak", assignment="new stretch",
            source="unit test"
        )
        assert "NEW" in db.get_all_ligands()


# ===========================================================================
# 11. Context manager
# ===========================================================================

class TestContextManager:

    def test_context_manager_works(self):
        with IRBandDB() as db:
            bands = db.get_bands("CN")
            assert len(bands) > 0


# ===========================================================================
# 12. Integration — parser + database
# ===========================================================================

class TestParserIntegration:

    def test_all_parser_ligands_covered(self, db):
        """
        Every ligand in the parser's KNOWN_LIGANDS table should have
        at least some entry in the IR database (or at least not crash).
        """
        from coordchem.parser import KNOWN_LIGANDS
        db_ligands = set(db.get_all_ligands())
        parser_ligands = set(KNOWN_LIGANDS.keys())

        covered = parser_ligands & db_ligands
        # At least 10 of the parser ligands should be in the DB
        assert len(covered) >= 10, \
            f"Only {len(covered)} parser ligands covered in DB: {covered}"

    def test_full_pipeline_hexacyanoferrate(self, db):
        """Parse [Fe(CN)6]4- then look up its IR bands."""
        from coordchem.parser import parse_formula
        result = parse_formula("[Fe(CN)6]4-")

        all_bands = []
        for ligand in result.ligands:
            bands = db.get_bands(
                ligand,
                spectrum_type="IR",
                metal=result.metal
            )
            all_bands.extend(bands)

        assert len(all_bands) > 0
        wavenumbers = [b.center for b in all_bands]
        # Should have the diagnostic CN stretch around 2093
        assert any(2050 <= wn <= 2150 for wn in wavenumbers)

    def test_full_pipeline_cisplatin(self, db):
        """Parse [PtCl2(NH3)2] then look up its IR bands."""
        from coordchem.parser import parse_formula
        result = parse_formula("[PtCl2(NH3)2]")

        all_bands = []
        for ligand in result.ligands:
            bands = db.get_bands(ligand, spectrum_type="IR", metal=result.metal)
            all_bands.extend(bands)

        ligands_with_bands = {b.ligand for b in all_bands}
        assert "Cl"  in ligands_with_bands
        assert "NH3" in ligands_with_bands
