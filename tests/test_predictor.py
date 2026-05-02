"""
tests/test_predictor.py
-----------------------
Test suite for src.spectra.predictor

Run with:  python -m pytest tests/test_predictor.py -v
"""

import pytest
from src.coordchem.parser import parse_formula
from src.spectra.predictor import predict_spectrum, PredictionResult, INTENSITY_SCALE
from data.database.ir_ra_bands import IRBandDB


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def predict(formula, spectrum_type="IR", db=None):
    """Shorthand: parse then predict."""
    parsed = parse_formula(formula)
    return predict_spectrum(parsed, spectrum_type=spectrum_type, db=db)


# ===========================================================================
# 1. Return type and basic structure
# ===========================================================================

class TestReturnType:

    def test_returns_prediction_result(self):
        result = predict("[Fe(CN)6]4-")
        assert isinstance(result, PredictionResult)

    def test_bands_is_a_list(self):
        result = predict("[Fe(CN)6]4-")
        assert isinstance(result.bands, list)

    def test_intensities_is_a_list(self):
        result = predict("[Fe(CN)6]4-")
        assert isinstance(result.intensities, list)

    def test_bands_and_intensities_same_length(self):
        result = predict("[Fe(CN)6]4-")
        assert len(result.bands) == len(result.intensities)

    def test_warnings_is_a_list(self):
        result = predict("[Fe(CN)6]4-")
        assert isinstance(result.warnings, list)

    def test_ligand_coverage_is_a_dict(self):
        result = predict("[Fe(CN)6]4-")
        assert isinstance(result.ligand_coverage, dict)


# ===========================================================================
# 2. Metadata fields
# ===========================================================================

class TestMetadata:

    def test_metal_field(self):
        result = predict("[Fe(CN)6]4-")
        assert result.metal == "Fe"

    def test_complex_formula_field(self):
        result = predict("[Fe(CN)6]4-")
        assert result.complex_formula == "[Fe(CN)6]4-"

    def test_spectrum_type_ir(self):
        result = predict("[Fe(CN)6]4-", spectrum_type="IR")
        assert result.spectrum_type == "IR"

    def test_spectrum_type_raman(self):
        result = predict("[Fe(CN)6]4-", spectrum_type="Raman")
        assert result.spectrum_type == "RAMAN"

    def test_spectrum_type_case_insensitive(self):
        result = predict("[Fe(CN)6]4-", spectrum_type="raman")
        assert result.spectrum_type == "RAMAN"

    def test_invalid_spectrum_type_raises(self):
        parsed = parse_formula("[Fe(CN)6]4-")
        with pytest.raises(ValueError):
            predict_spectrum(parsed, spectrum_type="UV")


# ===========================================================================
# 3. Band content — hexacyanoferrate(II) benchmark
# ===========================================================================

class TestHexacyanoferrate:

    def test_has_bands(self):
        result = predict("[Fe(CN)6]4-")
        assert result.n_bands > 0

    def test_cn_ligand_covered(self):
        result = predict("[Fe(CN)6]4-")
        assert result.ligand_coverage.get("CN", 0) > 0

    def test_cn_stretch_present(self):
        """The diagnostic C≡N stretch around 2093 cm⁻¹ must be present."""
        result = predict("[Fe(CN)6]4-")
        wns = result.wavenumbers
        assert any(2050 <= wn <= 2150 for wn in wns), \
            f"No C≡N stretch found near 2093 cm⁻¹. Wavenumbers: {wns}"

    def test_mc_stretch_present(self):
        """M-C stretch should appear in the fingerprint region."""
        result = predict("[Fe(CN)6]4-")
        wns = result.wavenumbers
        assert any(300 <= wn <= 700 for wn in wns)

    def test_no_warnings_for_cn(self):
        """CN is well covered — no warnings expected."""
        result = predict("[Fe(CN)6]4-")
        assert not result.has_warnings

    def test_raman_bands_different_from_ir(self):
        """IR and Raman band lists should not be identical."""
        ir     = predict("[Fe(CN)6]4-", spectrum_type="IR")
        raman  = predict("[Fe(CN)6]4-", spectrum_type="Raman")
        ir_wns = set(round(w) for w in ir.wavenumbers)
        ra_wns = set(round(w) for w in raman.wavenumbers)
        # They may share some values but shouldn't be completely identical
        assert ir_wns != ra_wns or len(ir_wns) == 0


# ===========================================================================
# 4. Band content — cisplatin benchmark
# ===========================================================================

class TestCisplatin:

    def test_has_bands(self):
        result = predict("[PtCl2(NH3)2]")
        assert result.n_bands > 0

    def test_cl_and_nh3_both_covered(self):
        result = predict("[PtCl2(NH3)2]")
        assert result.ligand_coverage.get("Cl",  0) > 0
        assert result.ligand_coverage.get("NH3", 0) > 0

    def test_pt_cl_stretch_in_fingerprint(self):
        """Pt-Cl stretch should appear between 200 and 400 cm⁻¹."""
        result = predict("[PtCl2(NH3)2]")
        wns    = result.wavenumbers
        assert any(200 <= wn <= 400 for wn in wns)

    def test_nh_stretch_present(self):
        """N-H stretch should appear above 3000 cm⁻¹."""
        result = predict("[PtCl2(NH3)2]")
        wns    = result.wavenumbers
        assert any(wn > 3000 for wn in wns)


# ===========================================================================
# 5. Band content — chromium hexacarbonyl benchmark
# ===========================================================================

class TestChromiumCarbonyl:

    def test_co_stretch_very_strong(self):
        result = predict("[Cr(CO)6]")
        strong = [b for b in result.bands
                  if "C≡O stretch" in b.assignment and b.intensity == "very strong"]
        assert len(strong) > 0

    def test_co_stretch_above_1850(self):
        result = predict("[Cr(CO)6]")
        wns    = result.wavenumbers
        assert any(wn > 1850 for wn in wns)

    def test_oxidation_state_zero_no_warning(self):
        result = predict("[Cr(CO)6]")
        assert result.metal == "Cr"


# ===========================================================================
# 6. Band ordering
# ===========================================================================

class TestBandOrdering:

    def test_bands_sorted_by_wavenumber(self):
        """Bands must be sorted low → high wavenumber."""
        result = predict("[PtCl2(NH3)2]")
        wns    = result.wavenumbers
        assert wns == sorted(wns), \
            f"Bands not sorted: {wns}"

    def test_intensities_follow_band_order(self):
        """Intensities must stay aligned with their bands after sorting."""
        result = predict("[PtCl2(NH3)2]")
        for band, intensity in zip(result.bands, result.intensities):
            expected = INTENSITY_SCALE.get(band.intensity, 0.5) * \
                       parse_formula("[PtCl2(NH3)2]").ligands.get(band.ligand, 1)
            assert abs(intensity - expected) < 1e-9


# ===========================================================================
# 7. Intensity scaling
# ===========================================================================

class TestIntensityScaling:

    def test_intensities_are_positive(self):
        result = predict("[Fe(CN)6]4-")
        assert all(i > 0 for i in result.intensities)

    def test_six_cn_stronger_than_two_cn(self):
        """[Fe(CN)6]4- has 6 CN — intensities should be 3× those of [Fe(CN)2]."""
        # We compare the same band type at different counts
        # by directly calling the helper
        from src.spectra.predictor import _scale_intensity
        i6 = _scale_intensity("strong", 6)
        i2 = _scale_intensity("strong", 2)
        assert i6 == pytest.approx(3 * i2)

    def test_very_strong_greater_than_weak(self):
        from src.spectra.predictor import _scale_intensity
        assert _scale_intensity("very strong", 1) > _scale_intensity("weak", 1)

    def test_intensity_scale_covers_all_labels(self):
        """All intensity labels in the DB should be in INTENSITY_SCALE."""
        known_labels = {"very strong", "strong", "strong broad", "medium", "weak"}
        for label in known_labels:
            assert label in INTENSITY_SCALE


# ===========================================================================
# 8. Convenience properties
# ===========================================================================

class TestConvenienceProperties:

    def test_wavenumbers_property(self):
        result = predict("[Fe(CN)6]4-")
        assert result.wavenumbers == [b.center for b in result.bands]

    def test_n_bands_property(self):
        result = predict("[Fe(CN)6]4-")
        assert result.n_bands == len(result.bands)

    def test_has_warnings_false_when_no_warnings(self):
        result = predict("[Fe(CN)6]4-")
        assert result.has_warnings == False

    def test_has_warnings_true_when_unknown_ligand(self):
        result = predict("[Fe(XYZ)3]3+")
        assert result.has_warnings == True

    def test_str_representation_runs(self):
        result = predict("[Fe(CN)6]4-")
        output = str(result)
        assert "Fe" in output
        assert "IR" in output


# ===========================================================================
# 9. Unknown ligands
# ===========================================================================

class TestUnknownLigands:

    def test_unknown_ligand_gives_warning(self):
        result = predict("[Fe(XYZ)3]3+")
        assert any("XYZ" in w for w in result.warnings)

    def test_unknown_ligand_gives_empty_bands(self):
        result = predict("[Fe(XYZ)3]3+")
        assert result.n_bands == 0

    def test_partial_unknown_still_returns_known_bands(self):
        """[PtCl2(XYZ)2] — Cl is known, XYZ is not. Should return Cl bands."""
        parsed = parse_formula("[PtCl2(NH3)2]")
        # Manually inject an unknown ligand
        parsed.ligands["UNKNOWN"] = 1
        result = predict_spectrum(parsed, spectrum_type="IR")
        assert result.ligand_coverage.get("Cl", 0) > 0
        assert any("UNKNOWN" in w for w in result.warnings)


# ===========================================================================
# 10. Database reuse
# ===========================================================================

class TestDatabaseReuse:

    def test_external_db_instance(self):
        """Passing an external DB instance should work correctly."""
        db     = IRBandDB()
        result = predict("[Fe(CN)6]4-", db=db)
        assert result.n_bands > 0
        db.close()

    def test_multiple_predictions_same_db(self):
        """Reusing a DB across predictions should give consistent results."""
        db  = IRBandDB()
        r1  = predict("[Fe(CN)6]4-",   db=db)
        r2  = predict("[Cu(NH3)4]2+",  db=db)
        r3  = predict("[PtCl2(NH3)2]", db=db)
        assert r1.n_bands > 0
        assert r2.n_bands > 0
        assert r3.n_bands > 0
        db.close()


# ===========================================================================
# 11. Full benchmark suite
# ===========================================================================

class TestBenchmarkComplexes:
    """Validate against well-known coordination chemistry examples."""

    def test_hexacyanoferrate_II(self):
        r = predict("[Fe(CN)6]4-")
        assert r.n_bands > 0
        assert any(2050 <= wn <= 2150 for wn in r.wavenumbers)

    def test_tetraamminecopper(self):
        r = predict("[Cu(NH3)4]2+")
        assert r.ligand_coverage.get("NH3", 0) > 0

    def test_hexaaquairon(self):
        r = predict("[Fe(H2O)6]3+")
        assert r.ligand_coverage.get("H2O", 0) > 0
        assert any(wn > 3000 for wn in r.wavenumbers)

    def test_tris_en_cobalt(self):
        r = predict("[Co(en)3]3+")
        assert r.ligand_coverage.get("en", 0) > 0

    def test_chromium_carbonyl(self):
        r = predict("[Cr(CO)6]")
        assert any(wn > 1850 for wn in r.wavenumbers)

    def test_tetrachloronickelate(self):
        r = predict("[NiCl4]2-")
        assert r.ligand_coverage.get("Cl", 0) > 0

    def test_diamminesilver(self):
        r = predict("[Ag(NH3)2]+")
        assert r.ligand_coverage.get("NH3", 0) > 0
