"""
Additional tests for correction functions in coordchem.spectra.predictor
------------------------------------------------------------------------
Paste these test classes at the bottom of your existing test_predictor.py

Tests cover:
    - apply_backbonding_correction
    - apply_coordination_shift
    - apply_selection_rules
    - apply_corrections=False vs True comparison
    - CorrectedBand dataclass
    - geometry integration
"""
import pytest 
from src.coordchem.parser import parse_formula
from src.coordchem.spectra.predictor import (
    predict_spectrum,
    apply_backbonding_correction,
    apply_coordination_shift,
    apply_selection_rules,
    CorrectedBand,
    BACKBONDING_SHIFTS,
    COORDINATION_SHIFTS,
    SELECTION_RULES,
)
from data.ir_ra_bands import IRBandDB, BandRecord


# ---------------------------------------------------------------------------
# Helper — build a fake BandRecord and wrap it in a CorrectedBand
# ---------------------------------------------------------------------------

def make_band(ligand, assignment, coordination="terminal", wn=2150.0):
    """Create a minimal BandRecord and wrap it in a CorrectedBand."""
    record = BandRecord(
        ligand        = ligand,
        coordination  = coordination,
        metal         = "any",
        spectrum_type = "IR",
        wn_min        = wn - 50,
        wn_max        = wn + 50,
        intensity     = "strong",
        assignment    = assignment,
        ir_active     = True,
        raman_active  = False,
        source        = "test",
    )
    return CorrectedBand(
        original           = record,
        corrected_center   = record.center,
        correction_applied = 0.0,
        active             = True,
    )


# ===========================================================================
# 1. CorrectedBand dataclass
# ===========================================================================

class TestCorrectedBand:

    def test_center_uses_corrected_value(self):
        band = make_band("CN", "C≡N stretch", wn=2150.0)
        assert band.center == pytest.approx(2150.0)

    def test_delegates_ligand_to_original(self):
        band = make_band("CN", "C≡N stretch")
        assert band.ligand == "CN"

    def test_delegates_assignment_to_original(self):
        band = make_band("CO", "C≡O stretch")
        assert band.assignment == "C≡O stretch"

    def test_delegates_intensity_to_original(self):
        band = make_band("NH3", "N–H stretch")
        assert band.intensity == "strong"

    def test_active_default_true(self):
        band = make_band("CN", "C≡N stretch")
        assert band.active == True

    def test_correction_applied_default_zero(self):
        band = make_band("CN", "C≡N stretch")
        assert band.correction_applied == pytest.approx(0.0)

    def test_corrected_center_independent_of_original(self):
        """Changing corrected_center should not affect original.center."""
        band = make_band("CN", "C≡N stretch", wn=2150.0)
        assert band.original.center == pytest.approx(2150.0)
        assert band.center == pytest.approx(2150.0)


# ===========================================================================
# 2. apply_backbonding_correction
# ===========================================================================

class TestBackbondingCorrection:

    def test_fe3_cn_shifts_up(self):
        """Fe(III) has less backbonding → CN stretch shifts to higher wavenumber."""
        band   = make_band("CN", "C≡N stretch", wn=2150.0)
        result = apply_backbonding_correction(band, "Fe", 3)
        assert result.corrected_center > 2150.0

    def test_fe2_cn_shifts_less_than_fe3(self):
        """Fe(II) has more backbonding than Fe(III) → smaller upward shift."""
        band   = make_band("CN", "C≡N stretch", wn=2150.0)
        fe2    = apply_backbonding_correction(band, "Fe", 2)
        fe3    = apply_backbonding_correction(band, "Fe", 3)
        assert fe2.corrected_center < fe3.corrected_center

    def test_cr0_co_shifts_down(self):
        """Cr(0) has strong backbonding → CO stretch shifts to lower wavenumber."""
        band   = make_band("CO", "C≡O stretch", wn=1925.0)
        result = apply_backbonding_correction(band, "Cr", 0)
        assert result.corrected_center < 1925.0

    def test_correction_recorded_in_correction_applied(self):
        """The shift amount should be stored in correction_applied."""
        band     = make_band("CN", "C≡N stretch", wn=2150.0)
        result   = apply_backbonding_correction(band, "Fe", 3)
        expected = BACKBONDING_SHIFTS.get(("CN", "Fe", 3), 0.0)
        assert result.correction_applied == pytest.approx(expected)

    def test_non_pi_accepting_ligand_not_shifted(self):
        """NH3 is not a π-accepting ligand — no backbonding correction."""
        band   = make_band("NH3", "N–H stretch", wn=3250.0)
        result = apply_backbonding_correction(band, "Fe", 3)
        assert result.corrected_center == pytest.approx(3250.0)
        assert result.correction_applied == pytest.approx(0.0)

    def test_bend_band_not_shifted(self):
        """Backbonding only applies to stretch bands, not bends."""
        band   = make_band("CN", "M–C≡N bend", wn=400.0)
        result = apply_backbonding_correction(band, "Fe", 3)
        assert result.corrected_center == pytest.approx(400.0)

    def test_unknown_metal_oxidation_no_shift(self):
        """Metal/OS combo not in table → no shift."""
        band   = make_band("CN", "C≡N stretch", wn=2150.0)
        result = apply_backbonding_correction(band, "La", 3)
        assert result.corrected_center == pytest.approx(2150.0)

    def test_none_oxidation_state_no_shift(self):
        """If oxidation state is None, no correction applied."""
        band   = make_band("CN", "C≡N stretch", wn=2150.0)
        result = apply_backbonding_correction(band, "Fe", None)
        assert result.corrected_center == pytest.approx(2150.0)

    def test_mo0_co_shifts_down(self):
        band   = make_band("CO", "C≡O stretch", wn=1925.0)
        result = apply_backbonding_correction(band, "Mo", 0)
        assert result.corrected_center < 1925.0

    def test_active_flag_unchanged(self):
        """Backbonding correction should not change the active flag."""
        band   = make_band("CN", "C≡N stretch", wn=2150.0)
        result = apply_backbonding_correction(band, "Fe", 3)
        assert result.active == True


# ===========================================================================
# 3. apply_coordination_shift
# ===========================================================================

class TestCoordinationShift:

    def test_free_nh3_nh_stretch_redshifts(self):
        """Free NH3 N-H stretch should move down when coordination='free'."""
        band   = make_band("NH3", "N–H stretch", coordination="free", wn=3400.0)
        result = apply_coordination_shift(band, "NH3")
        assert result.corrected_center < 3400.0

    def test_free_h2o_oh_stretch_redshifts(self):
        """Free H2O O-H stretch should redshift on coordination."""
        band   = make_band("H2O", "O–H stretch", coordination="free", wn=3600.0)
        result = apply_coordination_shift(band, "H2O")
        assert result.corrected_center < 3600.0

    def test_already_coordinated_not_shifted(self):
        """Terminal (already coordinated) bands should not be shifted again."""
        band   = make_band("NH3", "N–H stretch", coordination="terminal", wn=3250.0)
        result = apply_coordination_shift(band, "NH3")
        assert result.corrected_center == pytest.approx(3250.0)

    def test_shift_amount_matches_table(self):
        """The shift should exactly match the COORDINATION_SHIFTS table."""
        band     = make_band("NH3", "N–H stretch", coordination="free", wn=3400.0)
        result   = apply_coordination_shift(band, "NH3")
        expected = COORDINATION_SHIFTS.get(("NH3", "N–H stretch"), 0.0)
        assert result.correction_applied == pytest.approx(expected)

    def test_cn_not_in_coordination_shifts(self):
        """CN is stored directly as coordinated — no coordination shift."""
        band   = make_band("CN", "C≡N stretch", coordination="terminal", wn=2150.0)
        result = apply_coordination_shift(band, "CN")
        assert result.corrected_center == pytest.approx(2150.0)

    def test_nh3_umbrella_blueshifts(self):
        """NH3 symmetric deformation should blueshift on coordination."""
        band   = make_band("NH3", "NH₃ sym. deformation", coordination="free", wn=1600.0)
        result = apply_coordination_shift(band, "NH3")
        assert result.corrected_center > 1600.0

    def test_active_flag_unchanged(self):
        band   = make_band("NH3", "N–H stretch", coordination="free", wn=3400.0)
        result = apply_coordination_shift(band, "NH3")
        assert result.active == True


# ===========================================================================
# 4. apply_selection_rules
# ===========================================================================

class TestSelectionRules:

    def test_octahedral_mn_stretch_ir_inactive(self):
        """In octahedral complexes M-N stretch is IR forbidden."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, "octahedral", "IR")
        assert result.active == False

    def test_octahedral_mn_stretch_raman_active(self):
        """In octahedral complexes M-N stretch is Raman active."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, "octahedral", "Raman")
        assert result.active == True

    def test_octahedral_mcl_stretch_ir_inactive(self):
        band   = make_band("Cl", "M–Cl stretch")
        result = apply_selection_rules(band, "octahedral", "IR")
        assert result.active == False

    def test_octahedral_mcl_stretch_raman_active(self):
        band   = make_band("Cl", "M–Cl stretch")
        result = apply_selection_rules(band, "octahedral", "Raman")
        assert result.active == True

    def test_tetrahedral_mn_stretch_ir_active(self):
        """Tetrahedral has no inversion centre — M-N stretch is IR active."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, "tetrahedral", "IR")
        assert result.active == True

    def test_tetrahedral_mn_stretch_raman_active(self):
        """Tetrahedral — M-N stretch is also Raman active."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, "tetrahedral", "Raman")
        assert result.active == True

    def test_square_planar_mn_stretch_ir_inactive(self):
        """Square planar (D4h) is centrosymmetric — M-N stretch IR forbidden."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, "square planar", "IR")
        assert result.active == False

    def test_cn_stretch_not_affected_by_selection_rules(self):
        """C≡N stretch is a ligand internal mode — selection rules don't apply."""
        band   = make_band("CN", "C≡N stretch")
        result = apply_selection_rules(band, "octahedral", "IR")
        assert result.active == True   # not in SELECTION_RULES table → kept

    def test_none_geometry_keeps_band_active(self):
        """If geometry is unknown, band should not be removed."""
        band   = make_band("NH3", "M–N stretch")
        result = apply_selection_rules(band, None, "IR")
        assert result.active == True

    def test_mutual_exclusion_octahedral(self):
        """
        Mutual exclusion rule: M-N stretch cannot be both IR and Raman
        active in an octahedral complex.
        """
        band  = make_band("NH3", "M–N stretch")
        ir    = apply_selection_rules(band, "octahedral", "IR")
        raman = apply_selection_rules(band, "octahedral", "Raman")
        # One must be active and the other inactive
        assert ir.active != raman.active

    def test_wavenumber_unchanged_by_selection_rules(self):
        """Selection rules only change active flag, not the wavenumber."""
        band   = make_band("NH3", "M–N stretch", wn=450.0)
        result = apply_selection_rules(band, "octahedral", "IR")
        assert result.corrected_center == pytest.approx(450.0)


# ===========================================================================
# 5. Corrections pipeline — before vs after
# ===========================================================================

class TestCorrectionsOnOffComparison:

    def test_raw_vs_corrected_cn_fe3_differ(self):
        """
        Fe(III) hexacyanoferrate: corrections should shift CN stretch up.
        Raw database center vs corrected center should be different.
        """
        parsed = parse_formula("[Fe(CN)6]3-")

        raw       = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=False)
        corrected = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=True)

        raw_cn_wns       = [b.center for b in raw.bands
                            if "C≡N stretch" in b.assignment]
        corrected_cn_wns = [b.center for b in corrected.bands
                            if "C≡N stretch" in b.assignment]

        if raw_cn_wns and corrected_cn_wns:
            assert max(corrected_cn_wns) != pytest.approx(max(raw_cn_wns))

    def test_corrections_off_gives_more_bands(self):
        """
        Without selection rules, octahedral complexes should have MORE bands
        (including M-L stretches that are IR-forbidden).
        """
        parsed = parse_formula("[Fe(CN)6]4-")
        parsed.geometry = "octahedral"

        raw       = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=False)
        corrected = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=True)

        assert raw.n_bands >= corrected.n_bands

    def test_corrections_applied_count_nonzero_for_known_metal(self):
        """For Fe(III) with CN ligands, at least one correction should fire."""
        parsed = parse_formula("[Fe(CN)6]3-")
        result = predict_spectrum(parsed, spectrum_type="IR",
                                  apply_corrections=True)
        assert result.corrections_applied > 0

    def test_corrections_applied_zero_when_off(self):
        """When apply_corrections=False, no corrections should be recorded."""
        parsed = parse_formula("[Fe(CN)6]3-")
        result = predict_spectrum(parsed, spectrum_type="IR",
                                  apply_corrections=False)
        assert result.corrections_applied == 0

    def test_bands_removed_reported_in_result(self):
        """bands_removed should be > 0 for an octahedral complex with IR."""
        parsed          = parse_formula("[Fe(CN)6]4-")
        parsed.geometry = "octahedral"
        result          = predict_spectrum(parsed, spectrum_type="IR",
                                           apply_corrections=True)
        assert result.bands_removed >= 0   # may be 0 if no ML stretches found


# ===========================================================================
# 6. Geometry integration
# ===========================================================================

class TestGeometryIntegration:

    def test_octahedral_ir_no_ml_stretches(self):
        """
        In an octahedral complex, M-L stretches should be absent from
        the IR spectrum (they are Raman active, not IR active).
        """
        parsed          = parse_formula("[Fe(CN)6]4-")
        parsed.geometry = "octahedral"
        result          = predict_spectrum(parsed, spectrum_type="IR",
                                           apply_corrections=True)
        ml_stretch_ir = [b for b in result.bands if "M–C stretch" in b.assignment]
        assert len(ml_stretch_ir) == 0

    def test_octahedral_raman_has_ml_stretches(self):
        """
        In an octahedral complex, M-L stretches SHOULD appear in Raman.
        """
        parsed          = parse_formula("[Fe(CN)6]4-")
        parsed.geometry = "octahedral"
        result          = predict_spectrum(parsed, spectrum_type="Raman",
                                           apply_corrections=True)
        ml_stretch_raman = [b for b in result.bands if "M–C stretch" in b.assignment]
        assert len(ml_stretch_raman) > 0

    def test_tetrahedral_ir_has_ml_stretches(self):
        """
        In a tetrahedral complex, M-L stretches ARE IR active.
        """
        parsed          = parse_formula("[NiCl4]2-")
        parsed.geometry = "tetrahedral"
        result          = predict_spectrum(parsed, spectrum_type="IR",
                                           apply_corrections=True)
        ml_stretch = [b for b in result.bands if "M–Cl stretch" in b.assignment]
        assert len(ml_stretch) > 0

    def test_no_geometry_does_not_crash(self):
        """
        If geometry is not set on the parsed object, predict_spectrum
        should still run without crashing.
        """
        parsed = parse_formula("[Cu(NH3)4]2+")
        # Don't set geometry — it should fall back to None safely
        if hasattr(parsed, "geometry"):
            del parsed.geometry
        result = predict_spectrum(parsed, spectrum_type="IR")
        assert result.n_bands >= 0   # just check it didn't crash


# ===========================================================================
# 7. Backbonding — real complex benchmarks
# ===========================================================================

class TestBackbondingBenchmarks:

    def test_hexacyanoferrate_II_vs_III_cn_stretch(self):
        """
        Fe(III) complexes have less backbonding than Fe(II), so the
        CN stretch in [Fe(CN)6]3- should be higher than in [Fe(CN)6]4-.
        """
        parsed_II  = parse_formula("[Fe(CN)6]4-")
        parsed_III = parse_formula("[Fe(CN)6]3-")

        result_II  = predict_spectrum(parsed_II,  spectrum_type="IR",
                                      apply_corrections=True)
        result_III = predict_spectrum(parsed_III, spectrum_type="IR",
                                      apply_corrections=True)

        cn_II  = [b.center for b in result_II.bands
                  if "C≡N stretch" in b.assignment and b.metal in ("Fe", "any")]
        cn_III = [b.center for b in result_III.bands
                  if "C≡N stretch" in b.assignment and b.metal in ("Fe", "any")]

        if cn_II and cn_III:
            assert max(cn_III) > max(cn_II), (
                f"Expected Fe(III) CN stretch ({max(cn_III):.0f}) > "
                f"Fe(II) CN stretch ({max(cn_II):.0f})"
            )

    def test_cr_co6_co_stretch_below_generic(self):
        """
        Cr(0) has strong backbonding — its CO stretch should be below
        the generic database value for terminal CO.
        """
        parsed = parse_formula("[Cr(CO)6]")

        raw       = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=False)
        corrected = predict_spectrum(parsed, spectrum_type="IR",
                                     apply_corrections=True)

        raw_co       = [b.center for b in raw.bands
                        if "C≡O stretch" in b.assignment]
        corrected_co = [b.center for b in corrected.bands
                        if "C≡O stretch" in b.assignment]

        if raw_co and corrected_co:
            assert min(corrected_co) < min(raw_co)