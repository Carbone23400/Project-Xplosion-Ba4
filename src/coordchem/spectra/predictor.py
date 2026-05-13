"""
coordchem/spectra/predictor.py
------------------------------
Assembles a list of IR or Raman bands for a coordination complex
by querying the band database for each ligand, then applying
chemistry-based corrections to improve wavenumber accuracy.

Pipeline
--------
    parse_formula("[Fe(CN)6]4-")   →   ParsedComplex
            ↓
    predict_spectrum(parsed)       →   PredictionResult
            ↓
    render_spectrum(result)        →   matplotlib figure

Three correction functions are applied in order after the database query:

    1. apply_backbonding_correction   — shifts π-accepting ligand bands
                                        based on metal and oxidation state
    2. apply_coordination_shift       — shifts free ligand bands to their
                                        coordinated positions
    3. apply_selection_rules          — removes symmetry-forbidden bands
                                        based on complex geometry

Usage
-----
    from coordchem.parser import parse_formula
    from coordchem.spectra.predictor import predict_spectrum

    parsed = parse_formula("[Fe(CN)6]4-")
    result = predict_spectrum(parsed, spectrum_type="IR")

    print(result)
    result.bands         # list of CorrectedBand objects
    result.intensities   # scaled numerical intensities
    result.warnings      # any missing ligands or correction notes
"""
import sys
import os

# Navigate up from src/coordchem/spectra/predictor.py to the repo root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

from dataclasses import dataclass, field
from typing import Optional

from data.ir_ra_bands import IRBandDB, BandRecord
from coordchem.parser import ParsedComplex
from coordchem.geometry import geometry_report

# ---------------------------------------------------------------------------
# Intensity mapping
# ---------------------------------------------------------------------------

INTENSITY_SCALE: dict[str, float] = {
    "very strong":  1.00,
    "strong":       0.75,
    "strong broad": 0.65,
    "medium":       0.50,
    "weak":         0.25,
}


# ---------------------------------------------------------------------------
# Correction data tables
# ---------------------------------------------------------------------------

# Backbonding shifts for π-accepting ligands (CN, CO, NO)
# Values in cm⁻¹ relative to the generic database center value.
# Source: Nakamoto Vol.2 pp. 155–172 (CN), pp. 120–154 (CO)

BACKBONDING_SHIFTS: dict[tuple, float] = {
    # (ligand, metal, oxidation_state) → shift in cm⁻¹
    # CN backbonding shifts
    ("CN", "Fe", 2): +10,
    ("CN", "Fe", 3): +35,
    ("CN", "Co", 2): +15,
    ("CN", "Co", 3): +30,
    ("CN", "Cr", 0): -20,
    ("CN", "Cr", 2): +10,
    ("CN", "Cr", 3): +25,
    ("CN", "Mn", 2): +15,
    ("CN", "Mn", 3): +30,
    ("CN", "Ni", 2): +5,
    ("CN", "Cu", 1): -10,
    ("CN", "Cu", 2): +10,
    ("CN", "Mo", 0): -25,
    ("CN", "W",  0): -30,
    ("CN", "Ru", 2): +5,
    ("CN", "Os", 2): +5,
    # CO backbonding shifts — more sensitive than CN
    ("CO", "Cr", 0): -30,
    ("CO", "Mo", 0): -25,
    ("CO", "W",  0): -20,
    ("CO", "Fe", 0): -35,
    ("CO", "Ni", 0): -40,
    ("CO", "Mn", 1): -15,
    ("CO", "Re", 1): -10,
    ("CO", "Ru", 0): -30,
    ("CO", "Os", 0): -25,
    # NO shifts
    ("NO", "Fe", 2): -20,
    ("NO", "Co", 3): -10,
    ("NO", "Cr", 0): -30,
}

# Coordination shifts — how much a band moves when a free ligand coordinates.
# Negative = redshift (moves to lower wavenumber).
# Source: Nakamoto Vol.2 Chapter 2

COORDINATION_SHIFTS: dict[tuple, float] = {
    # (ligand, assignment_keyword) → shift in cm⁻¹
    ("NH3", "N–H stretch"):          -150,
    ("NH3", "NH₃ sym. deformation"):  +50,
    ("NH3", "NH₃ rock"):              +80,
    ("H2O", "O–H stretch"):          -200,
    ("H2O", "H–O–H bend"):            +30,
    ("CO",  "C≡O stretch"):          -100,
    ("NO",  "N≡O stretch"):           -80,
    ("NO",  "N=O stretch"):           -60,
    ("dmso","S=O stretch"):          -100,
}

# Selection rules — which M-L band types are IR/Raman active per geometry.
# Based on the mutual exclusion rule:
# In centrosymmetric complexes (Oh, D4h), a band cannot be
# both IR and Raman active.
# Source: Nakamoto Vol.1 Chapter 1

SELECTION_RULES: dict[tuple, bool] = {
    # (geometry, assignment_keyword, spectrum_type) → is_active
    # Octahedral (Oh) — centrosymmetric
    ("octahedral", "M–C stretch",   "IR"):    False,
    ("octahedral", "M–C stretch",   "Raman"): True,
    ("octahedral", "M–N stretch",   "IR"):    False,
    ("octahedral", "M–N stretch",   "Raman"): True,
    ("octahedral", "M–O stretch",   "IR"):    False,
    ("octahedral", "M–O stretch",   "Raman"): True,
    ("octahedral", "M–Cl stretch",  "IR"):    False,
    ("octahedral", "M–Cl stretch",  "Raman"): True,
    ("octahedral", "M–Br stretch",  "IR"):    False,
    ("octahedral", "M–Br stretch",  "Raman"): True,
    # Square planar (D4h) — also centrosymmetric
    ("square planar", "M–N stretch",  "IR"):    False,
    ("square planar", "M–N stretch",  "Raman"): True,
    ("square planar", "M–Cl stretch", "IR"):    False,
    ("square planar", "M–Cl stretch", "Raman"): True,
    ("square planar", "M–C stretch",  "IR"):    False,
    ("square planar", "M–C stretch",  "Raman"): True,
    # Tetrahedral (Td) — no inversion center, all active
    ("tetrahedral", "M–N stretch",  "IR"):    True,
    ("tetrahedral", "M–N stretch",  "Raman"): True,
    ("tetrahedral", "M–Cl stretch", "IR"):    True,
    ("tetrahedral", "M–Cl stretch", "Raman"): True,
    ("tetrahedral", "M–O stretch",  "IR"):    True,
    ("tetrahedral", "M–O stretch",  "Raman"): True,
}


# ---------------------------------------------------------------------------
# CorrectedBand — a BandRecord with an adjusted wavenumber
# ---------------------------------------------------------------------------

@dataclass
class CorrectedBand:
    """
    A band after chemistry corrections have been applied.
    Wraps the original BandRecord and stores the corrected center.
    Delegates all attribute access to the original so it can be used
    anywhere BandRecord is expected.
    """
    original:           BandRecord
    corrected_center:   float
    correction_applied: float
    active:             bool

    @property
    def ligand(self):        return self.original.ligand
    @property
    def coordination(self):  return self.original.coordination
    @property
    def metal(self):         return self.original.metal
    @property
    def spectrum_type(self): return self.original.spectrum_type
    @property
    def wn_min(self):        return self.original.wn_min
    @property
    def wn_max(self):        return self.original.wn_max
    @property
    def intensity(self):     return self.original.intensity
    @property
    def assignment(self):    return self.original.assignment
    @property
    def source(self):        return self.original.source

    @property
    def center(self) -> float:
        """Always use the corrected center for plotting."""
        return self.corrected_center


# ---------------------------------------------------------------------------
# PredictionResult dataclass
# ---------------------------------------------------------------------------

@dataclass
class PredictionResult:
    """
    The output of predict_spectrum().

    Attributes
    ----------
    bands : list[CorrectedBand]
        All predicted bands after corrections, sorted by wavenumber.
    intensities : list[float]
        Scaled numerical intensity for each band (same order as bands).
    spectrum_type : str
        "IR" or "Raman".
    metal : str
        Metal symbol from the parsed complex.
    complex_formula : str
        Original formula string the user typed.
    ligand_coverage : dict[str, int]
        How many bands were found per ligand.
    warnings : list[str]
        Missing ligands, unusual corrections, etc.
    corrections_applied : int
        Total number of bands that received a wavenumber correction.
    bands_removed : int
        Number of bands removed by selection rules.
    """
    bands:               list
    intensities:         list[float]
    spectrum_type:       str
    metal:               str
    complex_formula:     str
    ligand_coverage:     dict[str, int] = field(default_factory=dict)
    warnings:            list[str]      = field(default_factory=list)
    corrections_applied: int            = 0
    bands_removed:       int            = 0

    @property
    def wavenumbers(self) -> list[float]:
        return [b.center for b in self.bands]

    @property
    def has_warnings(self) -> bool:
        return len(self.warnings) > 0

    @property
    def n_bands(self) -> int:
        return len(self.bands)

    def __str__(self) -> str:
        lines = [
            f"",
            f"  Predicted {self.spectrum_type} spectrum — {self.complex_formula}",
            f"  Metal: {self.metal}   |   {self.n_bands} bands predicted   |   "
            f"{self.corrections_applied} corrections applied   |   "
            f"{self.bands_removed} bands removed by selection rules",
            f"  {'─'*75}",
            f"  {'Ligand':<8} {'Coord.':<12} {'Center (cm⁻¹)':>14}  "
            f"{'Correction':>12}  {'Intensity':<14} Assignment",
            f"  {'─'*75}",
        ]
        for band in self.bands:
            corr_str = (f"{band.correction_applied:+.0f}"
                        if band.correction_applied != 0 else "—")
            lines.append(
                f"  {band.ligand:<8} {band.coordination:<12} "
                f"{band.center:>14.0f}  "
                f"{corr_str:>12}  "
                f"{band.intensity:<14} {band.assignment}"
            )
        if self.warnings:
            lines.append(f"")
            lines.append(f"  Warnings:")
            for w in self.warnings:
                lines.append(f"  ⚠  {w}")
        lines.append(f"")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main public function
# ---------------------------------------------------------------------------
 

def predict_spectrum(
    parsed:            ParsedComplex,
    spectrum_type:     str = "IR",
    db:                Optional[IRBandDB] = None,
    apply_corrections: bool = True,
) -> PredictionResult:
    """
    Predict the IR or Raman spectrum of a coordination complex.

    Parameters
    ----------
    parsed : ParsedComplex
        Output of parse_formula().
    spectrum_type : str
        "IR" or "Raman".
    db : IRBandDB, optional
        Existing database instance. If None a fresh one is created.
    apply_corrections : bool
        If True (default), apply backbonding, coordination shift, and
        selection rule corrections. Set to False to get raw database
        values for comparison or debugging.

    Returns
    -------
    PredictionResult
    """
    # Normalise spectrum type to match database values exactly
    if spectrum_type.upper() == "IR":
        spectrum_type = "IR"
    elif spectrum_type.upper() == "RAMAN":
        spectrum_type = "Raman"
    else:
        raise ValueError(
            f"spectrum_type must be 'IR' or 'Raman', got '{spectrum_type}'"
        )

    close_db = False
    if db is None:
        db       = IRBandDB()
        close_db = True

    all_bands:       list[CorrectedBand] = []
    all_intensities: list[float]         = []
    warnings:        list[str]           = []
    ligand_coverage: dict[str, int]      = {}
    bands_removed                        = 0

    # Geometry needed for selection rules
    try:
        report   = geometry_report(parsed)
        geometry = report["geometry"]
    except Exception:
        geometry = None   # if geometry fails, selection rules are skipped

    for ligand_formula, ligand_count in parsed.ligands.items():

        bands = db.get_bands(
            ligand        = ligand_formula,
            spectrum_type = spectrum_type,
            metal         = parsed.metal,
        )

        if not bands:
            warnings.append(
                f"No {spectrum_type} bands found for ligand "
                f"'{ligand_formula}' — it may not be in the database yet."
            )
            ligand_coverage[ligand_formula] = 0
            continue

        for band in bands:

            corrected = CorrectedBand(
                original           = band,
                corrected_center   = band.center,
                correction_applied = 0.0,
                active             = True,
            )

            if apply_corrections:
                corrected = apply_backbonding_correction(
                    corrected, parsed.metal, parsed.oxidation_state
                )
                corrected = apply_coordination_shift(
                    corrected, ligand_formula
                )
                corrected = apply_selection_rules(
                    corrected, geometry, spectrum_type
                )

            if not corrected.active:
                bands_removed += 1
                continue

            all_bands.append(corrected)
            all_intensities.append(_scale_intensity(band.intensity, ligand_count))

        ligand_coverage[ligand_formula] = len(bands)

    corrections_applied = sum(1 for b in all_bands if b.correction_applied != 0)

    # Sort by wavenumber
    if all_bands:
        paired = sorted(
            zip(all_bands, all_intensities), key=lambda x: x[0].center
        )
        all_bands, all_intensities = zip(*paired)
        all_bands       = list(all_bands)
        all_intensities = list(all_intensities)

    if close_db:
        db.close()

    return PredictionResult(
        bands               = all_bands,
        intensities         = all_intensities,
        spectrum_type       = spectrum_type,
        metal               = parsed.metal,
        complex_formula     = parsed.raw_formula,
        ligand_coverage     = ligand_coverage,
        warnings            = warnings,
        corrections_applied = corrections_applied,
        bands_removed       = bands_removed,
    )


# ---------------------------------------------------------------------------
# Correction function 1 — Backbonding
# ---------------------------------------------------------------------------

def apply_backbonding_correction(
    band:            CorrectedBand,
    metal:           str,
    oxidation_state: Optional[int],
) -> CorrectedBand:
    """
    Shift the wavenumber of π-accepting ligand stretch bands based on
    the degree of metal-to-ligand π-backbonding.

    π-accepting ligands (CN, CO, NO) accept electron density from the
    metal into their π* antibonding orbitals. This weakens the internal
    bond and lowers the stretching frequency. The effect is stronger for
    electron-rich metals (low oxidation state) and weaker for electron-
    poor metals (high oxidation state).

    Only applies to stretch bands of CN, CO, and NO ligands.
    """
    if oxidation_state is None:
        return band

    pi_accepting = {"CN", "CO", "NO"}
    if band.ligand not in pi_accepting:
        return band
    if "stretch" not in band.assignment.lower():
        return band

    shift = BACKBONDING_SHIFTS.get((band.ligand, metal, oxidation_state), 0.0)
    if shift == 0.0:
        return band

    return CorrectedBand(
        original           = band.original,
        corrected_center   = band.corrected_center + shift,
        correction_applied = band.correction_applied + shift,
        active             = band.active,
    )


# ---------------------------------------------------------------------------
# Correction function 2 — Coordination shift
# ---------------------------------------------------------------------------

def apply_coordination_shift(
    band:           CorrectedBand,
    ligand_formula: str,
) -> CorrectedBand:
    """
    Shift bands stored as free-ligand values to their coordinated positions.

    When a ligand binds to a metal, its vibrational frequencies shift
    from the free-ligand values. For example:
    - NH3 N-H stretch: ~3400 cm⁻¹ (free) → ~3250 cm⁻¹ (coordinated)
    - H2O O-H stretch: ~3600 cm⁻¹ (free) → ~3400 cm⁻¹ (coordinated)

    Only applies to bands with coordination mode "free" in the database,
    since coordinated values are already stored at the correct position.
    """
    if band.coordination != "free":
        return band

    shift = 0.0
    for (lig, keyword), value in COORDINATION_SHIFTS.items():
        if lig == ligand_formula and keyword in band.assignment:
            shift = value
            break

    if shift == 0.0:
        return band

    return CorrectedBand(
        original           = band.original,
        corrected_center   = band.corrected_center + shift,
        correction_applied = band.correction_applied + shift,
        active             = band.active,
    )


# ---------------------------------------------------------------------------
# Correction function 3 — Selection rules
# ---------------------------------------------------------------------------

def apply_selection_rules(
    band:          CorrectedBand,
    geometry:      Optional[str],
    spectrum_type: str,
) -> CorrectedBand:
    """
    Apply symmetry selection rules to determine whether a band is
    observable for a given geometry and spectrum type.

    The mutual exclusion principle: in complexes with a centre of
    inversion (octahedral Oh, square planar D4h), a vibrational mode
    cannot be both IR and Raman active simultaneously.

    Complexes without an inversion centre (tetrahedral Td, linear) have
    no such restriction.

    If geometry is unknown, the band is kept active (conservative approach —
    better to show too many bands than to silently remove real ones).
    """
    if geometry is None:
        return band

    geometry_lower = geometry.lower()

    for (geom, keyword, stype), is_active in SELECTION_RULES.items():
        if (geom == geometry_lower
                and keyword in band.assignment
                and stype == spectrum_type):
            return CorrectedBand(
                original           = band.original,
                corrected_center   = band.corrected_center,
                correction_applied = band.correction_applied,
                active             = is_active,
            )

    return band   # not in table — keep active


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _scale_intensity(intensity_label: str, ligand_count: int) -> float:
    """Convert intensity label to number and scale by ligand count."""
    base = INTENSITY_SCALE.get(intensity_label, 0.5)
    return base * ligand_count
