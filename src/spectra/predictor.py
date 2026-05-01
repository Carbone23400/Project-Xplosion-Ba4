"""
src/spectra/predictor.py
------------------------------
Assembles a list of IR or Raman bands for a coordination complex
by querying the band database for each ligand.

This module sits between the parser and the renderer:

    parse_formula("[Fe(CN)6]4-")   →   ParsedComplex
            ↓
    predict_spectrum(parsed)       →   PredictionResult
            ↓
    render_spectrum(result)        →   matplotlib figure

Usage
-----
    from src.coordchem.parser import parse_formula
    from src.spectra.predictor import predict_spectrum

    parsed = parse_formula("[Fe(CN)6]4-")
    result = predict_spectrum(parsed, spectrum_type="IR")

    print(result)              # summary table
    result.bands               # list of BandRecord objects
    result.intensities         # scaled numerical intensities
    result.warnings            # any missing ligands
"""

from dataclasses import dataclass, field
from typing import Optional

from data.database.ir_ra_bands import IRBandDB, BandRecord
from src.coordchem.parser import ParsedComplex


# ---------------------------------------------------------------------------
# Intensity mapping — converts text label to a 0–1 numerical scale
# ---------------------------------------------------------------------------

INTENSITY_SCALE: dict[str, float] = {
    "very strong": 1.00,
    "strong":      0.75,
    "strong broad":0.65,
    "medium":      0.50,
    "weak":        0.25,
}


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class PredictionResult:
    """
    The output of predict_spectrum().

    Attributes
    ----------
    bands : list[BandRecord]
        All predicted bands, one BandRecord per peak.
    intensities : list[float]
        Scaled numerical intensity for each band (same order as bands).
        Already accounts for how many of that ligand are in the complex
        (e.g. 6× CN gives stronger peaks than 2× CN).
    spectrum_type : str
        "IR" or "Raman".
    metal : str
        Metal symbol from the parsed complex.
    complex_formula : str
        Original formula string the user typed.
    ligand_coverage : dict[str, int]
        How many bands were found per ligand.
        e.g. {"CN": 5, "NH3": 4} — useful for debugging.
    warnings : list[str]
        One warning per ligand that had no bands in the database.
    """
    bands:            list[BandRecord]
    intensities:      list[float]
    spectrum_type:    str
    metal:            str
    complex_formula:  str
    ligand_coverage:  dict[str, int]  = field(default_factory=dict)
    warnings:         list[str]       = field(default_factory=list)

    # ------------------------------------------------------------------
    # Convenience properties
    # ------------------------------------------------------------------

    @property
    def wavenumbers(self) -> list[float]:
        """Center wavenumber of each band — shortcut for the renderer."""
        return [b.center for b in self.bands]

    @property
    def has_warnings(self) -> bool:
        return len(self.warnings) > 0

    @property
    def n_bands(self) -> int:
        return len(self.bands)

    # ------------------------------------------------------------------
    # Pretty printing
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        lines = [
            f"",
            f"  Predicted {self.spectrum_type} spectrum — {self.complex_formula}",
            f"  Metal: {self.metal}   |   {self.n_bands} bands predicted",
            f"  {'─'*65}",
            f"  {'Ligand':<8} {'Coord.':<12} {'Center (cm⁻¹)':>14}  "
            f"{'Intensity':<14} Assignment",
            f"  {'─'*65}",
        ]
        for band, intensity in zip(self.bands, self.intensities):
            lines.append(
                f"  {band.ligand:<8} {band.coordination:<12} "
                f"{band.center:>14.0f}  "
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
    parsed:        ParsedComplex,
    spectrum_type: str = "IR",
    db:            Optional[IRBandDB] = None,
) -> PredictionResult:
    """
    Predict the IR or Raman spectrum of a coordination complex.

    Parameters
    ----------
    parsed : ParsedComplex
        Output of parse_formula(). Contains the metal, ligands, and
        oxidation state needed to query the database.
    spectrum_type : str
        "IR" or "Raman". Default is "IR".
    db : IRBandDB, optional
        An existing database instance to use. If None, a fresh
        in-memory database is created. Pass your own instance if you
        want to reuse a connection across multiple calls (more efficient).

    Returns
    -------
    PredictionResult
        Contains the band list, scaled intensities, and any warnings.

    Examples
    --------
    >>> from coordchem.parser import parse_formula
    >>> from coordchem.spectra.predictor import predict_spectrum
    >>> parsed = parse_formula("[Fe(CN)6]4-")
    >>> result = predict_spectrum(parsed, spectrum_type="IR")
    >>> print(result)
    """
    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    spectrum_type = spectrum_type.upper()
    if spectrum_type not in ("IR", "RAMAN"):
        raise ValueError(
            f"spectrum_type must be 'IR' or 'Raman', got '{spectrum_type}'"
        )

    # ------------------------------------------------------------------
    # Open database if not provided
    # ------------------------------------------------------------------
    close_db = False
    if db is None:
        db       = IRBandDB()
        close_db = True   # we opened it, so we close it when done

    # ------------------------------------------------------------------
    # Loop over every ligand in the complex
    # ------------------------------------------------------------------
    all_bands:       list[BandRecord] = []
    all_intensities: list[float]      = []
    warnings:        list[str]        = []
    ligand_coverage: dict[str, int]   = {}

    for ligand_formula, ligand_count in parsed.ligands.items():

        # Query the database for this ligand
        bands = db.get_bands(
            ligand        = ligand_formula,
            spectrum_type = spectrum_type,
            metal         = parsed.metal,
        )

        # Warn if nothing was found
        if not bands:
            warnings.append(
                f"No {spectrum_type} bands found for ligand "
                f"'{ligand_formula}' — it may not be in the database yet."
            )
            ligand_coverage[ligand_formula] = 0
            continue

        # Scale each band's intensity by how many of this ligand there are
        # (6× CN should give stronger peaks than 2× CN)
        for band in bands:
            scaled = _scale_intensity(band.intensity, ligand_count)
            all_bands.append(band)
            all_intensities.append(scaled)

        ligand_coverage[ligand_formula] = len(bands)

    # ------------------------------------------------------------------
    # Sort everything by wavenumber (low → high) for clean plotting
    # ------------------------------------------------------------------
    if all_bands:
        paired     = sorted(zip(all_bands, all_intensities), key=lambda x: x[0].center)
        all_bands, all_intensities = zip(*paired)
        all_bands       = list(all_bands)
        all_intensities = list(all_intensities)

    # ------------------------------------------------------------------
    # Close db if we opened it
    # ------------------------------------------------------------------
    if close_db:
        db.close()

    # ------------------------------------------------------------------
    # Build and return the result
    # ------------------------------------------------------------------
    return PredictionResult(
        bands           = all_bands,
        intensities     = all_intensities,
        spectrum_type   = spectrum_type,
        metal           = parsed.metal,
        complex_formula = parsed.raw_formula,
        ligand_coverage = ligand_coverage,
        warnings        = warnings,
    )


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _scale_intensity(intensity_label: str, ligand_count: int) -> float:
    """
    Convert a text intensity label to a number and scale by ligand count.

    The base scale (INTENSITY_SCALE) goes from 0.25 (weak) to 1.0 (very strong).
    Multiplying by ligand_count means a complex with 6 CN ligands will show
    proportionally stronger CN peaks than one with only 2.

    The result is NOT normalised here — the renderer handles normalisation
    so that the tallest peak always reaches 1.0 on the plot.

    Parameters
    ----------
    intensity_label : str
        e.g. "strong", "medium", "very strong"
    ligand_count : int
        How many of this ligand are in the complex (from the parser).

    Returns
    -------
    float
        Raw scaled intensity value.
    """
    base = INTENSITY_SCALE.get(intensity_label, 0.5)
    return base * ligand_count