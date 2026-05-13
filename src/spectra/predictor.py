"""Compatibility wrapper for the historical src.spectra.predictor path."""

from coordchem.spectra.predictor import (
    INTENSITY_SCALE,
    PredictionResult,
    _scale_intensity,
    predict_ir,
    predict_raman,
    predict_spectrum,
)

__all__ = [
    "INTENSITY_SCALE",
    "PredictionResult",
    "_scale_intensity",
    "predict_ir",
    "predict_raman",
    "predict_spectrum",
]
