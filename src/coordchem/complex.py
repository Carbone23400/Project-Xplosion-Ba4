"""
coordchem/complex.py
--------------------
High-level user-facing API for coordination complexes.

This module provides the Complex class, which wraps a ParsedComplex object
and exposes convenient constructors, properties, and helper methods.
"""

from __future__ import annotations

from dataclasses import dataclass

from .geometry import get_d_count, get_geometry, geometry_report
from .name2 import parse_name
from .parser import ParsedComplex, parse_formula


@dataclass
class Complex:
    """High-level wrapper around a ParsedComplex object."""

    parsed: ParsedComplex

    @classmethod
    def from_formula(cls, formula: str) -> "Complex":
        """Build a Complex from a coordination formula string."""
        return cls(parse_formula(formula))

    @classmethod
    def from_name(cls, name: str) -> "Complex":
        """Build a Complex from a coordination compound name."""
        return cls(parse_name(name))

    @classmethod
    def from_input(cls, text: str) -> "Complex":
        """
        Build a Complex from either a formula or a name.

        Heuristic:
        - if brackets are present, treat as formula
        - otherwise try formula first, then fall back to name
        """
        if not isinstance(text, str):
            raise TypeError("text must be a string")

        text = text.strip()
        if not text:
            raise ValueError("text must not be empty")

        if "[" in text or "]" in text:
            return cls.from_formula(text)

        try:
            return cls.from_formula(text)
        except Exception:
            return cls.from_name(text)

    @property
    def metal(self) -> str:
        return self.parsed.metal

    @property
    def ligands(self) -> dict[str, int]:
        return self.parsed.ligands

    @property
    def ligand_names(self) -> dict[str, str]:
        return self.parsed.ligand_names

    @property
    def ligand_charges(self) -> dict[str, int]:
        return self.parsed.ligand_charges

    @property
    def ligand_denticity(self) -> dict[str, int]:
        return self.parsed.ligand_denticity

    @property
    def donor_atoms(self) -> dict[str, str]:
        return self.parsed.donor_atoms

    @property
    def oxidation_state(self) -> int | None:
        return self.parsed.oxidation_state

    @property
    def coordination_number(self) -> int:
        return self.parsed.coordination_number

    @property
    def complex_charge(self) -> int:
        return self.parsed.complex_charge

    @property
    def total_ligand_charge(self) -> int:
        return self.parsed.total_ligand_charge

    @property
    def counter_ions(self) -> dict[str, int]:
        return self.parsed.counter_ions

    @property
    def raw_formula(self) -> str:
        return self.parsed.raw_formula

    @property
    def warnings(self) -> list[str]:
        return self.parsed.warnings

    @property
    def errors(self) -> list[str]:
        return self.parsed.errors

    @property
    def geometry(self) -> str:
        """Predicted geometry."""
        return get_geometry(self.parsed)

    @property
    def d_count(self) -> int | None:
        """Metal d-electron count."""
        return get_d_count(self.parsed)

    def report(self) -> dict[str, object]:
        """Return a structured report for the complex."""
        return geometry_report(self.parsed)

    def draw_2d_svg(self, size: int = 700, title: str | None = None) -> str:
        """Return an SVG depiction of the complex."""
        from .viz.diagram_2d import diagram_2d_svg

        return diagram_2d_svg(self.parsed, size=size, title=title)

    def save_2d(
        self,
        output_path: str,
        size: int = 700,
        title: str | None = None,
    ):
        """Save a 2D SVG depiction to disk."""
        from .viz.diagram_2d import save_diagram_2d

        return save_diagram_2d(self.parsed, output_path=output_path, size=size, title=title)

    def build_3d(self, distance: float = 2.0):
        """Return an RDKit ``Mol`` with a 3D conformer for this complex."""
        from .viz.structure_3d import build_complex_3d

        return build_complex_3d(self.parsed, distance=distance)

    def draw_3d(self, width: int = 400, height: int = 400, distance: float = 2.0):
        """Return a ``py3Dmol.view`` displaying the 3D structure of the complex."""
        from .viz.structure_3d import view_complex_3d

        return view_complex_3d(
            self.parsed, width=width, height=height, distance=distance,
        )

    def draw_3d_html(
        self, width: int = 400, height: int = 400, distance: float = 2.0,
    ) -> str:
        """Return self-contained HTML for embedding the 3D view (e.g. in Streamlit)."""
        from .viz.structure_3d import complex_3d_html

        return complex_3d_html(
            self.parsed, width=width, height=height, distance=distance,
        )

    def predict_ir(self):
        """Return predicted IR spectrum object/data using the spectra module."""
        from .spectra.predictor import predict_ir

        return predict_ir(self.parsed)

    def predict_raman(self):
        """Return predicted Raman spectrum object/data using the spectra module."""
        from .spectra.predictor import predict_raman

        return predict_raman(self.parsed)

    def __str__(self) -> str:
        return str(self.parsed)