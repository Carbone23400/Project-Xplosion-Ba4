"""
coordchem/geometry.py
---------------------
Predict coordination geometry and d-electron count from a complex formula
or ParsedComplex object.
"""

from .parser import FormulaParseError, METALS, ParsedComplex, parse_formula


SQUARE_PLANAR_D8_METALS = {"Ni", "Pd", "Pt", "Rh", "Ir", "Au"}
COMMON_TETRAHEDRAL_METALS = {"Zn", "Cd", "Hg", "Co", "Ni", "Cu"}

METAL_GROUPS = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8,
    "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8,
    "Rh": 9, "Pd": 10, "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8,
    "Ir": 9, "Pt": 10, "Au": 11, "Hg": 12,
}


def get_geometry(complex_input: str | ParsedComplex) -> str:
    """Return the predicted coordination geometry."""
    parsed = _ensure_parsed_complex(complex_input)
    return predict_geometry(parsed)


def get_d_count(complex_input: str | ParsedComplex) -> int | None:
    """
    Return the d-electron count of the metal center.

    d-count = metal group number - oxidation state
    """
    parsed = _ensure_parsed_complex(complex_input)
    group_number = METAL_GROUPS.get(parsed.metal)

    if group_number is None or parsed.oxidation_state is None:
        return None

    return group_number - parsed.oxidation_state


def predict_geometry(parsed: ParsedComplex) -> str:
    """Predict geometry from coordination number, metal, oxidation state, and ligands."""
    cn = parsed.coordination_number
    metal = parsed.metal
    ox = parsed.oxidation_state
    ligands = parsed.ligands

    if cn in {1, 2}:
        return "linear"

    if cn == 3:
        return "trigonal planar"

    if cn == 4:
        return _predict_cn4_geometry(metal, ox, ligands)

    if cn == 5:
        return "trigonal bipyramidal or square pyramidal"

    if cn == 6:
        return "octahedral"

    if cn == 7:
        return "pentagonal bipyramidal or capped octahedral"

    if cn == 8:
        return "square antiprismatic or dodecahedral"

    return "unknown geometry"


def geometry_report(complex_input: str | ParsedComplex) -> dict[str, object]:
    """Return a report with parsed information, d-count, and predicted geometry."""
    parsed = _ensure_parsed_complex(complex_input)

    return {
        "formula": parsed.raw_formula,
        "metal": parsed.metal,
        "oxidation_state": parsed.oxidation_state,
        "d_count": get_d_count(parsed),
        "coordination_number": parsed.coordination_number,
        "ligands": parsed.ligands,
        "geometry": predict_geometry(parsed),
        "warnings": parsed.warnings,
        "errors": parsed.errors,
    }


def _ensure_parsed_complex(complex_input: str | ParsedComplex) -> ParsedComplex:
    """Convert a formula string, compound name, or ParsedComplex if needed."""
    if isinstance(complex_input, ParsedComplex):
        return complex_input

    if isinstance(complex_input, str):
        if _has_complex_brackets(complex_input):
            return parse_formula(complex_input)

        from .name2 import parse_name
        if not _looks_like_unbracketed_formula(complex_input):
            return parse_name(complex_input)

        try:
            return parse_formula(complex_input)
        except FormulaParseError:
            return parse_name(complex_input)

    raise TypeError("complex_input must be a formula string, compound name, or ParsedComplex object")


def _has_complex_brackets(complex_input: str) -> bool:
    """Return True when the input contains coordination-complex brackets."""
    return "[" in complex_input or "]" in complex_input


def _looks_like_unbracketed_formula(complex_input: str) -> bool:
    """Return True for formulas like PtCl2(NH3)2 without complex brackets."""
    text = complex_input.strip()
    if not text or text[0].islower() or " " in text:
        return False

    for length in (2, 1):
        metal = text[:length]
        remainder = text[length:]
        if metal in METALS and remainder:
            return remainder[0].isupper() or remainder[0] == "("

    return False


def _predict_cn4_geometry(
    metal: str,
    oxidation_state: int | None,
    ligands: dict[str, int],
) -> str:
    """Predict geometry for coordination number 4."""
    strong_field_ligands = {"CN", "CO", "NO", "phen", "bipy", "bpy", "PPh3"}

    if oxidation_state == 2 and metal in {"Pt", "Pd"}:
        return "square planar"

    if oxidation_state == 2 and metal == "Ni":
        if any(ligand in strong_field_ligands for ligand in ligands):
            return "square planar"
        return "tetrahedral or square planar"

    if oxidation_state == 2 and metal == "Cu":
        return "distorted square planar or tetrahedral"

    if metal in {"Zn", "Cd", "Hg"}:
        return "tetrahedral"

    return "tetrahedral or square planar"
