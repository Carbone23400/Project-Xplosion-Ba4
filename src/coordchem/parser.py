"""
coordchem/parser.py
-------------------
Parses coordination complex formulas like [Fe(CN)6]4- into structured data.

Supported formats:
    [Fe(CN)6]4-
    [Cu(NH3)4]2+
    [Co(en)3]3+
    [PtCl2(NH3)2]
    [Fe(H2O)6]3+
    K4[Fe(CN)6]          (with counter ions)
    [CoCl2(en)2]+
"""

import re
from dataclasses import dataclass, field
from typing import Optional
from .hsab import choose_ambidentate_donor


# ---------------------------------------------------------------------------
# Known ligands: formula → (name, charge, denticity, donor_atom)
# ---------------------------------------------------------------------------
KNOWN_LIGANDS: dict[str, tuple[str, int, int, str]] = {
    # formula        : (common_name,          charge, denticity, donor)
    "CN"            : ("cyano",                 -1,    1,      "C"),
    "NC"            : ("isocyano",              -1,    1,      "N"),
    "CO"            : ("carbonyl",                 0,    1,      "C"),
    "NO"            : ("nitrosyl",                 0,    1,      "N"),
    "Cl"            : ("chloro",                -1,    1,      "Cl"),
    "Br"            : ("bromo",                 -1,    1,      "Br"),
    "I"             : ("iodo",                  -1,    1,      "I"),
    "F"             : ("fluoro",                -1,    1,      "F"),
    "OH"            : ("hydroxo",               -1,    1,      "O"),
    "O"             : ("oxo",                   -2,    1,      "O"),
    "S"             : ("sulfido",                 -2,    1,      "S"),
    "NH3"           : ("ammine",                   0,    1,      "N"),
    "H2O"           : ("aqua",                     0,    1,      "O"),
    "NO2"           : ("nitro",                 -1,    1,      "N"),
    "ONO"           : ("nitrito",               -1,    1,      "O"),
    "SCN"           : ("thiocyano",           -1,    1,      "S"),
    "NCS"           : ("isothiocyano",           -1,    1,      "N"),
    "N3"            : ("azido",                   -1,    1,      "N"),
    "en"            : ("ethylenediamine",           0,    2,      "N"),
    "phen"          : ("1,10-phenanthroline",       0,    2,      "N"),
    "bipy"          : ("2,2'-bipyridine",           0,    2,      "N"),
    "bpy"           : ("2,2'-bipyridine",           0,    2,      "N"),
    "ox"            : ("oxalato",                 -2,    2,      "O"),
    "acac"          : ("acetylacetonato",          -1,    2,      "O"),
    "EDTA"          : ("ethylenediaminetetraacetato", -4, 6,     "N/O"),
    "edta"          : ("ethylenediaminetetraacetato", -4, 6,     "N/O"),
    "Cp"            : ("cyclopentadienyl",         -1,    1,      "C"),
    "tpy"           : ("terpyridine",               0,    3,      "N"),
    "terpy"         : ("terpyridine",               0,    3,      "N"),
    "py"            : ("pyridine",                  0,    1,      "N"),
    "dmso"          : ("dimethylsulfoxide",          0,    1,      "S"),
    "PPh3"          : ("triphenylphosphine",         0,    1,      "P"),
    "PMe3"          : ("trimethylphosphine",         0,    1,      "P"),
    "PEt3"          : ("triethylphosphine",          0,    1,      "P"),
    "CH3"           :("methyl",                     -1,     1,      "C"),
}

# Counter ions (outside the brackets)
COUNTER_IONS: dict[str, int] = {
    "K"   :  1, "Na" :  1, "Li" :  1, "Rb" :  1, "Cs" :  1,
    "Ca"  :  2, "Ba" :  2, "Mg" :  2, "Sr" :  2,
    "NH4" :  1,
    "Cl"  : -1, "Br" : -1, "I"  : -1, "F"  : -1,
    "NO3" : -1, "ClO4": -1, "BF4": -1, "PF6": -1,
    "SO4" : -2, "CO3": -2, "C2O4": -2,
    "PO4" : -3,
}

COUNTER_ION_NAMES: dict[str, str] = {
    "K"   : "potassium",
    "Na"  : "sodium",
    "Li"  : "lithium",
    "Rb"  : "rubidium",
    "Cs"  : "cesium",
    "Ca"  : "calcium",
    "Ba"  : "barium",
    "Mg"  : "magnesium",
    "Sr"  : "strontium",
    "NH4" : "ammonium",
    "Cl"  : "chloride",
    "Br"  : "bromide",
    "I"   : "iodide",
    "F"   : "fluoride",
    "NO3" : "nitrate",
    "ClO4": "perchlorate",
    "BF4" : "tetrafluoroborate",
    "PF6" : "hexafluorophosphate",
    "SO4" : "sulfate",
    "CO3" : "carbonate",
    "C2O4": "oxalate",
    "PO4" : "phosphate",
}

# All known metal symbols (subset most relevant to coordination chemistry)
METALS: set[str] = {
    "Li", "Be", "Na", "Mg", "Al", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr",
    "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu",
    "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Ac", "Th", "U",
}
                  

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ParsedComplex:
    """Structured representation of a coordination complex."""
    metal: str
    ligands: dict[str, int]              # {ligand_formula: count}
    complex_charge: int                  # charge on the complex ion
    counter_ions: dict[str, int]         # {ion: count} outside brackets
    raw_formula: str                     # original input string

    # Derived fields (populated by post-processing)
    iupac_name: Optional[str] = None  
    ligand_charges: dict[str, int] = field(default_factory=dict)   # per ligand
    ligand_names: dict[str, str]   = field(default_factory=dict)   # IUPAC names
    ligand_denticity: dict[str, int] = field(default_factory=dict) # bite number
    donor_atoms: dict[str, str]    = field(default_factory=dict)   # donor element
    coordination_number: int       = 0
    total_ligand_charge: int       = 0
    oxidation_state: Optional[int] = None
    errors: list[str]              = field(default_factory=list)
    warnings: list[str]            = field(default_factory=list)

    def __str__(self) -> str:
        charge_str = _format_charge(self.complex_charge)
        lines = [
            f"Formula          : {self.raw_formula}",
            f"Metal            : {self.metal}",
            f"Oxidation state  : {'+' if self.oxidation_state and self.oxidation_state > 0 else ''}"
                                f"{self.oxidation_state}",
            f"Coordination no. : {self.coordination_number}",
            f"Complex charge   : {charge_str}",
            f"Ligands          :",
        ]
        for lig, count in self.ligands.items():
            name    = self.ligand_names.get(lig, "unknown")
            charge  = self.ligand_charges.get(lig, 0)
            dent    = self.ligand_denticity.get(lig, 1)
            donor   = self.donor_atoms.get(lig, "?")
            lines.append(
                f"   {count}x {lig:<8} ({name:<30}) "
                f"charge={_format_charge(charge)}  "
                f"denticity={dent}  donor={donor}"
            )
        if self.counter_ions:
            lines.append(f"Counter ions     : {self.counter_ions}")
        if self.warnings:
            for w in self.warnings:
                lines.append(f"⚠  {w}")
        if self.errors:
            for e in self.errors:
                lines.append(f"✗  {e}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_formula(formula: str) -> ParsedComplex:
    """
    Parse a coordination complex formula string.

    Parameters
    ----------
    formula : str
        e.g. "[Fe(CN)6]4-", "K4[Fe(CN)6]", "[Co(en)3]3+", "[PtCl2(NH3)2]"

    Returns
    -------
    ParsedComplex
        Fully populated data object.

    Raises
    ------
    FormulaParseError
        If the formula is structurally invalid.
    """
    formula = formula.strip()
    raw     = formula

    # 1. Split counter ions from the complex bracket
    counter_ions, inner = _extract_counter_ions(formula)

    # 2. Extract complex charge from suffix notation, or infer it from
    # counter-ions when the charge is omitted.
    inner, complex_charge, has_explicit_charge = _extract_complex_charge(inner)
    if not has_explicit_charge and counter_ions:
        complex_charge = _infer_complex_charge_from_counter_ions(counter_ions)

    # 3. Strip outer brackets  [ ... ]
    inner = _strip_outer_brackets(inner)

    # 4. Extract metal symbol
    metal, remainder = _extract_metal(inner)

    # 5. Parse ligands from the remainder
    ligands = _parse_ligands(remainder)

    # 6. Build the result object
    result = ParsedComplex(
        metal        = metal,
        ligands      = ligands,
        complex_charge = complex_charge,
        counter_ions = counter_ions,
        raw_formula  = raw,
    )

    # 7. Enrich with ligand metadata + derived fields
    _enrich(result)

    return result


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _extract_counter_ions(formula: str) -> tuple[dict[str, int], str]:
    """
    Separate counter ions written outside the brackets.

    K4[Fe(CN)6]   →  counter_ions={"K": 4},  inner="[Fe(CN)6]"
    [Cu(NH3)4]2+  →  counter_ions={},         inner="[Cu(NH3)4]2+"
    """
    counter_ions: dict[str, int] = {}

    # Pattern: optional leading cations + bracket complex + optional trailing anions
    # Leading ions: e.g. K4, Na2, (NH4)2
    # Trailing ions: e.g. Cl2, SO4, (NO3)2

    bracket_start = formula.find("[")
    bracket_end   = formula.rfind("]")

    if bracket_start == -1:
        # No brackets — assume the whole thing is the complex (unusual input)
        return {}, formula

    # Before the bracket
    before = formula[:bracket_start]
    # After the bracket (may include charge AND counter ion, handle carefully)
    after_bracket = formula[bracket_end + 1:]

    # Parse leading counter ions
    if before:
        ions, _ = _parse_ion_string(before)
        counter_ions.update(ions)

    # The charge suffix on the complex comes first after ]
    # e.g. ]4- or ]2+ — we leave that for _extract_complex_charge
    # Trailing counter ions would be something like ]Cl2 (rare) — try to detect
    # For now we leave 'after_bracket' untouched; charge parser handles it.

    charge_suffix, trailing_ion_text = _split_charge_suffix_and_trailing_ions(
        after_bracket
    )
    if trailing_ion_text:
        ions, _ = _parse_ion_string(trailing_ion_text)
        for ion, count in ions.items():
            counter_ions[ion] = counter_ions.get(ion, 0) + count

    inner = formula[bracket_start: bracket_end + 1] + charge_suffix
    return counter_ions, inner


def _parse_ion_string(s: str) -> tuple[dict[str, int], str]:
    """Parse ion symbols with optional counts, e.g. 'K4', 'Cl3', '(NH4)2'."""
    ions: dict[str, int] = {}
    pos = 0
    s = s.strip()
    symbols = sorted(COUNTER_IONS.keys(), key=len, reverse=True)

    while pos < len(s):
        if s[pos].isspace():
            pos += 1
            continue

        symbol = None
        if s[pos] == "(":
            close = s.find(")", pos + 1)
            if close != -1:
                bracketed = s[pos + 1:close]
                if bracketed in COUNTER_IONS:
                    symbol = bracketed
                    pos = close + 1

        if symbol is None:
            for candidate in symbols:
                if s.startswith(candidate, pos):
                    symbol = candidate
                    pos += len(candidate)
                    break

        if symbol is None:
            return ions, s[pos:]

        count_m = re.match(r"\d+", s[pos:])
        count = int(count_m.group(0)) if count_m else 1
        if count_m:
            pos += len(count_m.group(0))
        ions[symbol] = ions.get(symbol, 0) + count

    return ions, ""


def _split_charge_suffix_and_trailing_ions(after_bracket: str) -> tuple[str, str]:
    """
    Split text after ']' into an optional charge suffix and trailing ions.

    Examples:
        "4-"    -> ("4-", "")
        "Cl3"   -> ("", "Cl3")
        "3+Cl3" -> ("3+", "Cl3")
    """
    if not after_bracket:
        return "", ""

    match = re.match(r"^(\d*[+-])(.*)$", after_bracket)
    if match:
        return match.group(1), match.group(2)

    return "", after_bracket


def _extract_complex_charge(formula: str) -> tuple[str, int, bool]:
    """
    Extract the charge suffix from a bracketed formula.

    "[Fe(CN)6]4-"  →  ("[Fe(CN)6]", -4)
    "[Cu(NH3)4]2+" →  ("[Cu(NH3)4]", +2)
    "[PtCl2(NH3)2]"→  ("[PtCl2(NH3)2]", 0)
    """
    # Match ] followed by optional digit and +/-
    pattern = re.compile(r'(\])(\d*)([+-]?)$')
    match   = pattern.search(formula)

    if not match:
        return formula, 0, False

    bracket_close = match.group(1)   # ]
    magnitude     = match.group(2)   # e.g. "4" or ""
    sign          = match.group(3)   # "+" or "-" or ""

    if not sign:
        charge = 0
    else:
        mag    = int(magnitude) if magnitude else 1
        charge = mag if sign == "+" else -mag

    # Remove the charge suffix from the formula string
    clean = formula[:match.start()] + bracket_close
    return clean, charge, bool(sign)


def _infer_complex_charge_from_counter_ions(counter_ions: dict[str, int]) -> int:
    """Infer the complex charge needed to neutralize the parsed counter-ions."""
    counter_charge = sum(
        COUNTER_IONS.get(ion, 0) * count
        for ion, count in counter_ions.items()
    )
    return -counter_charge


def _strip_outer_brackets(formula: str) -> str:
    """Remove the outermost [ and ] brackets."""
    formula = formula.strip()
    if formula.startswith("[") and formula.endswith("]"):
        return formula[1:-1]
    return formula


def _extract_metal(inner: str) -> tuple[str, str]:
    """
    Extract the metal symbol from the start of the inner formula.

    "Fe(CN)6"       →  ("Fe", "(CN)6")
    "CoCl2(NH3)4"   →  ("Co", "Cl2(NH3)4")
    "Pt(NH3)2Cl2"   →  ("Pt", "(NH3)2Cl2")

    Raises FormulaParseError if no metal found.
    """
    # Try two-letter metals first, then one-letter
    for length in (2, 1):
        candidate = inner[:length]
        if candidate in METALS:
            return candidate, inner[length:]

    raise FormulaParseError(
        f"Could not identify a metal symbol at the start of '{inner}'. "
        f"Make sure the metal comes first inside the brackets."
    )


def _parse_ligands(remainder: str) -> dict[str, int]:
    """
    Parse the ligand portion of an inner formula string.

    "(CN)6"              →  {"CN": 6}
    "Cl2(NH3)4"          →  {"Cl": 2, "NH3": 4}
    "(en)3"              →  {"en": 3}
    "Cl2(en)2"           →  {"Cl": 2, "en": 2}
    "(NH3)2(H2O)2Cl2"    →  {"NH3": 2, "H2O": 2, "Cl": 2}
    """
    ligands: dict[str, int] = {}
    s = remainder

    while s:
        # --- Bracketed ligand: (XX)n ---
        m = re.match(r'\(([^)]+)\)(\d*)', s)
        if m:
            lig_formula = m.group(1)
            count       = int(m.group(2)) if m.group(2) else 1
            ligands[lig_formula] = ligands.get(lig_formula, 0) + count
            s = s[m.end():]
            continue

        # --- Known multi-character ligand abbreviations (en, acac, etc.) ---
        matched_abbrev = None
        for abbrev in sorted(KNOWN_LIGANDS.keys(), key=len, reverse=True):
            if s.startswith(abbrev):
                # Make sure it's not the start of a longer token
                after = s[len(abbrev):]
                # Valid if followed by digit, (, [, end, or another uppercase
                if not after or re.match(r'[\d([\]A-Z]', after):
                    matched_abbrev = abbrev
                    break

        if matched_abbrev:
            after = s[len(matched_abbrev):]
            count_m = re.match(r'(\d+)', after)
            count   = int(count_m.group(1)) if count_m else 1
            skip    = len(matched_abbrev) + (len(count_m.group(1)) if count_m else 0)
            ligands[matched_abbrev] = ligands.get(matched_abbrev, 0) + count
            s = s[skip:]
            continue

        # --- Unbracketed ligand: parse greedily using known ligand table ---
        # Try longest match first
        found = False
        for lig in sorted(KNOWN_LIGANDS.keys(), key=len, reverse=True):
            if s.startswith(lig):
                after   = s[len(lig):]
                count_m = re.match(r'(\d+)', after)
                count   = int(count_m.group(1)) if count_m else 1
                skip    = len(lig) + (len(count_m.group(1)) if count_m else 0)
                ligands[lig] = ligands.get(lig, 0) + count
                s    = s[skip:]
                found = True
                break

        if not found:
            # Try to parse as a raw element symbol + count (fallback)
            m = re.match(r'([A-Z][a-z]?)(\d*)', s)
            if m:
                lig   = m.group(1)
                count = int(m.group(2)) if m.group(2) else 1
                ligands[lig] = ligands.get(lig, 0) + count
                s = s[m.end():]
            else:
                # Cannot parse — skip one character to avoid infinite loop
                s = s[1:]

    return ligands


# ---------------------------------------------------------------------------
# Enrichment — oxidation state, coordination number, ligand metadata
# ---------------------------------------------------------------------------

def _enrich(result: ParsedComplex) -> None:
    """Populate derived fields on a ParsedComplex."""

    total_ligand_charge = 0
    total_coord_number  = 0

    for lig_formula, count in result.ligands.items():
        info = KNOWN_LIGANDS.get(lig_formula)

        if info:
            name, charge, denticity, donor = info
            result.ligand_names[lig_formula]    = name
            result.ligand_charges[lig_formula]  = charge
            result.ligand_denticity[lig_formula]= denticity
            result.donor_atoms[lig_formula]     = donor

            total_ligand_charge += charge * count
            total_coord_number  += denticity * count
        else:
            result.warnings.append(
                f"Ligand '{lig_formula}' not found in the known ligands table. "
                f"Charge assumed 0, denticity assumed 1."
            )
            result.ligand_names[lig_formula]    = lig_formula
            result.ligand_charges[lig_formula]  = 0
            result.ligand_denticity[lig_formula]= 1
            result.donor_atoms[lig_formula]     = "?"
            total_coord_number += count

    result.total_ligand_charge  = total_ligand_charge
    result.coordination_number  = total_coord_number

    # Oxidation state: metal_charge + ligand_charge = complex_charge
    # → metal_charge = complex_charge − ligand_charge
    result.oxidation_state = result.complex_charge - total_ligand_charge

    _apply_ambidentate_donor_assignments(result)

    # input checks
    if result.coordination_number == 0:
        result.errors.append("Coordination number is 0 — no ligands were parsed.")

    if result.oxidation_state < -4 or result.oxidation_state > 9:
        result.warnings.append(
            f"Unusual oxidation state: {result.oxidation_state}. "
            f"Double-check the formula and charge."
        )
    
    # Generate IUPAC name if not already set (i.e. input was a formula)
    if result.iupac_name is None:
        try:
            result.iupac_name = get_iupac_name(result)
        except Exception:
            result.iupac_name = "unknown"


def _apply_ambidentate_donor_assignments(result: ParsedComplex) -> None:
    """Apply metal-dependent donor choices for ambidentate ligands."""
    result.warnings = [
        warning for warning in result.warnings
        if not warning.startswith("DMSO ")
    ]

    for lig_formula in result.ligands:
        donor, warning = choose_ambidentate_donor(
            lig_formula,
            result.metal,
            result.oxidation_state,
        )
        if donor is not None:
            result.donor_atoms[lig_formula] = donor
        if warning is not None:
            result.warnings.append(warning)

def get_iupac_name(parsed) -> str:
    """
    Build a basic IUPAC name from a ParsedComplex object.
    e.g. [Fe(CN)6]4- → hexacyanidoferrate(II)
    """
    # Number prefixes
    PREFIXES = {
        1: "mono", 2: "di",    3: "tri",   4: "tetra",
        5: "penta", 6: "hexa", 7: "hepta", 8: "octa",
        9: "nona",  10: "deca"
    }

    # Build ligand part — alphabetical order by ligand name (IUPAC rule)
    ligand_parts = []
    for lig, count in sorted(
        parsed.ligands.items(),
        key=lambda x: parsed.ligand_names.get(x[0], x[0])
    ):
        name   = parsed.ligand_names.get(lig, lig)
        if count == 1:
            ligand_parts.append(name)
            continue

        prefix = PREFIXES.get(count, str(count))
        # Use bis/tris/tetrakis for complex ligand names containing a number
        if any(p in name for p in ("di", "tri", "tetra", "bis", "tris")):
            complex_prefixes = {2: "bis", 3: "tris", 4: "tetrakis",
                                5: "pentakis", 6: "hexakis"}
            prefix = complex_prefixes.get(count, str(count))
            ligand_parts.append(f"{prefix}({name})")
        else:
            ligand_parts.append(f"{prefix}{name}")

    ligand_str = "".join(ligand_parts)

    # Metal name — lowercase
    ANION_METAL_NAMES = {
        "Fe": "ferrate", "Cu": "cuprate", "Co": "cobaltate",
        "Ni": "nickelate", "Cr": "chromate", "Mn": "manganate",
        "Pt": "platinate", "Pd": "palladate", "Au": "aurate",
        "Ag": "argentate", "Zn": "zincate",  "Mo": "molybdate",
        "W":  "tungstate",  "Ru": "ruthenate","Os": "osmate",
        "Rh": "rhodate",    "Ir": "iridate",  "V":  "vanadate",
        "Ti": "titanate",   "Cd": "cadmate",   "Zr": "zirconate",
    }

    CATION_METAL_NAMES = {
        "Fe": "iron", "Cu": "copper", "Co": "cobalt",
        "Ni": "nickel", "Cr": "chromium", "Mn": "manganese",
        "Pt": "platinum", "Pd": "palladium", "Au": "gold",
        "Ag": "silver", "Zn": "zinc", "Mo": "molybdenum",
        "W": "tungsten", "Ru": "ruthenium", "Os": "osmium",
        "Rh": "rhodium", "Ir": "iridium", "V": "vanadium",
        "Ti": "titanium", "Cd": "cadmium", "Zr": "zirconium",
    }

    ox = parsed.oxidation_state
    ox_str     = f"({_roman(ox)})" if ox is not None else ""
    metal_name = ANION_METAL_NAMES.get(parsed.metal, parsed.metal.lower() + "ate")

    # If complex charge is positive or zero, use metal name directly (cation)
    if parsed.complex_charge >= 0:
        metal_name = CATION_METAL_NAMES.get(parsed.metal, parsed.metal.lower())

    complex_name = f"{ligand_str}{metal_name}{ox_str}"
    if not parsed.counter_ions:
        return complex_name

    leading_counter_ions = []
    trailing_counter_ions = []
    for ion in parsed.counter_ions:
        ion_name = COUNTER_ION_NAMES.get(ion, ion.lower())
        if COUNTER_IONS.get(ion, 0) > 0:
            leading_counter_ions.append(ion_name)
        else:
            trailing_counter_ions.append(ion_name)

    return " ".join(leading_counter_ions + [complex_name] + trailing_counter_ions)


def _roman(n: int) -> str:
    """Convert an integer to a Roman numeral string."""
    if n is None:
        return ""
    if n == 0:
        return "0"
    vals  = [10,"X", 9,"IX", 5,"V", 4,"IV", 1,"I"]
    vals  = [(10,"X"),(9,"IX"),(5,"V"),(4,"IV"),(1,"I")]
    neg   = n < 0
    n     = abs(n)
    result = ""
    for value, numeral in vals:
        while n >= value:
            result += numeral
            n      -= value
    return f"−{result}" if neg else result

# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _format_charge(charge: int) -> str:
    if charge == 0:
        return "0"
    return f"+{charge}" if charge > 0 else str(charge)


# ---------------------------------------------------------------------------
# Custom exception
# ---------------------------------------------------------------------------

class FormulaParseError(ValueError):
    """Raised when a formula string cannot be parsed."""
    pass

