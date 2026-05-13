"""
coordchem.hsab
--------------

Simple Hard-Soft Acid-Base logic for ambidentate ligands.

This module is not intended to be a full quantitative HSAB model.
It provides simple rule-based decisions for ligands whose donor atom
can change depending on the metal centre.
"""

from __future__ import annotations


# Very simplified classification.
# This can be extended later.
HARD_METALS = {
    "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co",  # mostly when oxidation state is high
    "Al", "Ga",
    "Mg", "Ca", "Sr", "Ba",
    "Ln",
}

SOFT_METALS = {
    "Pd", "Pt", "Ag", "Au", "Hg", "Cd",
    "Cu",  # especially Cu(I)
}

BORDERLINE_METALS = {
    "Fe", "Co", "Ni", "Cu", "Zn", "Ru", "Rh", "Ir",
}


def classify_metal_hardness(metal: str, oxidation_state: int | None = None) -> str:
    """
    Return a simple HSAB classification for a metal centre.

    Returns
    -------
    str
        "hard", "soft" or "borderline".
    """

    if metal in {"Pd", "Pt", "Ag", "Au", "Hg", "Cd"}:
        return "soft"

    if metal == "Cu":
        if oxidation_state == 1:
            return "soft"
        return "borderline"

    if metal in {"Fe", "Co", "Cr", "Mn"}:
        if oxidation_state is not None and oxidation_state >= 3:
            return "hard"
        return "borderline"

    if metal in HARD_METALS:
        return "hard"

    if metal in SOFT_METALS:
        return "soft"

    return "borderline"


def choose_dmso_donor(metal: str, oxidation_state: int | None = None) -> tuple[str, str | None]:
    """
    Choose whether DMSO should bind through sulfur or oxygen.

    DMSO is an ambidentate ligand:
    - O-bound DMSO is favoured by hard metal centres.
    - S-bound DMSO is favoured by soft metal centres.
    - Borderline cases are chemically ambiguous.

    Returns
    -------
    tuple[str, str | None]
        donor atom, warning message
    """

    hardness = classify_metal_hardness(metal, oxidation_state)

    if hardness == "hard":
        return "O", (
            f"DMSO assigned as O-bound using HSAB theory: "
            f"{metal}({oxidation_state}) is treated as a hard metal centre."
        )

    if hardness == "soft":
        return "S", (
            f"DMSO assigned as S-bound using HSAB theory: "
            f"{metal}({oxidation_state}) is treated as a soft metal centre."
        )

    return "S/O", (
        f"DMSO binding mode is ambiguous for {metal}({oxidation_state}). "
        "The metal is treated as borderline by the simple HSAB model."
    )


def choose_ambidentate_donor(
    ligand: str,
    metal: str,
    oxidation_state: int | None = None,
) -> tuple[str | None, str | None]:
    """
    Choose the donor atom for ambidentate ligands.

    Currently implemented:
    - DMSO: S-bound or O-bound depending on HSAB logic.
    """

    key = ligand.lower()

    if key == "dmso":
        return choose_dmso_donor(metal, oxidation_state)

    return None, None