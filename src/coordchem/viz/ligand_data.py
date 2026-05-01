"""
coordchem/viz/ligand_data.py
----------------------------
Ligand-specific drawing data for coordination complex 2D depictions.
"""

from __future__ import annotations


LIGAND_SMILES: dict[str, str] = {
    "NH3": "N",
    "H2O": "O",
    "Cl": "[Cl-]",
    "Br": "[Br-]",
    "I": "[I-]",
    "F": "[F-]",
    "OH": "[OH-]",
    "O": "[O-2]",
    "S": "[S-2]",
    "CN": "[C-]#N",
    "NC": "N#[C-]",
    "CO": "[C-]#[O+]",
    "NO": "[N]=O",
    "NO2": "O=[N+]([O-])",
    "ONO": "O[N+](=O)[O-]",
    "SCN": "[S-]C#N",
    "NCS": "N=C=[S-]",
    "N3": "[N-]=[N+]=N",
    "CN2H2": "NC#N",
    "en": "NCCN",
    "phen": "c1ccc2nc3ccccc3nc2c1",
    "bipy": "n1ccccc1-c1ncccc1",
    "bpy": "n1ccccc1-c1ncccc1",
    "ox": "[O-]C(=O)C(=O)[O-]",
    "acac": "CC(=O)C=C(O)C",
    "dien": "NCCNCCN",
    "tren": "N(CCN)CCN",
    "EDTA": "N(CC(=O)O)(CC(=O)O)CCN(CC(=O)O)CC(=O)O",
    "edta": "N(CC(=O)O)(CC(=O)O)CCN(CC(=O)O)CC(=O)O",
    "Cp": "c1cccc1",
    "tpy": "n1ccccc1-c1nc(ccc1)-c1ncccc1",
    "terpy": "n1ccccc1-c1nc(ccc1)-c1ncccc1",
    "py": "n1ccccc1",
    "dmso": "CS(=O)C",
    "PPh3": "P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "PMe3": "P(C)(C)C",
    "PEt3": "P(CC)(CC)CC",
}


LIGAND_DONOR_INDEX_OVERRIDES: dict[str, tuple[int, ...]] = {
    "NH3": (0,),
    "H2O": (0,),
    "OH": (0,),
    "CN": (0,),
    "NC": (0,),
    "CO": (0,),
    "NO": (0,),
    "NO2": (1,),
    "ONO": (0,),
    "SCN": (0,),
    "NCS": (0,),
    "py": (0,),
    "en": (0, 3),
    "bipy": (0, 7),
    "bpy": (0, 7),
    "phen": (4, 11),
    "ox": (0, 5),
    "acac": (2, 5),
    "dien": (0, 3, 6),
    "EDTA": (0, 4, 8, 11, 15, 19),
    "edta": (0, 4, 8, 11, 15, 19),
    "tpy": (0, 7, 13),
    "terpy": (0, 7, 13),
    "dmso": (1,),
    "PPh3": (0,),
    "PMe3": (0,),
    "PEt3": (0,),
}


EXPLICIT_H_LIGANDS: set[str] = set()


ABBREVIATED_MONODENTATE_LIGANDS: set[str] = {
    "NH3",
    "H2O",
    "Cl",
    "Br",
    "I",
    "F",
    "OH",
    "O",
    "S",
    "CN",
    "NC",
    "CO",
    "NO",
    "NO2",
    "ONO",
    "SCN",
    "NCS",
    "N3",
}

LIGAND_DISPLAY_LABELS: dict[str, tuple[str, str]] = {
    "NH3": ("NH3", "H3N"),
    "H2O": ("OH2", "H2O"),
    "OH": ("OH", "HO"),
    "en": ("NH2", "H2N"),
    "dien": ("NH2", "H2N"),
    "tren": ("NH2", "H2N"),
    "EDTA": ("NH2", "H2N"),
    "edta": ("NH2", "H2N"),
}

MONODENTATE_DISPLAY_LABELS: dict[str, tuple[str, str]] = {
    "NH3": ("NH3", "H3N"),
    "H2O": ("OH2", "H2O"),
    "OH": ("OH", "HO"),
    "CN": ("CN", "NC"),
    "NC": ("NC", "CN"),
    "CO": ("CO", "OC"),
    "NO": ("NO", "ON"),
    "NO2": ("NO2", "O2N"),
    "ONO": ("ONO", "ONO"),
    "SCN": ("SCN", "NCS"),
    "NCS": ("NCS", "SCN"),
    "N3": ("N3", "N3"),
}

POLYDENTATE_DONOR_DISPLAY_LABELS: dict[str, dict[int, tuple[str, str]]] = {
    "en": {
        0: ("NH2", "H2N"),
        3: ("NH2", "H2N"),
    },
    "dien": {
        0: ("NH2", "H2N"),
        3: ("NH", "HN"),
        6: ("NH2", "H2N"),
    },
    "tren": {
        0: ("N", "N"),
        3: ("NH2", "H2N"),
        6: ("NH2", "H2N"),
        9: ("NH2", "H2N"),
    },
}


INTERRUPTED_LIGAND_BONDS: dict[str, set[tuple[int, int]]] = {
    "EDTA": {
        (5, 6),
        (16, 17),
    },
    "edta": {
        (5, 6),
        (16, 17),
    },
}
