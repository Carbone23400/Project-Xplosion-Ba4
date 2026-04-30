"""
coordchem/viz/diagram_2d.py
---------------------------
Advanced 2D drawing for coordination complexes using RDKit.

This module parses a coordination complex with the existing package logic,
builds a composite RDKit molecule manually, assigns 2D coordinates, and draws
the result as SVG. It is designed for teaching/report depictions rather than
crystallographic accuracy.
"""

from __future__ import annotations

from html import escape
from importlib import import_module
from math import atan2, cos, pi, sin, sqrt
from pathlib import Path
from typing import Callable
from dataclasses import dataclass

from ..geometry import get_geometry
from ..parser import FormulaParseError, ParsedComplex, parse_formula


try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "coordchem.viz.diagram_2d requires RDKit. Install rdkit first."
    ) from exc


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
    "acac": "CC(=O)CC(=O)C",
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
    "acac": (1, 4),
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
    # (donor atom visually first, donor atom visually last)
    "H2O": ("OH₂", "H₂O"),
    "NH3": ("NH₃", "H₃N"),
    "OH": ("O-H", "HO-"),
}


MONODENTATE_DISPLAY_LABELS: dict[str, tuple[str, str]] = {
    # (donor atom visually first, donor atom visually last)
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


@dataclass(frozen=True)
class Site:
    """One coordination site around the metal."""

    x: float
    y: float
    style: str  # plain, wedge, dash


def _get_name_parser() -> Callable[[str], ParsedComplex] | None:
    """Find parse_name() in the evolving package structure."""
    for module_name in ("..name2", "..name", "..geometry"):
        try:
            module = import_module(module_name, package=__package__)
        except Exception:
            continue

        parse_name = getattr(module, "parse_name", None)
        if callable(parse_name):
            return parse_name

    return None


def parse_complex_input(complex_input: str | ParsedComplex) -> ParsedComplex:
    """Accept a formula string, compound name, or ParsedComplex."""
    if isinstance(complex_input, ParsedComplex):
        return complex_input

    if not isinstance(complex_input, str):
        raise TypeError("complex_input must be a string or ParsedComplex")

    try:
        return parse_formula(complex_input)
    except FormulaParseError:
        parse_name = _get_name_parser()
        if parse_name is None:
            raise ValueError(
                "Could not parse the input as a formula, and no parse_name() "
                "function could be located in the package."
            )
        return parse_name(complex_input)


def _expand_ligands(parsed: ParsedComplex) -> list[str]:
    """Expand {'CN': 6} into ['CN', 'CN', ...]."""
    ligands: list[str] = []
    for ligand, count in parsed.ligands.items():
        ligands.extend([ligand] * count)
    return ligands


def _coordination_sites(geometry: str, n_sites: int) -> list[Site]:
    """Return idealized 2D positions plus bond style cues."""
    g = geometry.lower().strip()

    if g == "linear":
        return [Site(-3.0, 0.0, "plain"), Site(3.0, 0.0, "plain")][:n_sites]

    if g == "trigonal planar":
        return [
            Site(0.0, -3.1, "plain"),
            Site(2.7, 1.55, "plain"),
            Site(-2.7, 1.55, "plain"),
        ][:n_sites]

    if g == "square planar":
        return [
            Site(0.0, -3.2, "plain"),
            Site(3.2, 0.0, "plain"),
            Site(0.0, 3.2, "plain"),
            Site(-3.2, 0.0, "plain"),
        ][:n_sites]

    if g == "tetrahedral":
        return [
            Site(-3.1, 0.2, "plain"),
            Site(3.1, 0.2, "plain"),
            Site(0.0, -2.8, "wedge"),
            Site(0.0, 2.9, "dash"),
        ][:n_sites]

    if g == "octahedral":
        return [
            Site(0.0, -3.4, "plain"),
            Site(3.4, 0.0, "plain"),
            Site(0.0, 3.4, "plain"),
            Site(-3.4, 0.0, "plain"),
            Site(2.4, -2.4, "wedge"),
            Site(-2.4, 2.4, "dash"),
        ][:n_sites]

    if "square planar" in g:
        return _coordination_sites("square planar", n_sites)

    if "tetrahedral" in g:
        return _coordination_sites("tetrahedral", n_sites)

    if "trigonal bipyramidal" in g:
        return [
            Site(0.0, -3.5, "plain"),
            Site(2.8, -0.7, "plain"),
            Site(1.7, 2.3, "plain"),
            Site(-1.7, 2.3, "plain"),
            Site(-2.8, -0.7, "plain"),
        ][:n_sites]

    if "square pyramidal" in g:
        return [
            Site(0.0, -3.2, "plain"),
            Site(3.2, 0.0, "plain"),
            Site(0.0, 3.2, "plain"),
            Site(-3.2, 0.0, "plain"),
            Site(2.2, -2.2, "wedge"),
        ][:n_sites]

    return _regular_polygon_sites(n_sites)


def _chelate_octahedral_sites(n_ligands: int) -> list[Site]:
    """Return spread-out site pairs for tris-bidentate octahedral complexes."""
    if n_ligands != 3:
        return _coordination_sites("octahedral", n_ligands * 2)

    return [
        Site(-3.2, -1.8, "dash"),
        Site(-1.2, -3.5, "plain"),
        Site(3.2, -1.8, "plain"),
        Site(3.2, 1.8, "wedge"),
        Site(-1.2, 3.5, "plain"),
        Site(-3.2, 1.8, "plain"),
    ]


def _tridentate_octahedral_sites(n_ligands: int) -> list[Site]:
    """Return site triplets for bis-tridentate octahedral complexes."""
    if n_ligands != 2:
        return _coordination_sites("octahedral", n_ligands * 3)

    radius = 3.4
    angles = [150, 90, 30, -30, -90, -150]
    styles = ["dash", "plain", "plain", "wedge", "plain", "plain"]

    return [
        Site(
            radius * cos(angle * pi / 180),
            radius * sin(angle * pi / 180),
            style,
        )
        for angle, style in zip(angles, styles)
    ]


def _edta_octahedral_sites() -> list[Site]:
    """Return a compact octahedral-like donor layout for EDTA."""
    return [
        Site(-1.45, 1.05, "dash"),
        Site(0.85, 1.55, "plain"),
        Site(-1.05, -0.95, "plain"),
        Site(-1.45, -1.05, "plain"),
        Site(0.85, -1.55, "wedge"),
        Site(1.8, 0.0, "plain"),
    ]


def _should_use_chelate_layout(parsed: ParsedComplex, ligand_items: list[str], geometry: str) -> bool:
    """Return True when a tris-bidentate octahedral layout will be clearer."""
    if geometry.lower().strip() != "octahedral":
        return False

    if len(ligand_items) != 3:
        return False

    return all(parsed.ligand_denticity.get(ligand, 1) == 2 for ligand in ligand_items)


def _should_use_tridentate_layout(parsed: ParsedComplex, ligand_items: list[str], geometry: str) -> bool:
    """Return True when a bis-tridentate octahedral layout will be clearer."""
    if geometry.lower().strip() != "octahedral":
        return False

    if len(ligand_items) != 2:
        return False

    return all(parsed.ligand_denticity.get(ligand, 1) == 3 for ligand in ligand_items)


def _should_use_edta_layout(
    parsed: ParsedComplex,
    ligand_items: list[str],
    geometry: str,
) -> bool:
    """Return True for the special one-ligand EDTA octahedral drawing."""
    if geometry.lower().strip() != "octahedral":
        return False

    if len(ligand_items) != 1:
        return False

    ligand = ligand_items[0]
    return ligand in {"EDTA", "edta"} and parsed.ligand_denticity.get(ligand) == 6


def _regular_polygon_sites(n_sites: int, radius: float = 3.2) -> list[Site]:
    """Fallback placement on a regular polygon."""
    if n_sites <= 0:
        return []

    sites: list[Site] = []
    for i in range(n_sites):
        angle = -pi / 2 + 2 * pi * i / n_sites
        sites.append(Site(radius * cos(angle), radius * sin(angle), "plain"))
    return sites


def _dist(x: float, y: float) -> float:
    return sqrt(x * x + y * y)


def _angle(x: float, y: float) -> float:
    return atan2(y, x)


def _max_distance_from_point(
    coords: dict[int, tuple[float, float]],
    center: tuple[float, float],
) -> float:
    """Return the maximum distance of any coordinate from center."""
    if not coords:
        return 1.0

    return max(_dist(x - center[0], y - center[1]) for x, y in coords.values())


def _rotate_point(x: float, y: float, theta: float) -> tuple[float, float]:
    c = cos(theta)
    s = sin(theta)
    return (c * x - s * y, s * x + c * y)


def _reflect_point_across_line(
    point: tuple[float, float],
    line_a: tuple[float, float],
    line_b: tuple[float, float],
) -> tuple[float, float]:
    """Reflect one point across the line passing through line_a and line_b."""
    px, py = point
    ax, ay = line_a
    bx, by = line_b
    vx = bx - ax
    vy = by - ay
    length_sq = vx * vx + vy * vy
    if length_sq == 0:
        return point

    t = ((px - ax) * vx + (py - ay) * vy) / length_sq
    proj_x = ax + t * vx
    proj_y = ay + t * vy
    return (2 * proj_x - px, 2 * proj_y - py)


def _make_ligand_mol(ligand_symbol: str) -> Chem.Mol:
    """Build one ligand molecule with 2D coordinates."""
    smiles = LIGAND_SMILES.get(ligand_symbol)
    if smiles is None:
        raise KeyError(
            f"No ligand SMILES registered for '{ligand_symbol}'. "
            "Add it to LIGAND_SMILES."
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid ligand SMILES for '{ligand_symbol}': {smiles}")

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        pass

    if ligand_symbol in EXPLICIT_H_LIGANDS:
        mol = Chem.AddHs(mol)

    rdDepictor.Compute2DCoords(mol)
    return mol


def _copy_atom(atom: Chem.Atom, atom_label: str | None = None) -> Chem.Atom:
    """Copy atom attributes needed for drawing."""
    new_atom = Chem.Atom(atom.GetAtomicNum())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
    new_atom.SetNoImplicit(atom.GetNoImplicit())
    new_atom.SetIsAromatic(atom.GetIsAromatic())
    if atom_label is not None:
        new_atom.SetProp("atomLabel", atom_label)
    return new_atom


def _display_atom_label(
    ligand_symbol: str,
    atom_idx: int,
    donor_indices: tuple[int, ...],
    donor_anchors: dict[int, Site],
) -> str | None:
    """Return a custom visible label for donor atoms that RDKit would simplify."""
    if atom_idx not in donor_indices:
        return None

    labels = LIGAND_DISPLAY_LABELS.get(ligand_symbol)
    if labels is None:
        return None

    donor_first, donor_last = labels
    anchor = donor_anchors.get(atom_idx)
    if anchor is not None and anchor.x < 0:
        return donor_last

    return donor_first


def _monodentate_label(
    ligand_symbol: str,
    anchor: Site,
) -> str:
    """Return a formula label with the donor side facing the metal."""
    donor_first, donor_last = MONODENTATE_DISPLAY_LABELS.get(
        ligand_symbol,
        (ligand_symbol, ligand_symbol),
    )

    if anchor.x < 0:
        return donor_last

    return donor_first


def _monodentate_label_atom(
    ligand_symbol: str,
    parsed: ParsedComplex,
    anchor: Site,
) -> Chem.Atom:
    """Build a single labeled donor atom for a monodentate ligand."""
    donor_symbol = str(parsed.donor_atoms.get(ligand_symbol, "")).split("/")[0]
    if not donor_symbol or donor_symbol == "?":
        donor_symbol = LIGAND_SMILES.get(ligand_symbol, ligand_symbol)
        donor_symbol = donor_symbol.strip("[]+-0123456789")

    if donor_symbol not in {"B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}:
        donor_symbol = "C"

    atom = Chem.Atom(donor_symbol)
    atom.SetNoImplicit(True)
    atom.SetProp("atomLabel", _monodentate_label(ligand_symbol, anchor))
    return atom


def _get_2d_coords(mol: Chem.Mol) -> dict[int, tuple[float, float]]:
    """Return atom index to 2D coordinate mapping."""
    conf = mol.GetConformer()
    return {
        i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    }


def _centroid(
    coords: dict[int, tuple[float, float]],
    exclude: set[int] | None = None,
) -> tuple[float, float]:
    """Centroid of all coordinates except optional excluded atoms."""
    exclude = exclude or set()
    points = [xy for i, xy in coords.items() if i not in exclude]
    if not points:
        return (0.0, 0.0)

    return (
        sum(point[0] for point in points) / len(points),
        sum(point[1] for point in points) / len(points),
    )


def _match_donor_indices(
    parsed: ParsedComplex,
    ligand_symbol: str,
    mol: Chem.Mol,
) -> tuple[int, ...]:
    """Determine donor atom indices for a ligand molecule."""
    override = LIGAND_DONOR_INDEX_OVERRIDES.get(ligand_symbol)
    if override is not None:
        return override

    donor_info = parsed.donor_atoms.get(ligand_symbol, "?")
    if donor_info in {"?", "", None}:
        raise ValueError(f"Could not determine donor atom for '{ligand_symbol}'.")

    candidate_symbols = tuple(
        part.strip() for part in str(donor_info).split("/") if part.strip()
    )

    matches = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetSymbol() in candidate_symbols
    ]
    if not matches:
        raise ValueError(
            f"No donor atom {candidate_symbols} found in '{ligand_symbol}'."
        )

    denticity = parsed.ligand_denticity.get(ligand_symbol, 1)
    if denticity == 1:
        matches.sort(key=lambda idx: mol.GetAtomWithIdx(idx).GetDegree())
        return (matches[0],)

    if len(matches) < denticity:
        raise ValueError(
            f"Ligand '{ligand_symbol}' has denticity {denticity}, "
            f"but only {len(matches)} donor atoms were found."
        )

    return tuple(matches[:denticity])


def _transform_monodentate(
    coords: dict[int, tuple[float, float]],
    donor_idx: int,
    anchor: Site,
) -> dict[int, tuple[float, float]]:
    """Place a monodentate ligand with donor atom on the coordination site."""
    donor_x, donor_y = coords[donor_idx]
    cx, cy = _centroid(coords, exclude={donor_idx})

    current_vx = cx - donor_x
    current_vy = cy - donor_y
    if _dist(current_vx, current_vy) < 1e-8:
        current_vx, current_vy = 1.0, 0.0

    theta = _angle(anchor.x, anchor.y) - _angle(current_vx, current_vy)

    new_coords: dict[int, tuple[float, float]] = {}
    for idx, (x, y) in coords.items():
        rx, ry = _rotate_point(x - donor_x, y - donor_y, theta)
        new_coords[idx] = (rx + anchor.x, ry + anchor.y)

    return new_coords


def _transform_polydentate(
    coords: dict[int, tuple[float, float]],
    donor_indices: tuple[int, ...],
    anchors: tuple[Site, ...],
) -> dict[int, tuple[float, float]]:
    """Place a polydentate ligand on multiple coordination sites."""
    if len(donor_indices) != len(anchors):
        raise ValueError("Number of donor indices and anchors must match.")

    if len(donor_indices) == 1:
        return _transform_monodentate(coords, donor_indices[0], anchors[0])

    d1 = donor_indices[0]
    d2 = donor_indices[-1]
    a1 = anchors[0]
    a2 = anchors[-1]

    x1, y1 = coords[d1]
    x2, y2 = coords[d2]

    cur_mid_x = (x1 + x2) / 2
    cur_mid_y = (y1 + y2) / 2
    cur_vx = x2 - x1
    cur_vy = y2 - y1
    cur_len = _dist(cur_vx, cur_vy) or 1.0

    tar_mid_x = (a1.x + a2.x) / 2
    tar_mid_y = (a1.y + a2.y) / 2
    tar_vx = a2.x - a1.x
    tar_vy = a2.y - a1.y
    tar_len = _dist(tar_vx, tar_vy) or 1.0

    theta = _angle(tar_vx, tar_vy) - _angle(cur_vx, cur_vy)
    scale = tar_len / cur_len

    new_coords: dict[int, tuple[float, float]] = {}
    for idx, (x, y) in coords.items():
        x0 = (x - cur_mid_x) * scale
        y0 = (y - cur_mid_y) * scale
        xr, yr = _rotate_point(x0, y0, theta)
        new_coords[idx] = (xr + tar_mid_x, yr + tar_mid_y)

    ligand_radius = _max_distance_from_point(new_coords, (tar_mid_x, tar_mid_y))
    target_radius = 3.2 if len(donor_indices) >= 3 else 2.8
    if len(donor_indices) < 3 and ligand_radius > target_radius:
        shrink = target_radius / ligand_radius
        new_coords = {
            idx: (
                tar_mid_x + (x - tar_mid_x) * shrink,
                tar_mid_y + (y - tar_mid_y) * shrink,
            )
            for idx, (x, y) in new_coords.items()
        }

    centroid = _centroid(new_coords, exclude=set(donor_indices))
    outward_x = tar_mid_x
    outward_y = tar_mid_y
    ligand_side_x = centroid[0] - tar_mid_x
    ligand_side_y = centroid[1] - tar_mid_y

    if ligand_side_x * outward_x + ligand_side_y * outward_y < 0:
        line_a = (a1.x, a1.y)
        line_b = (a2.x, a2.y)
        new_coords = {
            idx: _reflect_point_across_line(point, line_a, line_b)
            for idx, point in new_coords.items()
        }

    return new_coords


def _transform_edta(anchors: tuple[Site, ...]) -> dict[int, tuple[float, float]]:
    """Place EDTA in a wrapped hexadentate drawing around the metal."""
    if len(anchors) != 6:
        raise ValueError("EDTA layout requires six coordination sites.")

    n_top, o_top, o_left_bottom, n_bottom, o_bottom, o_right = anchors

    return {
        0: (n_top.x, n_top.y),
        1: (-0.95, 1.95),
        2: (0.05, 2.20),
        3: (0.10, 3.05),
        4: (o_top.x, o_top.y),
        5: (-2.25, 0.55),
        6: (-2.10, -0.45),
        7: (-3.00, -0.75),
        8: (o_left_bottom.x, o_left_bottom.y),
        9: (-2.35, 1.35),
        10: (-2.35, -1.35),
        11: (n_bottom.x, n_bottom.y),
        12: (-0.95, -1.95),
        13: (0.05, -2.20),
        14: (0.10, -3.05),
        15: (o_bottom.x, o_bottom.y),
        16: (0.55, -0.80),
        17: (1.35, -0.30),
        18: (2.10, -0.85),
        19: (o_right.x, o_right.y),
    }


def build_coordination_mol(
    complex_input: str | ParsedComplex,
    geometry_override: str | None = None,
) -> Chem.Mol:
    """Build a composite RDKit molecule for the complex."""
    parsed = parse_complex_input(complex_input)
    geometry = geometry_override or get_geometry(parsed)
    ligand_items = _expand_ligands(parsed)
    n_sites = sum(parsed.ligand_denticity.get(lig, 1) for lig in ligand_items)
    if _should_use_edta_layout(parsed, ligand_items, geometry):
        sites = _edta_octahedral_sites()
    elif _should_use_tridentate_layout(parsed, ligand_items, geometry):
        sites = _tridentate_octahedral_sites(len(ligand_items))
    elif _should_use_chelate_layout(parsed, ligand_items, geometry):
        sites = _chelate_octahedral_sites(len(ligand_items))
    else:
        sites = _coordination_sites(geometry, n_sites)

    rw = Chem.RWMol()
    metal_atom = Chem.Atom(parsed.metal)
    metal_atom.SetNoImplicit(True)
    metal_idx = rw.AddAtom(metal_atom)

    global_coords: dict[int, tuple[float, float]] = {metal_idx: (0.0, 0.0)}
    site_cursor = 0

    for ligand_symbol in ligand_items:
        denticity_from_parser = parsed.ligand_denticity.get(ligand_symbol, 1)
        if (
            denticity_from_parser == 1
            and ligand_symbol in ABBREVIATED_MONODENTATE_LIGANDS
        ):
            anchors = tuple(sites[site_cursor: site_cursor + 1])
            if len(anchors) != 1:
                raise ValueError(f"Not enough coordination sites for '{ligand_symbol}'.")
            site_cursor += 1

            anchor = anchors[0]
            ligand_idx = rw.AddAtom(
                _monodentate_label_atom(ligand_symbol, parsed, anchor)
            )
            rw.AddBond(metal_idx, ligand_idx, Chem.BondType.SINGLE)

            bond = rw.GetBondBetweenAtoms(metal_idx, ligand_idx)
            if anchor.style == "wedge":
                bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif anchor.style == "dash":
                bond.SetBondDir(Chem.BondDir.BEGINDASH)

            global_coords[ligand_idx] = (anchor.x, anchor.y)
            continue

        lig_mol = _make_ligand_mol(ligand_symbol)
        lig_coords = _get_2d_coords(lig_mol)
        donor_indices = _match_donor_indices(parsed, ligand_symbol, lig_mol)
        denticity = len(donor_indices)

        anchors = tuple(sites[site_cursor: site_cursor + denticity])
        if len(anchors) != denticity:
            raise ValueError(f"Not enough coordination sites for '{ligand_symbol}'.")
        site_cursor += denticity

        if ligand_symbol in {"EDTA", "edta"} and denticity == 6:
            transformed = _transform_edta(anchors)
        else:
            transformed = _transform_polydentate(lig_coords, donor_indices, anchors)
        donor_anchors = dict(zip(donor_indices, anchors))

        idx_map: dict[int, int] = {}
        for atom in lig_mol.GetAtoms():
            atom_label = _display_atom_label(
                ligand_symbol,
                atom.GetIdx(),
                donor_indices,
                donor_anchors,
            )
            idx_map[atom.GetIdx()] = rw.AddAtom(_copy_atom(atom, atom_label))

        for bond in lig_mol.GetBonds():
            rw.AddBond(
                idx_map[bond.GetBeginAtomIdx()],
                idx_map[bond.GetEndAtomIdx()],
                bond.GetBondType(),
            )

        for donor_local_idx, anchor in zip(donor_indices, anchors):
            donor_global_idx = idx_map[donor_local_idx]
            rw.AddBond(metal_idx, donor_global_idx, Chem.BondType.SINGLE)

            bond = rw.GetBondBetweenAtoms(metal_idx, donor_global_idx)
            if anchor.style == "wedge":
                bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif anchor.style == "dash":
                bond.SetBondDir(Chem.BondDir.BEGINDASH)

        for local_idx, global_idx in idx_map.items():
            global_coords[global_idx] = transformed[local_idx]

    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)

    conf = Chem.Conformer(mol.GetNumAtoms())
    for idx in range(mol.GetNumAtoms()):
        x, y = global_coords[idx]
        conf.SetAtomPosition(idx, (x, y, 0.0))

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def _geometry_options(geometry: str) -> list[str]:
    """Split ambiguous geometry labels into drawable alternatives."""
    if " or " not in geometry:
        return [geometry]

    options: list[str] = []
    for option in geometry.split(" or "):
        clean = option.strip()
        if clean == "distorted square planar":
            clean = "square planar"
        if clean not in options:
            options.append(clean)

    return options or [geometry]


def _add_svg_labels(
    svg: str,
    size: int,
    title: str | None,
) -> str:
    """Add a full white background and an optional centered title."""
    svg = svg.replace(
        "<!-- END OF HEADER -->",
        (
            "<!-- END OF HEADER -->\n"
            f"<rect width='{size}' height='{size}' fill='#FFFFFF'/>\n"
        ),
        1,
    )

    if title:
        title_markup = (
            f"<text x='{size / 2:.1f}' y='28' text-anchor='middle' "
            "font-family='Arial, Helvetica, sans-serif' font-size='18' "
            "font-weight='700' fill='#111111'>"
            f"{escape(title)}</text>\n"
        )
        svg = svg.replace("</svg>", f"{title_markup}</svg>", 1)

    return svg


def diagram_2d_svg(
    complex_input: str | ParsedComplex,
    size: int = 700,
    *,
    title: str | None = None,
) -> str:
    """Draw the coordination complex as SVG using RDKit."""
    if size <= 0:
        raise ValueError("size must be a positive integer")

    parsed = parse_complex_input(complex_input)
    geometry = get_geometry(parsed)
    geometry_options = _geometry_options(geometry)
    mols = [
        build_coordination_mol(parsed, geometry_override=option)
        for option in geometry_options
    ]

    if len(mols) == 1:
        drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    else:
        panel_width = size // len(mols)
        drawer = rdMolDraw2D.MolDraw2DSVG(size, size, panel_width, size)
    options = drawer.drawOptions()
    options.useBWAtomPalette()
    options.clearBackground = True
    options.bondLineWidth = 2.2
    options.fixedBondLength = 32
    options.addStereoAnnotation = False
    options.prepareMolsBeforeDrawing = False
    options.legendFontSize = 18

    if len(mols) == 1:
        drawer.DrawMolecule(mols[0], legend=geometry, confId=0)
    else:
        legends = geometry_options
        drawer.DrawMolecules(mols, legends=legends, confIds=[0] * len(mols))

    drawer.FinishDrawing()
    return _add_svg_labels(drawer.GetDrawingText(), size=size, title=title)


def save_diagram_2d(
    complex_input: str | ParsedComplex,
    output_path: str | Path,
    size: int = 700,
    *,
    title: str | None = None,
) -> Path:
    """Save the SVG depiction to disk."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    svg = diagram_2d_svg(complex_input=complex_input, size=size, title=title)
    output.write_text(svg, encoding="utf-8")
    return output


def rdkit_2d_svg(smiles: str, size: int = 420, legend: str | None = None) -> str:
    """Plain RDKit drawing from a SMILES string."""
    if size <= 0:
        raise ValueError("size must be a positive integer")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(size, size)
    options = drawer.drawOptions()
    options.clearBackground = True
    options.legendFontSize = 18

    drawer.DrawMolecule(mol, legend=legend or "")
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def save_rdkit_2d(
    smiles: str,
    output_path: str | Path,
    size: int = 420,
    legend: str | None = None,
) -> Path:
    """Save a plain RDKit SMILES depiction as SVG."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    svg = rdkit_2d_svg(smiles=smiles, size=size, legend=legend)
    output.write_text(svg, encoding="utf-8")
    return output


def draw_from_complex_data(
    complex_data: dict[str, object],
    size: int = 700,
    *,
    title: str | None = None,
) -> str:
    """Draw from a dict-like parsed complex representation."""
    parsed = ParsedComplex(
        metal=complex_data.get("metal"),
        ligands=complex_data.get("ligands", {}),
        complex_charge=complex_data.get("complex_charge", 0),
        counter_ions=complex_data.get("counter_ions", {}),
        raw_formula=complex_data.get("raw_formula", "manual_input"),
    )

    parsed.oxidation_state = complex_data.get("oxidation_state")
    parsed.coordination_number = complex_data.get(
        "coordination_number",
        sum(complex_data.get("ligands", {}).values()),
    )
    parsed.warnings = complex_data.get("warnings", [])
    parsed.errors = complex_data.get("errors", [])
    parsed.donor_atoms = complex_data.get("donor_atoms", {})
    parsed.ligand_denticity = complex_data.get("ligand_denticity", {})

    return diagram_2d_svg(parsed, size=size, title=title)


__all__ = [
    "parse_complex_input",
    "build_coordination_mol",
    "diagram_2d_svg",
    "save_diagram_2d",
    "rdkit_2d_svg",
    "save_rdkit_2d",
    "draw_from_complex_data",
]
