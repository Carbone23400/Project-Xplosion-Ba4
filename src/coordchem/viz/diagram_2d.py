"""
coordchem/viz/diagram_2d.py
---------------------------
Advanced 2D drawing for coordination complexes using RDKit.
"""

from __future__ import annotations

from html import escape
from importlib import import_module
from pathlib import Path
import re
from typing import Callable

from ..geometry import get_geometry
from ..parser import FormulaParseError, ParsedComplex, parse_formula
from .layout_2d import (
    Site,
    chelate_octahedral_sites,
    coordination_sites,
    edta_octahedral_sites,
    should_use_chelate_layout,
    should_use_edta_layout,
    should_use_tridentate_layout,
    tridentate_octahedral_sites,
)
from .ligand_data import (
    ABBREVIATED_MONODENTATE_LIGANDS,
    EXPLICIT_H_LIGANDS,
    INTERRUPTED_LIGAND_BONDS,
    LIGAND_DONOR_INDEX_OVERRIDES,
    LIGAND_SMILES,
    MONODENTATE_DISPLAY_LABELS,
    POLYDENTATE_DONOR_DISPLAY_LABELS,
)
from .transform_2d import (
    transform_acac,
    transform_edta,
    transform_monodentate,
    transform_polydentate,
)

try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "coordchem.viz.diagram_2d requires RDKit. Install rdkit first."
    ) from exc


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
    new_atom.SetIsAromatic(atom.GetIsAromatic())
    new_atom.SetFormalCharge(atom.GetFormalCharge())

    if atom_label is not None:
        new_atom.SetNumExplicitHs(0)
        new_atom.SetNoImplicit(True)
        new_atom.SetProp("atomLabel", atom_label)
        return new_atom

    new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
    new_atom.SetNoImplicit(atom.GetNoImplicit())
    return new_atom


def _monodentate_label(ligand_symbol: str, anchor: Site) -> str:
    """Return a formula label with donor side facing the metal."""
    donor_first, donor_last = MONODENTATE_DISPLAY_LABELS.get(
        ligand_symbol,
        (ligand_symbol, ligand_symbol),
    )
    return donor_last if anchor.x < 0 else donor_first


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


def _polydentate_donor_label(
    ligand_symbol: str,
    atom_idx: int,
    anchor: Site,
) -> str | None:
    """Return donor label for polydentate ligands."""
    ligand_labels = POLYDENTATE_DONOR_DISPLAY_LABELS.get(ligand_symbol, {})
    labels = ligand_labels.get(atom_idx)
    if labels is None:
        return None

    donor_first, donor_last = labels
    return donor_last if anchor.x < 0 else donor_first


def _get_2d_coords(mol: Chem.Mol) -> dict[int, tuple[float, float]]:
    """Return atom index to 2D coordinate mapping."""
    conf = mol.GetConformer()
    return {
        i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
        for i in range(mol.GetNumAtoms())
    }


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


def _is_interrupted_ligand_bond(
    ligand_symbol: str,
    begin_idx: int,
    end_idx: int,
) -> bool:
    """Return True when a ligand bond should be drawn as interrupted."""
    interrupted = INTERRUPTED_LIGAND_BONDS.get(ligand_symbol, set())
    pair = tuple(sorted((begin_idx, end_idx)))
    return pair in interrupted


def _interrupted_gap_fraction(
    coords: dict[int, tuple[float, float]],
    begin_idx: int,
    end_idx: int,
) -> float:
    """Return where an interrupted bond crosses the front vertical bond."""
    x1, _ = coords[begin_idx]
    x2, _ = coords[end_idx]

    if (x1 <= 0.0 <= x2) or (x2 <= 0.0 <= x1):
        denom = x2 - x1
        if abs(denom) > 1e-8:
            fraction = (0.0 - x1) / denom
            return max(0.15, min(0.85, fraction))

    return 0.5


def build_coordination_mol(
    complex_input: str | ParsedComplex,
    geometry_override: str | None = None,
) -> Chem.Mol:
    """Build a composite RDKit molecule for the complex."""
    parsed = parse_complex_input(complex_input)
    geometry = geometry_override or get_geometry(parsed)
    ligand_items = _expand_ligands(parsed)
    n_sites = sum(parsed.ligand_denticity.get(lig, 1) for lig in ligand_items)

    if should_use_edta_layout(parsed, ligand_items, geometry):
        sites = edta_octahedral_sites()
    elif should_use_tridentate_layout(parsed, ligand_items, geometry):
        sites = tridentate_octahedral_sites(len(ligand_items))
    elif should_use_chelate_layout(parsed, ligand_items, geometry):
        sites = chelate_octahedral_sites(len(ligand_items))
    else:
        sites = coordination_sites(geometry, n_sites)

    rw = Chem.RWMol()
    metal_atom = Chem.Atom(parsed.metal)
    metal_atom.SetNoImplicit(True)
    metal_idx = rw.AddAtom(metal_atom)

    global_coords: dict[int, tuple[float, float]] = {metal_idx: (0.0, 0.0)}
    site_cursor = 0

    for ligand_symbol in ligand_items:
        denticity_from_parser = parsed.ligand_denticity.get(ligand_symbol, 1)

        # Compact labels for monodentate ligands
        if (
            denticity_from_parser == 1
            and ligand_symbol in ABBREVIATED_MONODENTATE_LIGANDS
        ):
            anchors = tuple(sites[site_cursor: site_cursor + 1])
            if len(anchors) != 1:
                raise ValueError(
                    f"Not enough coordination sites for '{ligand_symbol}'."
                )
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

        # Full ligand drawing
        lig_mol = _make_ligand_mol(ligand_symbol)
        lig_coords = _get_2d_coords(lig_mol)
        donor_indices = _match_donor_indices(parsed, ligand_symbol, lig_mol)
        denticity = len(donor_indices)

        anchors = tuple(sites[site_cursor: site_cursor + denticity])
        if len(anchors) != denticity:
            raise ValueError(f"Not enough coordination sites for '{ligand_symbol}'.")
        site_cursor += denticity

        if ligand_symbol == "acac" and denticity == 2:
            transformed = transform_acac(
                lig_coords,
                donor_indices,
                anchors,
            )
        elif ligand_symbol in {"EDTA", "edta"} and denticity == 6:
            transformed = transform_edta(anchors)
        elif denticity == 1:
            transformed = transform_monodentate(
                lig_coords,
                donor_indices[0],
                anchors[0],
            )
        else:
            transformed = transform_polydentate(
                lig_coords,
                donor_indices,
                anchors,
                shrink_large_ligand=ligand_symbol not in {"acac", "bipy", "bpy"},
            )

        donor_anchors = dict(zip(donor_indices, anchors))

        idx_map: dict[int, int] = {}
        for atom in lig_mol.GetAtoms():
            local_idx = atom.GetIdx()
            atom_label = None

            if local_idx in donor_indices:
                atom_label = _polydentate_donor_label(
                    ligand_symbol=ligand_symbol,
                    atom_idx=local_idx,
                    anchor=donor_anchors[local_idx],
                )

            idx_map[local_idx] = rw.AddAtom(_copy_atom(atom, atom_label=atom_label))

        for bond in lig_mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            begin_global = idx_map[begin_idx]
            end_global = idx_map[end_idx]

            rw.AddBond(begin_global, end_global, bond.GetBondType())

            new_bond = rw.GetBondBetweenAtoms(begin_global, end_global)
            if _is_interrupted_ligand_bond(ligand_symbol, begin_idx, end_idx):
                new_bond.SetBoolProp("_coordchem_interrupted", True)
                new_bond.SetDoubleProp(
                    "_coordchem_gap_fraction",
                    _interrupted_gap_fraction(transformed, begin_idx, end_idx),
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


def _interrupted_bond_gap_fractions(mol: Chem.Mol) -> dict[int, float]:
    """Return interrupted bond indices and their SVG gap positions."""
    gaps: dict[int, float] = {}
    for bond in mol.GetBonds():
        if not bond.HasProp("_coordchem_interrupted"):
            continue

        gap = 0.5
        if bond.HasProp("_coordchem_gap_fraction"):
            gap = bond.GetDoubleProp("_coordchem_gap_fraction")

        gaps[bond.GetIdx()] = gap

    return gaps


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


def _get_name_module():
    """Load name2.py as the source of naming conventions."""
    try:
        return import_module("..name2", package=__package__)
    except Exception:
        return None


def _roman_numeral(value: int | None) -> str:
    """Return a Roman numeral using name2.py conventions."""
    name_module = _get_name_module()
    roman_number = getattr(name_module, "ROMAN_NUMBER", {}) if name_module else {}
    numerals = {number: roman for roman, number in roman_number.items()}
    numerals[0] = "0"
    return numerals.get(value, str(value))


def _ligand_count_prefix(count: int) -> str:
    """Return a coordination prefix using name2.py prefix values."""
    preferred_prefixes = {
        1: "",
        2: "di",
        3: "tri",
        4: "tetra",
        5: "penta",
        6: "hexa",
        7: "hepta",
        8: "octa",
    }

    prefix = preferred_prefixes.get(count)
    if prefix is None:
        return f"{count}-"

    if count == 1:
        return prefix

    name_module = _get_name_module()
    name2_prefixes = getattr(name_module, "PREFIXE", {}) if name_module else {}
    if name2_prefixes.get(prefix) == count:
        return prefix

    return f"{count}-"


def _metal_name(parsed: ParsedComplex) -> str:
    """Return the metal name, using -ate for anionic complex ions."""
    name_module = _get_name_module()
    if name_module is None:
        return parsed.metal

    METALS_NAME = getattr(name_module, "METALS_NAME", {})
    cation_name, anion_name = METALS_NAME.get(parsed.metal, (parsed.metal, parsed.metal))
    return anion_name if parsed.complex_charge < 0 else cation_name


def _coordination_compound_name(parsed: ParsedComplex) -> str:
    """Build a compact name from parser-enriched data and name2.py tables."""
    ligand_parts: list[str] = []
    for ligand, count in parsed.ligands.items():
        ligand_name = parsed.ligand_names.get(ligand, ligand)
        ligand_parts.append(f"{_ligand_count_prefix(count)}{ligand_name}")

    metal = _metal_name(parsed)
    oxidation = _roman_numeral(parsed.oxidation_state)
    return f"{''.join(ligand_parts)}{metal}({oxidation})"


def _add_svg_labels(svg: str, size: int, title: str | None) -> str:
    """Add white background and optional centered title."""
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


def _split_svg_bond_path(
    path_element: str,
    gap_center_fraction: float = 0.5,
    gap_fraction: float = 0.18,
) -> str:
    """Split a simple SVG bond path into two pieces with a local gap."""
    match = re.search(
        r"d='M ([\d.-]+),([\d.-]+) L ([\d.-]+),([\d.-]+)'",
        path_element,
    )
    if match is None:
        return path_element

    x1, y1, x2, y2 = (float(value) for value in match.groups())
    gap_center_fraction = max(0.0, min(1.0, gap_center_fraction))
    mid_x = x1 + (x2 - x1) * gap_center_fraction
    mid_y = y1 + (y2 - y1) * gap_center_fraction
    dx = (x2 - x1) * gap_fraction / 2
    dy = (y2 - y1) * gap_fraction / 2

    left_end = (mid_x - dx, mid_y - dy)
    right_start = (mid_x + dx, mid_y + dy)

    first_d = f"d='M {x1:.1f},{y1:.1f} L {left_end[0]:.1f},{left_end[1]:.1f}'"
    second_d = f"d='M {right_start[0]:.1f},{right_start[1]:.1f} L {x2:.1f},{y2:.1f}'"

    first = re.sub(r"d='M [^']+'", first_d, path_element, count=1)
    second = re.sub(r"d='M [^']+'", second_d, path_element, count=1)
    return first + "\n" + second


def _add_interrupted_bond_styles(svg: str, mols: list[Chem.Mol]) -> str:
    """Render selected bonds with a small gap where they pass behind."""
    interrupted: dict[int, float] = {}
    for mol in mols:
        interrupted.update(_interrupted_bond_gap_fractions(mol))

    for bond_idx, gap_center_fraction in interrupted.items():
        pattern = re.compile(
            rf"<path class='bond-{bond_idx} atom-\d+ atom-\d+'[^>]*/>"
        )
        svg = pattern.sub(
            lambda match: _split_svg_bond_path(
                match.group(0),
                gap_center_fraction=gap_center_fraction,
            ),
            svg,
        )

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
    display_title = title if title is not None else _coordination_compound_name(parsed)
    geometry_options = _geometry_options(geometry)
    ligand_items = _expand_ligands(parsed)
    is_edta_drawing = should_use_edta_layout(parsed, ligand_items, geometry)
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
    options.fixedBondLength = 42 if is_edta_drawing else 32
    options.addStereoAnnotation = False
    options.prepareMolsBeforeDrawing = False
    options.legendFontSize = 18

    if len(mols) == 1:
        drawer.DrawMolecule(mols[0], legend=geometry, confId=0)
    else:
        drawer.DrawMolecules(mols, legends=geometry_options, confIds=[0] * len(mols))

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = _add_interrupted_bond_styles(svg, mols)
    return _add_svg_labels(svg, size=size, title=display_title)


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
