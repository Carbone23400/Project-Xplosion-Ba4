"""
coordchem/viz/structure_3d.py
-----------------------------
3D structure construction for coordination complexes.

The module turns a :class:`coordchem.parser.ParsedComplex` into an RDKit
:class:`Mol` object carrying a 3D conformer that can be displayed with
``py3Dmol`` (e.g. inside a Streamlit app).

Public API
~~~~~~~~~~
``build_ligand_3d(smiles)``
    Build an RDKit ``Mol`` for a single ligand with explicit hydrogens
    and a single 3D conformer.

``find_donor_atom(mol, donor_symbol, override=None)``
    Locate the donor atom index inside a ligand ``Mol``.

``geometry_positions(geometry, n, distance=2.0)``
    Return ``n`` coordination-site positions arranged according to
    ``geometry`` ("octahedral", "tetrahedral", "square planar", ...).

``build_complex_3d(parsed, distance=2.0)``
    Assemble the metal + ligands into a single RDKit ``Mol`` with a
    conformer. Suitable for display.

``view_complex_3d(parsed, width=400, height=400)``
    Return a ``py3Dmol.view`` ready to be displayed.

``complex_3d_html(parsed, width=400, height=400)``
    Return a self-contained HTML snippet — convenient for embedding in
    Streamlit through ``st.components.v1.html``.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from ..parser import ParsedComplex, parse_formula
from .ligand_data import LIGAND_DONOR_INDEX_OVERRIDES, LIGAND_SMILES


Position = Tuple[float, float, float]


# ---------------------------------------------------------------------------
# Geometry placeholders
# ---------------------------------------------------------------------------

def octahedral_positions(distance: float = 2.0) -> list[Position]:
    """Return six orthogonal positions on the +/- x, y, z axes."""
    return [
        ( distance,  0.0,       0.0),
        (-distance,  0.0,       0.0),
        ( 0.0,       distance,  0.0),
        ( 0.0,      -distance,  0.0),
        ( 0.0,       0.0,       distance),
        ( 0.0,       0.0,      -distance),
    ]


def tetrahedral_positions(distance: float = 2.0) -> list[Position]:
    """Return four positions pointing to the corners of a tetrahedron."""
    a = distance / (3 ** 0.5)
    return [
        ( a,  a,  a),
        ( a, -a, -a),
        (-a,  a, -a),
        (-a, -a,  a),
    ]


def square_planar_positions(distance: float = 2.0) -> list[Position]:
    """Return four positions in the xy plane."""
    return [
        ( distance, 0.0, 0.0),
        (-distance, 0.0, 0.0),
        ( 0.0,  distance, 0.0),
        ( 0.0, -distance, 0.0),
    ]


def linear_positions(distance: float = 2.0) -> list[Position]:
    """Return two collinear positions along x."""
    return [( distance, 0.0, 0.0), (-distance, 0.0, 0.0)]


def trigonal_planar_positions(distance: float = 2.0) -> list[Position]:
    """Return three positions in the xy plane, 120° apart."""
    import math
    return [
        (distance * math.cos(math.radians(angle)),
         distance * math.sin(math.radians(angle)),
         0.0)
        for angle in (0, 120, 240)
    ]


def trigonal_bipyramidal_positions(distance: float = 2.0) -> list[Position]:
    """Three equatorial + two axial positions."""
    return trigonal_planar_positions(distance) + [
        (0.0, 0.0,  distance),
        (0.0, 0.0, -distance),
    ]


# Mapping from geometry label (as produced by ``coordchem.geometry``) to
# the corresponding position generator. ``geometry_positions`` does a
# best-effort lookup so that ambiguous labels like
# "trigonal bipyramidal or square pyramidal" still produce something
# reasonable.
_GEOMETRY_BUILDERS = {
    "linear": linear_positions,
    "trigonal planar": trigonal_planar_positions,
    "tetrahedral": tetrahedral_positions,
    "square planar": square_planar_positions,
    "trigonal bipyramidal": trigonal_bipyramidal_positions,
    "octahedral": octahedral_positions,
}


def geometry_positions(
    geometry: str | None,
    n: int,
    distance: float = 2.0,
) -> list[Position]:
    """
    Return ``n`` coordination-site positions matching ``geometry``.

    Falls back to an octahedral arrangement (truncated/extended to ``n``)
    when the geometry label is unknown — useful as a TODO placeholder for
    geometries we have not implemented yet.
    """
    if geometry:
        for label, builder in _GEOMETRY_BUILDERS.items():
            if label in geometry.lower():
                positions = builder(distance)
                if len(positions) >= n:
                    return positions[:n]
                # not enough sites: extend with octahedral ones
                extra = octahedral_positions(distance)[: n - len(positions)]
                return positions + extra

    # Default: octahedral, padded if needed
    base = octahedral_positions(distance)
    if n <= len(base):
        return base[:n]
    # TODO: implement explicit positions for CN > 6
    return base + [(0.0, 0.0, 0.0)] * (n - len(base))


# ---------------------------------------------------------------------------
# Ligand-level helpers
# ---------------------------------------------------------------------------

def build_ligand_3d(smiles: str) -> Chem.Mol:
    """
    Build an RDKit ``Mol`` from ``smiles`` with explicit Hs and a
    single 3D conformer.

    Raises
    ------
    ValueError
        If the SMILES cannot be parsed or a 3D conformer cannot be embedded.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles!r}")

    mol = Chem.AddHs(mol)

    # Ions like [Cl-] or [O-2] have a single atom — nothing to embed.
    if mol.GetNumAtoms() <= 1:
        conf = Chem.Conformer(mol.GetNumAtoms())
        if mol.GetNumAtoms() == 1:
            conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
        mol.AddConformer(conf, assignId=True)
        return mol

    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status != 0:
        # Retry with random coordinates as fallback
        status = AllChem.EmbedMolecule(
            mol, randomSeed=42, useRandomCoords=True
        )
        if status != 0:
            raise ValueError(f"3D embedding failed for SMILES: {smiles!r}")

    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        # Optimization is best-effort: keep the embedded geometry
        pass

    return mol


def find_donor_atom(
    mol: Chem.Mol,
    donor_symbol: str,
    override: Optional[int] = None,
) -> int:
    """
    Return the atom index of the donor inside a ligand ``Mol``.

    ``override`` is the explicit index from
    :data:`coordchem.viz.ligand_data.LIGAND_DONOR_INDEX_OVERRIDES`
    when available. Otherwise the first atom matching ``donor_symbol``
    is returned.

    Multi-donor symbols like ``"N/O"`` (used by EDTA) try each candidate
    in order.

    Raises
    ------
    ValueError
        When no atom with the requested element symbol exists.
    """
    if override is not None:
        if 0 <= override < mol.GetNumAtoms():
            return override

    candidates: Iterable[str]
    if donor_symbol == "?" or not donor_symbol:
        # Unknown: take the first heavy atom
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "H":
                return atom.GetIdx()
        raise ValueError("Ligand has no non-H atoms")

    candidates = [s.strip() for s in donor_symbol.split("/") if s.strip()]
    for symbol in candidates:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == symbol:
                return atom.GetIdx()

    raise ValueError(f"No donor atom matching {donor_symbol!r} found in ligand")


# ---------------------------------------------------------------------------
# Complex-level builder
# ---------------------------------------------------------------------------

def build_complex_3d(
    parsed: ParsedComplex,
    distance: float = 2.0,
    geometry: str | None = None,
) -> Chem.Mol:
    """
    Assemble a 3D RDKit ``Mol`` for the given :class:`ParsedComplex`.

    The metal sits at the origin and ligand donor atoms are placed on
    coordination sites generated by :func:`geometry_positions`.

    Parameters
    ----------
    parsed:
        Parsed coordination complex.
    distance:
        Metal–donor distance in Å.
    geometry:
        Optional override for the geometry label. If ``None``, the
        prediction returned by :func:`coordchem.geometry.predict_geometry`
        is used. Falls back to "octahedral" when unavailable.
    """
    if geometry is None:
        try:
            from ..geometry import predict_geometry
            geometry = predict_geometry(parsed)
        except Exception:
            geometry = "octahedral"

    cn = parsed.coordination_number or 6
    sites = geometry_positions(geometry, cn, distance=distance)

    rw = Chem.RWMol()
    coords: dict[int, Position] = {}

    # Metal at origin
    metal_atom = Chem.Atom(parsed.metal)
    metal_atom.SetNoImplicit(True)
    metal_idx = rw.AddAtom(metal_atom)
    coords[metal_idx] = (0.0, 0.0, 0.0)

    site_index = 0

    for ligand_symbol, count in parsed.ligands.items():
        smiles = LIGAND_SMILES.get(ligand_symbol)
        if smiles is None:
            # Unknown ligand: skip, but record a placeholder atom on the site
            # so the user still sees something where the ligand should be.
            for _ in range(count):
                if site_index >= len(sites):
                    break
                placeholder = Chem.Atom("X")  # dummy
                idx = rw.AddAtom(placeholder)
                coords[idx] = sites[site_index]
                rw.AddBond(metal_idx, idx, Chem.BondType.DATIVE)
                site_index += 1
            continue

        donor_symbol = parsed.donor_atoms.get(ligand_symbol, "?")
        donor_overrides: Sequence[int] = LIGAND_DONOR_INDEX_OVERRIDES.get(
            ligand_symbol, ()
        )

        for _ in range(count):
            if site_index >= len(sites):
                break

            try:
                ligand_mol = build_ligand_3d(smiles)
            except ValueError:
                continue  # skip a ligand that fails to embed

            # Pick the primary donor — the first override if any,
            # else the first atom matching ``donor_symbol``.
            primary_override = (
                donor_overrides[0] if donor_overrides else None
            )
            try:
                donor_idx_local = find_donor_atom(
                    ligand_mol, donor_symbol, override=primary_override
                )
            except ValueError:
                continue

            target = sites[site_index]
            ligand_conf = ligand_mol.GetConformer()
            donor_pos = ligand_conf.GetAtomPosition(donor_idx_local)

            offset = rw.GetNumAtoms()

            # Copy atoms with translated coordinates so that the donor
            # atom lands on ``target``.
            for atom in ligand_mol.GetAtoms():
                global_idx = rw.AddAtom(atom)
                old_pos = ligand_conf.GetAtomPosition(atom.GetIdx())
                coords[global_idx] = (
                    old_pos.x - donor_pos.x + target[0],
                    old_pos.y - donor_pos.y + target[1],
                    old_pos.z - donor_pos.z + target[2],
                )

            # Copy intra-ligand bonds
            for bond in ligand_mol.GetBonds():
                a1 = bond.GetBeginAtomIdx() + offset
                a2 = bond.GetEndAtomIdx() + offset
                rw.AddBond(a1, a2, bond.GetBondType())

            # Metal–donor dative bond for the primary site
            rw.AddBond(
                metal_idx,
                donor_idx_local + offset,
                Chem.BondType.DATIVE,
            )
            site_index += 1

            # TODO: properly distribute remaining donors of polydentate
            # ligands across multiple sites. For now we just consume one
            # site per ligand instance and rely on the geometry/2D
            # diagram to communicate the denticity.

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_idx, position in coords.items():
        conf.SetAtomPosition(atom_idx, position)
    mol.AddConformer(conf, assignId=True)

    return mol


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def _to_parsed(complex_or_formula) -> ParsedComplex:
    """Best-effort coercion to a ``ParsedComplex``."""
    if isinstance(complex_or_formula, ParsedComplex):
        return complex_or_formula
    if isinstance(complex_or_formula, str):
        return parse_formula(complex_or_formula)
    parsed = getattr(complex_or_formula, "parsed", None)
    if isinstance(parsed, ParsedComplex):
        return parsed
    raise TypeError(
        "Expected a ParsedComplex, Complex, or formula string"
    )


def view_complex_3d(
    complex_or_formula,
    width: int = 400,
    height: int = 400,
    distance: float = 2.0,
):
    """
    Return a ``py3Dmol.view`` for the complex.

    The caller is responsible for displaying the view (``view.show()``
    in a notebook or ``view._make_html()`` for embedding).
    """
    import py3Dmol  # local import: keep py3Dmol optional

    parsed = _to_parsed(complex_or_formula)
    mol = build_complex_3d(parsed, distance=distance)

    block = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=width, height=height)
    view.removeAllModels()
    view.addModel(block, "sdf")
    view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
    view.setBackgroundColor("white")
    view.zoomTo()
    return view


def complex_3d_html(
    complex_or_formula,
    width: int = 400,
    height: int = 400,
    distance: float = 2.0,
) -> str:
    """Return self-contained HTML for embedding the 3D view."""
    view = view_complex_3d(
        complex_or_formula,
        width=width,
        height=height,
        distance=distance,
    )
    return view._make_html()
