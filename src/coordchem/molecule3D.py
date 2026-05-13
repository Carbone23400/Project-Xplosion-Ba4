"""
Backward-compatible 3D helpers.

The maintained implementation lives in ``coordchem.viz.structure_3d``.  This
module keeps older imports such as ``coordchem.molecule3D.ligand_3D`` working.
"""

from __future__ import annotations

from rdkit import Chem

from .viz.structure_3d import (
    build_complex_3d,
    build_ligand_3d,
    complex_3d_html,
    find_donor_atom,
    geometry_positions,
    octahedral_positions,
    tetrahedral_positions,
    view_complex_3d,
)


def ligand_3D(ligand_smile: str) -> Chem.Mol:
    """Compatibility alias for ``build_ligand_3d``."""
    return build_ligand_3d(ligand_smile)


def drawit(mol: Chem.Mol, p=None, confId: int = -1):
    """Display an RDKit molecule with py3Dmol, matching the old helper name."""
    import py3Dmol

    mol_block = Chem.MolToMolBlock(mol, confId=confId)
    if p is None:
        p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()
    p.addModel(mol_block, "sdf")
    p.setStyle({"stick": {}})
    p.setBackgroundColor("white")
    p.zoomTo()
    return p.show()


__all__ = [
    "build_complex_3d",
    "build_ligand_3d",
    "complex_3d_html",
    "drawit",
    "find_donor_atom",
    "geometry_positions",
    "ligand_3D",
    "octahedral_positions",
    "tetrahedral_positions",
    "view_complex_3d",
]
