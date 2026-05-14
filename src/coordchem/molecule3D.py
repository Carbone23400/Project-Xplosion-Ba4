"""
Backward-compatible 3D helpers.

The maintained implementation lives in ``coordchem.viz.structure_3d``.  This
module keeps older imports such as ``coordchem.molecule3D.ligand_3D`` working.
"""

from __future__ import annotations
from rdkit import Chem
from coordchem.viz.diagram_2d import parse_complex_input

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
from typing import Iterable, Optional, Sequence, Tuple
import math
from rdkit import Chem
from rdkit.Chem import AllChem

from coordchem.parser import ParsedComplex, parse_formula, KNOWN_LIGANDS
from coordchem.viz.ligand_data import LIGAND_DONOR_INDEX_OVERRIDES, LIGAND_SMILES
from coordchem.name import parse_name

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
def octahedral_positions_bidentate(sites):
    return [(sites[0],sites[2]),(sites[1],sites[4]),(sites[3],sites[5])]


def octahedral_positions_tridentate(sites):
    #return [(sites[0],sites[2],sites[1]),(sites[3],sites[5],sites[4])]
    return [(sites[0],sites[2],sites[1]),(sites[3],sites[4],sites[5])]

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


def find_donor_atom2(
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

#on peut prendre directement les infos données dans KNOWNLIGANDS a partir de .donor_atoms


def find_donor_atom(mol, donor_symbol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == donor_symbol:
            return atom.GetIdx()

    raise ValueError(f"No donor atom {donor_symbol} found in ligand")

        
 #function for the disosition of bidentate ligands 


def vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

def vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])

def vec_scale(v, s):
    return (v[0] * s, v[1] * s, v[2] * s)

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return (
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0],
    )

def norm(v):
    return (v[0]**2 + v[1]**2 + v[2]**2) ** 0.5

def unit(v):
    n = norm(v)
    if n == 0:
        return (1.0, 0.0, 0.0)
    return (v[0]/n, v[1]/n, v[2]/n)

def rotate_vector(v, axis, angle):
    import math

    axis = unit(axis)
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)

    term1 = vec_scale(v, cos_a)
    term2 = vec_scale(cross(axis, v), sin_a)
    term3 = vec_scale(axis, dot(axis, v) * (1 - cos_a))

    return vec_add(vec_add(term1, term2), term3)

def translate_bidentate_ligand(ligand_mol, donor_indices, target_sites):
    ligand_conf = ligand_mol.GetConformer()

    donor1 = ligand_conf.GetAtomPosition(donor_indices[0])
    donor2 = ligand_conf.GetAtomPosition(donor_indices[1])

    donor1 = (donor1.x, donor1.y, donor1.z)
    donor2 = (donor2.x, donor2.y, donor2.z)

    ligand_mid = (
    (donor1[0] + donor2[0]) / 2,
    (donor1[1] + donor2[1]) / 2,
    (donor1[2] + donor2[2]) / 2,
)

    target1, target2 = target_sites

    target_mid = (
    (target1[0] + target2[0]) / 2,
    (target1[1] + target2[1]) / 2,
    (target1[2] + target2[2]) / 2,
)
    push= 2.0
    outward=unit(target_mid)
    
    target_mid = (
        target_mid[0] + outward[0] * push,
        target_mid[1] + outward[1] * push,
        target_mid[2] + outward[2] * push,
    )

    ligand_axis = unit(vec_sub(donor2, donor1))
    target_axis = unit(vec_sub(target2, target1))

    rotation_axis=cross(ligand_axis, target_axis)
    rotation_axis_norm=norm(rotation_axis)

    if rotation_axis_norm == 0:
        angle = 0.0
        rotation_axis = (1.0, 0.0, 0.0)
    else:
        cos_angle = max(-1.0, min(1.0, dot(ligand_axis, target_axis)))
        angle = math.acos(cos_angle)

    coords = {}
    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        old_pos = ligand_conf.GetAtomPosition(idx)
        old_pos = (old_pos.x, old_pos.y, old_pos.z)

        centered = vec_sub(old_pos, ligand_mid)
        rotated = rotate_vector(centered, rotation_axis, angle)
        final_pos = vec_add(rotated, target_mid)

        coords[idx] = final_pos
    return coords
       
def translate_tridentate_ligand(ligand_mol, donor_indices, target_sites, reverse):
    ligand_conf = ligand_mol.GetConformer()

    donor_positions=[]
    for idx in donor_indices:
        pos=ligand_conf.GetAtomPosition(idx)
        donor_positions.append((pos.x,pos.y,pos.z))

    ligand_mid=(sum(pos[0] for pos in donor_positions)/3,
                sum(pos[1] for pos in donor_positions)/3,
                sum(pos[2] for pos in donor_positions)/3
                    )
    target_mid=(sum(pos[0] for pos in target_sites)/3,
                sum(pos[1] for pos in target_sites)/3,
                sum(pos[2] for pos in target_sites)/3
                    )
    outward=unit(target_mid)
    push=0.8

    target_mid=(target_mid[0] + outward[0] * push,
        target_mid[1] + outward[1] * push,
        target_mid[2] + outward[2] * push,)
    #shift= (target_mid[0] - ligand_mid[0],
     #   target_mid[1] - ligand_mid[1],
      #  target_mid[2] - ligand_mid[2],)
    coords={}
    for atom in ligand_mol.GetAtoms():
        idx=atom.GetIdx()
        old_pos=ligand_conf.GetAtomPosition(idx)
        old_pos=(old_pos.x, old_pos.y, old_pos.z)
        centered=vec_sub(old_pos, ligand_mid)
        if reverse:
            centered=(-centered[0], -centered[1], -centered[2])

        final_pos=vec_add(centered, target_mid)
        coords[idx]=final_pos

    return coords


# ---------------------------------------------------------------------------
# Complex-level builder
# ---------------------------------------------------------------------------
#à garder ou non car ne change pas grand chose
def thiocyanato_ligand_positions(ligand_mol, donor_idx,target,ligand_symbol):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    n_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "N"]
    s_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "S"] 
    c_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "C"] 
    
    if ligand_symbol=="NCS":
        n_idx=donor_idx
        coords[n_idx]=target
        outward=vec_scale(direction,0.8)
        side=vec_scale(u,0.7)

        c_idx=c_indices[0]
        c_pos=vec_add(target,vec_scale(direction,1.2))
        coords[c_idx]=c_pos
            
        s_idx=s_indices[0]
        s_pos=vec_add(c_pos,vec_scale(direction,1.6))
        coords[s_idx]=s_pos

    elif ligand_symbol=="SCN":
        s_idx=donor_idx
        coords[s_idx]=target
        outward=vec_scale(direction,0.8)
        side=vec_scale(u,0.7)

        c_idx=c_indices[0]
        c_pos=vec_add(target,vec_scale(direction,1.6))
        coords[c_idx]=c_pos
            
        n_idx=n_indices[0]
        n_pos=vec_add(target,vec_scale(direction,1.2))
        coords[n_idx]=n_pos

    for atom in ligand_mol.GetAtoms():
        idx= atom.GetIdx()
        if idx not in coords:
            coords[idx]=target
    return coords
def nitrito_ligand_positions(ligand_mol, donor_idx, target, ligand_symbol):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    n_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "N"]
    o_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "O"]
    if ligand_symbol=="NO2":
        n_idx=donor_idx
        coords[n_idx]=target
        outward=vec_scale(direction,0.8)
        side=vec_scale(u,0.7)

        if len(o_indices)>=2:
            coords[o_indices[0]]=vec_add(target,vec_add(outward,side))
            coords[o_indices[1]]=vec_add(target,vec_sub(outward,side))

    elif ligand_symbol=="ONO":
        o_idx=donor_idx
        coords[o_idx]=target
        outward=vec_scale(direction,0.8)
        side=vec_scale(u,0.7)

        coords[n_indices[0]]=vec_add(target,vec_add(outward,side))  
        other_o=[idx for idx in o_indices if idx!=donor_idx]
        if other_o:
            o_pos=vec_add(coords[n_indices[0]],vec_add(vec_scale(direction,0.8),vec_scale(u,0.7),),)
            coords[other_o[0]]=o_pos
    for atom in ligand_mol.GetAtoms():
        idx= atom.GetIdx()
        if idx not in coords:
            coords[idx]=target
    return coords


def ammonia_ligand_positions(ligand_mol, donor_idx, target, nh_distance=1.0, spread=0.8):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    v = unit(cross(direction, u))

    h_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "H"
    ]

    import math

    for i, h_idx in enumerate(h_indices[:3]):
        angle = 2 * math.pi * i / 3

        radial = vec_add(
            vec_scale(u, math.cos(angle) * spread),
            vec_scale(v, math.sin(angle) * spread),
        )

        outward = vec_scale(direction, 0.4)

        h_pos = vec_add(center, vec_add(radial, outward))
        coords[h_idx] = h_pos

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in coords:
            coords[idx] = center

    return coords



def pyridine_ligand_positions(ligand_mol, donor_idx_local,target,ligand_number=0):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx_local] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    v = unit(cross(direction, u))
    ring_radius=0.75
    outward_distance=1.3
    angle_offset=ligand_number*0.7
    all_heavy_atoms = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() != "H"
    ]
    heavy_atoms = [donor_idx_local] + [
    idx for idx in all_heavy_atoms
    if idx != donor_idx_local
]
    ring_center=vec_add(target,vec_scale(direction, ring_radius))
    coords[donor_idx_local]=target
    angles=[math.pi, 2*math.pi/3, math.pi/3,0.0,-math.pi/3,-2*math.pi/3]

    for idx, angle in zip(heavy_atoms, angles):
        ring_pos = vec_add(
                vec_scale(direction, math.cos(angle) * ring_radius),
                vec_scale(u, math.sin(angle) * ring_radius))
        coords[idx] = vec_add(ring_center, ring_pos)
    for atom in ligand_mol.GetAtoms():
        idx=atom.GetIdx()
        if idx in coords:
            continue 
        if atom.GetSymbol() == "H":
            neighbors = atom.GetNeighbors()
            if neighbors:
                heavy_idx = neighbors[0].GetIdx()
                base = coords.get(heavy_idx, target)
                outward=unit(vec_sub(base, ring_center))
                #base = coords.get(heavy_idx, target)
                #h_pos = vec_add(base, vec_scale(direction, 0.45))
                #coords[idx] = h_pos
                coords[idx]=vec_add(base, vec_scale(outward,0.45))
            else:
                #coords[idx] = vec_add(target, vec_scale(direction, outward_distance + 0.5))
                coords[idx]=target
    return coords



def trimethyl_ligand_positions(ligand_mol, donor_idx, target, nh_distance=1.0, spread=0.8):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    v = unit(cross(direction, u))

    methyl_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "C"
    ]



    for i, c_idx in enumerate(methyl_indices[:3]):
        angle = 2 * math.pi * i / 3

        radial = vec_add(
            vec_scale(u, math.cos(angle) * spread),
            vec_scale(v, math.sin(angle) * spread),
        )

        outward = vec_scale(direction, 0.8)

        ch3_pos = vec_add(center, vec_add(radial, outward))
        coords[c_idx] = ch3_pos
        carbon_atom=ligand_mol.GetAtomWithIdx(c_idx)
        h_neighbors= [neighbor.GetIdx() for neighbor in carbon_atom.GetNeighbors() if neighbor.GetSymbol() == "H"]
        for j, h_idx in enumerate(h_neighbors):
            h_angle = angle + 2 * math.pi * j / 3

            h_radial = vec_add(
        vec_scale(u, math.cos(h_angle) * 0.45),
        vec_scale(v, math.sin(h_angle) * 0.45),
    )

            h_outward = vec_scale(direction, 0.35)

            coords[h_idx] = vec_add(
        ch3_pos,
        vec_add(h_radial, h_outward),
    )


    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in coords:
            coords[idx] = center

    return coords

def triethyl_ligand_positions(ligand_mol, donor_idx, target, nh_distance=1.0, spread=0.8):
    direction = unit(target)

    center = target
    coords = {}
    coords[donor_idx] = center

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    v = unit(cross(direction, u))

    import math
    p_atom=ligand_mol.GetAtomWithIdx(donor_idx)

    first_carbons=[neighbor.GetIdx() for neighbor in p_atom.GetNeighbors() if neighbor.GetSymbol()=="C"]
    
    
    for i, c_idx in enumerate(first_carbons[:3]):
        second_carbons=[neighbor.GetIdx() for neighbor in ligand_mol.GetAtomWithIdx(c_idx).GetNeighbors() if neighbor.GetSymbol()=="C"]
        
        angle = 2 * math.pi * i / 3

        radial = vec_add(
            vec_scale(u, math.cos(angle) * spread),
            vec_scale(v, math.sin(angle) * spread),
        )

        outward = vec_scale(direction, 0.8)

        eth_pos = vec_add(center, vec_add(radial, outward))
        coords[c_idx] = eth_pos

        if second_carbons:
            c2_idx=second_carbons[0]
            eth2_pos=vec_add(center,vec_add(radial,vec_scale(direction,1.6)))
            coords[c2_idx]=eth2_pos
            carbon_atom2=ligand_mol.GetAtomWithIdx(c2_idx)
        
            h_neighbors2= [neighbor.GetIdx() for neighbor in carbon_atom2.GetNeighbors() if neighbor.GetSymbol() == "H"]
            for j, h2_idx in enumerate(h_neighbors2):
                h_angle2 = angle + math.pi / 3 + 2 * math.pi * j / 3

                h_radial2 = vec_add(
        vec_scale(u, math.cos(h_angle2) * 0.45),
        vec_scale(v, math.sin(h_angle2) * 0.45),
    )

                h_outward2 = vec_scale(direction, 0.45)

                coords[h2_idx] = vec_add(
        eth2_pos,
        vec_add(h_radial2, h_outward2),
    )

        carbon_atom=ligand_mol.GetAtomWithIdx(c_idx)
        
        h_neighbors= [neighbor.GetIdx() for neighbor in carbon_atom.GetNeighbors() if neighbor.GetSymbol() == "H"]
        for j, h_idx in enumerate(h_neighbors):
            h_angle = angle + 2 * math.pi * j / 3

            h_radial = vec_add(
        vec_scale(u, math.cos(h_angle) * 0.45),
        vec_scale(v, math.sin(h_angle) * 0.45),
    )

            h_outward = vec_scale(direction, 0.2)

            coords[h_idx] = vec_add(
        eth_pos,
        vec_add(h_radial, h_outward),
    )
        

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in coords:
            coords[idx] = center

    return coords

def terpyridine_ligand_positions(ligand_mol, donor_indices, tridentate_index):
    ligand_conf = ligand_mol.GetConformer()

    coords = {}

    donor_positions = []
    for idx in donor_indices:
        pos = ligand_conf.GetAtomPosition(idx)
        donor_positions.append((pos.x, pos.y, pos.z))

    ligand_mid = (
        sum(p[0] for p in donor_positions) / 3,
        sum(p[1] for p in donor_positions) / 3,
        sum(p[2] for p in donor_positions) / 3,
    )

    if tridentate_index == 0:
        center = (0.0, 0.0, 1.8)
        angle = 0.0
    else:
        center = (0.0, 0.0, -1.8)
        angle = 3.14159

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        pos = ligand_conf.GetAtomPosition(idx)
        old = (pos.x, pos.y, pos.z)

        centered = vec_sub(old, ligand_mid)

        rotated = (
            centered[0] * math.cos(angle) - centered[1] * math.sin(angle),
            centered[0] * math.sin(angle) + centered[1] * math.cos(angle),
            centered[2],
        )

        coords[idx] = vec_add(rotated, center)

    return coords

def terpyridine_ligand_positions2(ligand_mol, donor_indices, tridentate_index):
    coords = {}

    if tridentate_index == 0:
        z = 1.8
        y_shift = 0.0
        sign = 1.0
    else:
        z = -1.8
        y_shift = 0.0
        sign = -1.0

    donor_targets = [
        (-1.8, y_shift, z),
        (0.0, y_shift, z),
        (1.8, y_shift, z),
    ]

    if tridentate_index == 1:
        donor_targets = list(reversed(donor_targets))

    for donor_idx, target in zip(donor_indices, donor_targets):
        coords[donor_idx] = target

    ring_centers = [
    (-2.2, 0.9 * sign, z),
    (0.0, 1.1 * sign, z),
    (2.2, 0.9 * sign, z),
]

    ring_radius = 0.65
    heavy_count = 0

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()

        if idx in coords:
            continue

        symbol = atom.GetSymbol()

        ring_id = heavy_count // 5
        if ring_id > 2:
            ring_id = 2

        pos_in_ring = heavy_count % 5
        angle = 2 * math.pi * pos_in_ring / 5

        cx, cy, cz = ring_centers[ring_id]

        x = cx + ring_radius * math.cos(angle)
        y = cy + sign * ring_radius * math.sin(angle)

        if symbol == "H":
            x = cx + (ring_radius + 0.45) * math.cos(angle)
            y = cy + sign * (ring_radius + 0.45) * math.sin(angle)
        else:
            heavy_count += 1

        coords[idx] = (x, y, z)



    return coords


def build_complex_3d(
    parsed: ParsedComplex,
    distance: float = 2.0,
    geometry: str | None = None,
) -> Chem.Mol:

    if geometry is None:
        try:
            from coordchem.geometry import predict_geometry
            geometry = predict_geometry(parsed)
        except Exception:
            geometry = "octahedral"

    cn = parsed.coordination_number or 6
    sites = geometry_positions(geometry, cn, distance=distance)

    rw = Chem.RWMol()
    coords: dict[int, Position] = {}

    metal_atom = Chem.Atom(parsed.metal)
    metal_atom.SetNoImplicit(True)
    metal_idx = rw.AddAtom(metal_atom)
    coords[metal_idx] = (0.0, 0.0, 0.0)

    site_index = 0

    for ligand_symbol, count in parsed.ligands.items():
        smiles = LIGAND_SMILES.get(ligand_symbol)

        if smiles is None:
            for _ in range(count):
                if site_index >= len(sites):
                    break

                placeholder = Chem.Atom("X")
                idx = rw.AddAtom(placeholder)
                coords[idx] = sites[site_index]
                rw.AddBond(metal_idx, idx, Chem.BondType.DATIVE)
                site_index += 1

            continue

        donor_symbol = parsed.donor_atoms.get(ligand_symbol, "?")
        donor_overrides = LIGAND_DONOR_INDEX_OVERRIDES.get(ligand_symbol, ())
        denticity = parsed.ligand_denticity.get(ligand_symbol, 1)

        for _ in range(count):
            if site_index >= len(sites):
                break

            try:
                ligand_mol = build_ligand_3d(smiles)
            except ValueError:
                continue

            if denticity == 2 and len(donor_overrides) >= 2:
                bidentate_index = site_index // 2

                if "octahedral" in geometry.lower() and len(sites) >= 6:
                    site_pairs = octahedral_positions_bidentate(sites)
                else:
                    site_pairs = [
                        (sites[i], sites[i + 1])
                        for i in range(0, len(sites) - 1, 2)
                    ]

                if bidentate_index >= len(site_pairs):
                    break

                target_sites = (sites[site_index],sites[site_index+1])
                donor_indices = donor_overrides[:2]

                local_coords = translate_bidentate_ligand(
                    ligand_mol,
                    donor_indices,
                    target_sites,
                )

                offset = rw.GetNumAtoms()

                for atom in ligand_mol.GetAtoms():
                    global_idx = rw.AddAtom(atom)
                    coords[global_idx] = local_coords[atom.GetIdx()]

                for bond in ligand_mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx() + offset
                    a2 = bond.GetEndAtomIdx() + offset
                    rw.AddBond(a1, a2, bond.GetBondType())

                for donor_idx in donor_indices:
                    rw.AddBond(
                        metal_idx,
                        donor_idx + offset,
                        Chem.BondType.DATIVE,
                    )

                site_index += 2
                continue

            elif denticity==3 and len(donor_overrides)>=3:
                tridentate_index=site_index//3
                target_sites=octahedral_positions_tridentate(sites)[tridentate_index]
                donor_indices = donor_overrides[:3]
                if ligand_symbol in ("tpy", "terpy"):
                    local_coords = terpyridine_ligand_positions(
        ligand_mol,
        donor_overrides[:3],
        tridentate_index,
    )
                else:
                    local_coords = translate_tridentate_ligand(
        ligand_mol,
        donor_overrides[:3],
        target_sites,
        reverse=(tridentate_index % 2 == 1),
    )

                offset = rw.GetNumAtoms()

                for atom in ligand_mol.GetAtoms():
                    global_idx = rw.AddAtom(atom)
                    coords[global_idx] = local_coords[atom.GetIdx()]

                for bond in ligand_mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx() + offset
                    a2 = bond.GetEndAtomIdx() + offset
                    rw.AddBond(a1, a2, bond.GetBondType())

                for donor_idx in donor_indices:
                    rw.AddBond(metal_idx, donor_idx + offset, Chem.BondType.DATIVE,)

                site_index += 3
                continue
            
            primary_override = donor_overrides[0] if donor_overrides else None

            if primary_override is not None and primary_override < ligand_mol.GetNumAtoms():
                donor_idx_local = primary_override
            else:
                try:
                    donor_idx_local = find_donor_atom(ligand_mol, donor_symbol)
                except ValueError:
                    continue
#maintenant c'est le cas monodentate
            target = sites[site_index]

            if ligand_symbol=="NH3":
                local_coords=ammonia_ligand_positions(ligand_mol, donor_idx_local,target)
            elif ligand_symbol=="py":
                local_coords=pyridine_ligand_positions(ligand_mol, donor_idx_local,target, site_index)    
            elif ligand_symbol=="PMe3":
                local_coords=trimethyl_ligand_positions(ligand_mol, donor_idx_local,target)
            elif ligand_symbol=="PEt3":
                local_coords=triethyl_ligand_positions(ligand_mol,donor_idx_local,target)
            #elif ligand_symbol=="PPh3":
                #local_coords=triphenylphosphine_ligand_position(ligand_mol, donor_idx_local,target, site_index)
            elif ligand_symbol in ("NO2","ONO"):
                local_coords=nitrito_ligand_positions(ligand_mol, donor_idx_local,target,ligand_symbol)
            elif ligand_symbol in ("NCS","SCN"):
                local_coords=thiocyanato_ligand_positions(ligand_mol, donor_idx_local,target,ligand_symbol)    
            elif ligand_mol.GetNumAtoms() == 2:
                local_coords = diatomic_ligand_positions(
                    ligand_mol,
                    donor_idx_local,
                    target,
                    bond_length=1.2,
                )
            else:
                ligand_conf = ligand_mol.GetConformer()
                donor_pos = ligand_conf.GetAtomPosition(donor_idx_local)

                local_coords = {}

                for atom in ligand_mol.GetAtoms():
                    old_pos = ligand_conf.GetAtomPosition(atom.GetIdx())
                    local_coords[atom.GetIdx()] = (
                        old_pos.x - donor_pos.x + target[0],
                        old_pos.y - donor_pos.y + target[1],
                        old_pos.z - donor_pos.z + target[2],
                    )

            offset = rw.GetNumAtoms()

            for atom in ligand_mol.GetAtoms():
                global_idx = rw.AddAtom(atom)
                coords[global_idx] = local_coords[atom.GetIdx()]

            for bond in ligand_mol.GetBonds():
                a1 = bond.GetBeginAtomIdx() + offset
                a2 = bond.GetEndAtomIdx() + offset
                rw.AddBond(a1, a2, bond.GetBondType())

            rw.AddBond(
                metal_idx,
                donor_idx_local + offset,
                Chem.BondType.DATIVE,
            )

            site_index += 1

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())

    for atom_idx, position in coords.items():
        conf.SetAtomPosition(atom_idx, position)
    mol.AddConformer(conf, assignId=True)
    return mol

   



#fonctions qui display les ligands diatiomiques pour que le deuxieme atome soit bien positionné
def normalize(v):
    length=(v[0]**2+v[1]**2+v[2]**2)**0.5
    if length==0:
        return(1.0,0.0,0.0)
    return (v[0]/length,v[1]/length,v[2]/length)


def diatomic_ligand_positions(ligand_mol, donor_idx, target, bond_length=1.2):
    direction = normalize(target)
    coords = {}
    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx == donor_idx:
            coords[idx] = target
        else:
            coords[idx] = (
                target[0] + direction[0] * bond_length,
                target[1] + direction[1] * bond_length,
                target[2] + direction[2] * bond_length,
            )
    return coords


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def _to_parsed(complex_or_formula) -> ParsedComplex:
    """Best-effort coercion to a ``ParsedComplex``."""
    if isinstance(complex_or_formula, ParsedComplex):
        return complex_or_formula
    if isinstance(complex_or_formula, str):
        return parse_complex_input(complex_or_formula)
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

    block = Chem.MolToMolBlock(mol, kekulize=False)#pour ne pas forcer la conversion des cycles aromatiques
    view = py3Dmol.view(width=width, height=height)
    view.removeAllModels()
    view.addModel(block, "sdf")
    view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
    view.setStyle({}, {"stick": {"radius": 0.12}, "sphere": {"scale": 0.18}})

    view.setBackgroundColor("white")
    view.zoomTo()
    return view
#fonction utile pour display la molecule sur un notebook

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

