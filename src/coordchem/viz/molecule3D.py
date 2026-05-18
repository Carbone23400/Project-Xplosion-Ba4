
from typing import Iterable, Optional, Sequence, Tuple
import math
from rdkit import Chem
from rdkit.Chem import AllChem

from coordchem.parser import ParsedComplex, parse_formula, KNOWN_LIGANDS
from coordchem.viz.ligand_data import (
    LIGAND_DONOR_INDEX_OVERRIDES,
    LIGAND_SMILES,
    donor_index_overrides_for_ligand,
    is_short_bidentate_ligand,
)
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


def square_antiprismatic_positions(distance: float = 2.0) -> list[Position]:
    """Return eight sites arranged as two staggered squares."""
    z = distance * 0.45
    radius = (distance * distance - z * z) ** 0.5
    positions = []
    for z_value, offset in ((z, math.pi / 4), (-z, 0.0)):
        for i in range(4):
            angle = offset + i * math.pi / 2
            positions.append(
                (
                    radius * math.cos(angle),
                    radius * math.sin(angle),
                    z_value,
                )
            )
    return positions


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

def square_antiprismatic_positions(distance: float = 2.0) -> list[Position]:
    """Return eight positions for a square antiprismatic geometry."""
    import math

    z = distance * 0.55
    r = (distance**2 - z**2) ** 0.5

    bottom_square = []
    top_square = []

    for i in range(4):
        angle = math.pi / 2 * i

        bottom_square.append(
            (
                r * math.cos(angle),
                r * math.sin(angle),
                -z,
            )
        )

    for i in range(4):
        angle = math.pi / 4 + math.pi / 2 * i

        top_square.append(
            (
                r * math.cos(angle),
                r * math.sin(angle),
                z,
            )
        )

    return bottom_square + top_square


def dodecahedral_positions(distance: float = 2.0) -> list[Position]:
    """Return eight approximate positions for a dodecahedral CN=8 geometry."""
    import math

    z1 = distance * 0.75
    z2 = distance * 0.25

    r1 = (distance**2 - z1**2) ** 0.5
    r2 = (distance**2 - z2**2) ** 0.5

    return [
        ( r1,  0.0,  z1),
        ( 0.0,  r1,  z1),
        (-r1,  0.0,  z1),
        ( 0.0, -r1,  z1),

        ( r2,  r2, -z2),
        (-r2,  r2, -z2),
        (-r2, -r2, -z2),
        ( r2, -r2, -z2),
    ]
def pentagonal_bipyramidal_positions(distance: float = 2.0) -> list[Position]:
    """Return seven positions: five equatorial + two axial."""
    import math

    equatorial = [
        (
            distance * math.cos(2 * math.pi * i / 5),
            distance * math.sin(2 * math.pi * i / 5),
            0.0,
        )
        for i in range(5)
    ]

    axial = [
        (0.0, 0.0, distance),
        (0.0, 0.0, -distance),
    ]

    return equatorial + axial

def capped_octahedral_positions(distance: float = 2.0) -> list[Position]:
    """Return seven positions for an approximate capped octahedral geometry."""
    base = octahedral_positions(distance)

    cap = (
        distance / (3 ** 0.5),
        distance / (3 ** 0.5),
        distance / (3 ** 0.5),
    )

    return base + [cap]

def square_pyramidal_positions(distance: float = 2.0) -> list[Position]:
    """Return five positions for a square pyramidal geometry."""
    return [
        ( distance,  0.0, 0.0),
        (-distance,  0.0, 0.0),
        ( 0.0,  distance, 0.0),
        ( 0.0, -distance, 0.0),
        ( 0.0,  0.0, distance),
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
    "square antiprismatic": square_antiprismatic_positions,
    "trigonal bipyramidal": trigonal_bipyramidal_positions,
    "square pyramidal": square_pyramidal_positions,
    "pentagonal bipyramidal": pentagonal_bipyramidal_positions,
    "capped octahedral": capped_octahedral_positions,   
    "octahedral": octahedral_positions,
    "square antiprismatic": square_antiprismatic_positions,
    "antiprismatic": square_antiprismatic_positions,
    "dodecahedral":dodecahedral_positions,
   
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

    # if the geometry is unknown the octahedral one is adopted by default. 
    # For n ligands, if n is < to the number of ligand in the octahedral complex only the n first positions are attributed. 
    base = octahedral_positions(distance)
    if n <= len(base):
        return base[:n]
    return base + [(0.0, 0.0, 0.0)] * (n - len(base))
    #if n > o the number of ligand in the octahedral complex the 2 extra ligands fall in the center (0.0, 0.0, 0.0)


def bidentate_site_pair_indices(
    geometry: str | None,
    n_sites: int,
    *,
    short_bidentate: bool = True,
) -> list[tuple[int, int]]:
    """Return chemically reasonable site pairs for bidentate ligands."""
    if not short_bidentate:
        return [(i, i + 1) for i in range(0, n_sites - 1, 2)]

    normalized = (geometry or "").lower()

    if "octahedral" in normalized and n_sites >= 6:
        return [(0, 2), (1, 4), (3, 5)]

    if "square planar" in normalized and n_sites >= 4:
        return [(0, 2), (1, 3)]

    if "pentagonal bipyramidal" in normalized and n_sites >= 7:
        return [(3, 4), (5, 6), (2, 3), (6, 2)]

    if "square antiprismatic" in normalized and n_sites >= 8:
        return [(0, 1), (2, 3), (4, 5), (6, 7)]

    return [(i, i + 1) for i in range(0, n_sites - 1, 2)]


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


def find_donor_atom(mol, donor_symbol, override: Optional[int] = None):
    if override is not None and 0 <= override < mol.GetNumAtoms():
        return override

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == donor_symbol:
            return atom.GetIdx()

    raise ValueError(f"No donor atom {donor_symbol} found in ligand")

        
 #function for the ligands disposition around the metal


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
        target_mid[2] + outward[2] * push)
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

def translate_tetradentate_ligand(ligand_mol, donor_indices, target_sites):
    ligand_conf = ligand_mol.GetConformer()

    donor_positions=[]
    for idx in donor_indices:
        pos=ligand_conf.GetAtomPosition(idx)
        donor_positions.append((pos.x,pos.y,pos.z))

    ligand_mid=(sum(pos[0] for pos in donor_positions)/4,
                sum(pos[1] for pos in donor_positions)/4,
                sum(pos[2] for pos in donor_positions)/4
                    )
    target_mid=(sum(pos[0] for pos in target_sites)/4,
                sum(pos[1] for pos in target_sites)/4,
                sum(pos[2] for pos in target_sites)/4
                    )
    
    ligand_x=unit(vec_sub(donor_positions[1], donor_positions[0]))
    ligand_tmp=unit(vec_sub(donor_positions[2], donor_positions[0]))
    ligand_z=unit(cross(ligand_x, ligand_tmp))
    ligand_y=unit(cross(ligand_z, ligand_x))


    target_x = unit(vec_sub(target_sites[1], target_sites[0]))
    target_tmp = unit(vec_sub(target_sites[2], target_sites[0]))
    target_z = unit(cross(target_x, target_tmp))
    target_y = unit(cross(target_z, target_x))

    coords={}
    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        old_pos = ligand_conf.GetAtomPosition(idx)
        old_pos = (old_pos.x, old_pos.y, old_pos.z)

        centered = vec_sub(old_pos, ligand_mid)

        x = dot(centered, ligand_x)
        y = dot(centered, ligand_y)
        z = dot(centered, ligand_z)

        final_pos = vec_add(
            target_mid,
            vec_add(
                vec_add(vec_scale(target_x, x), vec_scale(target_y, y)),
                vec_scale(target_z, z),
            ),
        )

        coords[idx] = final_pos

    return coords

def translate_hexadentate_ligand(ligand_mol, donor_indices, target_sites):
    ligand_conf = ligand_mol.GetConformer()

    donor_positions = []
    for idx in donor_indices:
        pos = ligand_conf.GetAtomPosition(idx)
        donor_positions.append((pos.x, pos.y, pos.z))

    ligand_mid = (
        sum(pos[0] for pos in donor_positions) / 6,
        sum(pos[1] for pos in donor_positions) / 6,
        sum(pos[2] for pos in donor_positions) / 6,
    )

    target_mid = (
        sum(pos[0] for pos in target_sites) / 6,
        sum(pos[1] for pos in target_sites) / 6,
        sum(pos[2] for pos in target_sites) / 6,
    )
    donor_shift={}
    for donor_idx, old_donor_pos, target_pos in zip(donor_indices, donor_positions,target_sites):
        donor_shift[donor_idx]=vec_sub(target_pos,old_donor_pos)
    coords = {}

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        old_pos = ligand_conf.GetAtomPosition(idx)
        old_pos = (old_pos.x, old_pos.y, old_pos.z)
        if idx in donor_shift:
            coords[idx]=vec_add(old_pos,donor_shift[idx])
        else:
            centered = vec_sub(old_pos, ligand_mid)
            final_pos = vec_add(centered, target_mid)
            coords[idx] = final_pos

    return coords


# ---------------------------------------------------------------------------
# Complex-level builder
# ---------------------------------------------------------------------------

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

def methyl_ligand_positions(ligand_mol, donor_idx, target, nh_distance=1.0, spread=0.8):
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
                coords[idx]=vec_add(base, vec_scale(outward,0.45))
            else:
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
def oxalate_ligand_positions(ligand_mol, donor_indices, target_sites):
    coords={}
    o1_idx, o2_idx = donor_indices[:2]
    o1_pos, o2_pos = target_sites

    coords[o1_idx] = o1_pos
    coords[o2_idx] = o2_pos
    axis=unit(vec_sub(o2_pos, o1_pos))
    mid = (
        (o1_pos[0] + o2_pos[0]) / 2,
        (o1_pos[1] + o2_pos[1]) / 2,
        (o1_pos[2] + o2_pos[2]) / 2,
    )
    outward=unit(mid)

    if norm(outward) == 0:
        outward = (0.0, 0.0, 1.0)

    side = cross(axis, outward)

    if norm(side) == 0:
        side = (0.0, 0.0, 1.0)

    side = unit(side)
    o1_atom=ligand_mol.GetAtomWithIdx(o1_idx)
    o2_atom=ligand_mol.GetAtomWithIdx(o2_idx)

    c1_candidates = [
    nbr.GetIdx()
    for nbr in o1_atom.GetNeighbors()
    if nbr.GetSymbol() == "C"
]

    c2_candidates = [
    nbr.GetIdx()
    for nbr in o2_atom.GetNeighbors()
    if nbr.GetSymbol() == "C"
]
    c1_idx = c1_candidates[0]
    c2_idx = c2_candidates[0]

    if c1_candidates and c2_candidates:
        coords[c1_idx] = vec_add(o1_pos,vec_add(vec_scale(axis,0.45),vec_scale(outward,0.7)))
        coords[c2_idx] = vec_add(o2_pos,vec_add(vec_scale(axis,-0.45),vec_scale(outward,0.7)))

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()

        if idx in coords:
            continue

        if atom.GetSymbol() == "O":
            neighbors = atom.GetNeighbors()
            if neighbors:
                c_idx = neighbors[0].GetIdx()
                base=coords.get(c_idx,mid)
                local_out = unit(vec_sub(base, mid))

                if norm(local_out) == 0:
                    local_out = outward

                coords[idx] = vec_add(
                    base,
                    vec_add(
                        vec_scale(local_out, 0.35),
                        vec_scale(side, 0.25 if idx % 2 == 0 else -0.25),
                    ),
                )
            else:
                coords[idx] = mid
        else:
            coords[idx] = mid

    return coords

def ethylenediamine_ligand_positions(ligand_mol, donor_indices, target_sites):
    coords = {}

    n1_idx, n2_idx = donor_indices[:2]
    n1_pos, n2_pos = target_sites

    coords[n1_idx] = n1_pos
    coords[n2_idx] = n2_pos

    mid = (
        (n1_pos[0] + n2_pos[0]) / 2,
        (n1_pos[1] + n2_pos[1]) / 2,
        (n1_pos[2] + n2_pos[2]) / 2,
    )

    axis = unit(vec_sub(n2_pos, n1_pos))

    outward = unit(mid)
    if norm(outward) == 0:
        outward = (0.0, 0.0, 1.0)

    side = cross(axis, outward)
    if norm(side) == 0:
        side = (0.0, 0.0, 1.0)
    side = unit(side)

    n1_atom = ligand_mol.GetAtomWithIdx(n1_idx)
    n2_atom = ligand_mol.GetAtomWithIdx(n2_idx)

    c1_candidates = [
        nbr.GetIdx()
        for nbr in n1_atom.GetNeighbors()
        if nbr.GetSymbol() == "C"
    ]

    c2_candidates = [
        nbr.GetIdx()
        for nbr in n2_atom.GetNeighbors()
        if nbr.GetSymbol() == "C"
    ]

    if c1_candidates and c2_candidates:
        c1_idx = c1_candidates[0]
        c2_idx = c2_candidates[0]

        coords[c1_idx] = vec_add(
            mid,
            vec_add(
                vec_scale(axis, -0.45),
                vec_scale(outward, 1.2),
            ),
        )

        coords[c2_idx] = vec_add(
            mid,
            vec_add(
                vec_scale(axis, 0.45),
                vec_scale(outward, 1.2),
            ),
        )

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()

        if idx in coords:
            continue

        if atom.GetSymbol() == "H":
            neighbors = atom.GetNeighbors()
            if neighbors:
                base_idx = neighbors[0].GetIdx()
                base = coords.get(base_idx, mid)

                h_out = unit(vec_sub(base, mid))
                if norm(h_out) == 0:
                    h_out = outward

                coords[idx] = vec_add(
                    base,
                    vec_add(
                        vec_scale(h_out, 0.35),
                        vec_scale(side, 0.25 if idx % 2 == 0 else -0.25),
                    ),
                )
            else:
                coords[idx] = mid
        else:
            coords[idx] = mid

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

def bipyridine_ligand_positions(ligand_mol, donor_indices, target_sites, ligand_number=0):
    coords={}
    n1_idx, n2_idx = donor_indices[:2]
    n1_pos, n2_pos = target_sites

    coords[n1_idx] = n1_pos
    coords[n2_idx] = n2_pos


    n1_atom=ligand_mol.GetAtomWithIdx(n1_idx)
    n2_atom=ligand_mol.GetAtomWithIdx(n2_idx)
    axis=unit(vec_sub(n2_pos, n1_pos))
    mid = (
        (n1_pos[0] + n2_pos[0]) / 2,
        (n1_pos[1] + n2_pos[1]) / 2,
        (n1_pos[2] + n2_pos[2]) / 2,
    )
    outward=unit(mid)

    if norm(outward) == 0:
        outward = (0.0, 0.0, 1.0)

    side = cross(axis, outward)

    if norm(side) == 0:
        side = (0.0, 0.0, 1.0)

    side = unit(side)


    if ligand_number % 2 == 1:
        side = vec_scale(side, -1.0)

    ring_info = ligand_mol.GetRingInfo()
    rings = [list(ring) for ring in ring_info.AtomRings()]

    donor_rings = []
    for donor_idx in donor_indices:
        for ring in rings:
            if donor_idx in ring:
                donor_rings.append((donor_idx, ring))
                break

    ring_radius = 0.75

    for donor_idx, ring in donor_rings:
        donor_pos = coords[donor_idx]
        ring_center = vec_add(donor_pos, vec_scale(outward, ring_radius))

        ordered_ring = [donor_idx] + [idx for idx in ring if idx != donor_idx]

        angles = [
            math.pi,
            2 * math.pi / 3,
            math.pi / 3,
            0.0,
            -math.pi / 3,
            -2 * math.pi / 3,
        ]

        for idx, angle in zip(ordered_ring, angles):
            ring_pos = vec_add(
                vec_scale(outward, math.cos(angle) * ring_radius),
                vec_scale(side, math.sin(angle) * ring_radius),
            )
            coords[idx] = vec_add(ring_center, ring_pos)

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()

        if idx in coords:
            continue

        if atom.GetSymbol() == "H":
            neighbors = atom.GetNeighbors()
            if neighbors:
                base_idx = neighbors[0].GetIdx()
                base = coords.get(base_idx, mid)

                h_out = unit(vec_sub(base, mid))

                if norm(h_out) == 0:
                    h_out = outward

                coords[idx] = vec_add(
                    base,
                        vec_scale(h_out, 0.45),
                )
            else:
                coords[idx] = mid
        else:
            coords[idx] = mid

    return coords

def water_ligand_positions(ligand_mol, donor_idx, target):
    direction = unit(target)
    coords = {}

    coords[donor_idx] = target

    if abs(direction[0]) < 0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    side = unit(cross(direction, ref))

    h_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "H"
    ]

    angle = 0.55
    oh_distance = 0.75
    outward = vec_scale(direction, oh_distance)

    if len(h_indices) >= 2:
        coords[h_indices[0]] = vec_add(
            target,
            vec_add(outward, vec_scale(side, angle)),
        )

        coords[h_indices[1]] = vec_add(
            target,
            vec_add(outward, vec_scale(side, -angle)),
        )

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in coords:
            coords[idx] = target

    return coords

def dmso_ligand_positions(ligand_mol, donor_idx, target, ligand_number=0):
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

    angle_offset= ligand_number * math.pi

    s_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "S"]
    o_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "O"]
    c_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "C"]
    s_idx=s_indices[0]
    o_idx=o_indices[0]

    if donor_idx==s_idx:
        coords[s_idx]=target
        s_center=target
        o_pos=vec_add(s_center, vec_scale(direction, 1.2))
        coords[o_idx]=o_pos

        for i, c_idx in enumerate(c_indices):
            angle = angle_offset + 2 * math.pi * i / 2

            side = vec_add(
                vec_scale(u, math.cos(angle) * 1.0),
                vec_scale(v, math.sin(angle) * 1.0),
            )

            c_pos = vec_add(
                s_center,
                vec_add(vec_scale(direction, 0.7), side),
            )
            coords[c_idx] = c_pos

    elif donor_idx==o_idx:
        coords[o_idx]=target
        o_center=target
        s_pos=vec_add(o_center, vec_scale(direction, 1.2))
        coords[s_idx]=s_pos

        for i, c_idx in enumerate(c_indices):
            angle = angle_offset + 2 * math.pi * i / 2

            side = vec_add(
                vec_scale(u, math.cos(angle) * 1.0),
                vec_scale(v, math.sin(angle) * 1.0),
            )

            c_pos = vec_add(
                s_pos,
                vec_add(vec_scale(direction, 0.6), side),
            )
            coords[c_idx] = c_pos
    for atom in ligand_mol.GetAtoms():
        idx= atom.GetIdx()
        if idx in coords:
            continue
        if atom.GetSymbol() == "H":
            neighbors = atom.GetNeighbors()
            if neighbors:
                c_idx = neighbors[0].GetIdx()
                base = coords.get(c_idx, target)

                angle = angle_offset + idx * 2.1
                h_side = vec_add(
                    vec_scale(u, math.cos(angle) * 0.45),
                    vec_scale(v, math.sin(angle) * 0.45),
                )
                h_out = vec_scale(direction, 0.35)

                coords[idx] = vec_add(base, vec_add(h_side, h_out))
            else:
                coords[idx] = target
        else:
            coords[idx] = target
    return coords

def cyclopentadienyl_ligand_positions(ligand_mol,target,ligand_number=0):
    direction=unit(target)
    if abs(direction[0])<0.9:
        ref = (1.0, 0.0, 0.0)
    else:
        ref = (0.0, 1.0, 0.0)

    u = unit(cross(direction, ref))
    v = unit(cross(direction, u))

    coords = {}
    ring_center=target

    c_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "C"]
    h_indices = [
        atom.GetIdx()
        for atom in ligand_mol.GetAtoms()
        if atom.GetSymbol() == "H"]
    
    angle_offset = ligand_number * math.pi / 5
    ring_radius = 0.85

    for i, c_idx in enumerate(c_indices[:5]):
        angle = angle_offset + 2 * math.pi * i / 5

        ring_pos = vec_add(
            vec_scale(u, math.cos(angle) * ring_radius),
            vec_scale(v, math.sin(angle) * ring_radius),
        )

        coords[c_idx] = vec_add(ring_center, ring_pos)

    for h_idx in h_indices:
        atom = ligand_mol.GetAtomWithIdx(h_idx)
        neighbors = atom.GetNeighbors()

        if neighbors:
            c_idx = neighbors[0].GetIdx()
            base = coords.get(c_idx, ring_center)

            outward = unit(vec_sub(base, ring_center))
            coords[h_idx] = vec_add(base, vec_scale(outward, 0.45))
        else:
            coords[h_idx] = ring_center

    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx not in coords:
            coords[idx] = ring_center

    return coords

def edta_ligand_positions(ligand_mol,donor_indices, target_sites):
    ligand_conf=ligand_mol.GetConformer()
    coords={}
#place donors
    for donor_idx, target in zip(donor_indices, target_sites):
        coords[donor_idx]=target
    donor_set=set(donor_indices)
    n_donors = [
    idx for idx in donor_indices
    if ligand_mol.GetAtomWithIdx(idx).GetSymbol() == "N"
]

    if len(n_donors) >= 2:
        n1_idx, n2_idx = n_donors[:2]
        n1_pos = coords[n1_idx]
        n2_pos = coords[n2_idx]

        axis = unit(vec_sub(n2_pos, n1_pos))

        if abs(axis[0]) < 0.9:
            ref = (1.0, 0.0, 0.0)
        else:
            ref = (0.0, 1.0, 0.0)

        side = unit(cross(axis, ref))
    bridge_carbons = []

    for bond in ligand_mol.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()

            if a.GetSymbol() != "C" or b.GetSymbol() != "C":
                continue

            a_neighbors = [n.GetIdx() for n in a.GetNeighbors()]
            b_neighbors = [n.GetIdx() for n in b.GetNeighbors()]

            a_connected_to_n = any(n in n_donors for n in a_neighbors)
            b_connected_to_n = any(n in n_donors for n in b_neighbors)

            if a_connected_to_n and b_connected_to_n:
                bridge_carbons = [a.GetIdx(), b.GetIdx()]
                break

    if len(bridge_carbons) == 2:
            c1_idx, c2_idx = bridge_carbons

            coords[c1_idx] = vec_add(
    n1_pos,
    vec_add(
        vec_scale(axis, 0.65),
        vec_scale(side, 1.4),
    ),
)

            coords[c2_idx] = vec_add(
    n2_pos,
    vec_add(
        vec_scale(axis, -0.65),
        vec_scale(side, 1.4),
    ),
)

#place N-C 
    for donor_idx in donor_indices:
        donor_atom = ligand_mol.GetAtomWithIdx(donor_idx)

        if donor_atom.GetSymbol() != "N":
            continue

        donor_pos = coords[donor_idx]
        direction = unit(donor_pos)

        if abs(direction[0]) < 0.9:
            ref = (1.0, 0.0, 0.0)
        else:
            ref = (0.0, 1.0, 0.0)

        side = unit(cross(direction, ref))

        carbon_neighbors = [
        n.GetIdx()
        for n in donor_atom.GetNeighbors()
        if n.GetSymbol() == "C"
    ]
        for i,c_idx in enumerate(carbon_neighbors):
            if c_idx in coords:
                continue
            sign=1 if i%2==0 else -1
            coords[c_idx]=vec_add(donor_pos,vec_add(vec_scale(direction,0.65),vec_scale(side,0.75*sign)))
#place C-O           
    for donor_idx in donor_indices:
        donor_atom = ligand_mol.GetAtomWithIdx(donor_idx)

        if donor_atom.GetSymbol() != "O":
            continue

        donor_pos = coords[donor_idx]
        direction=unit(donor_pos)

        if abs(direction[0]) < 0.9:
            ref = (1.0, 0.0, 0.0)
        else:
            ref = (0.0, 1.0, 0.0)

        side = unit(cross(direction, ref))

        carbonyl_carbons = [
        n.GetIdx()
        for n in donor_atom.GetNeighbors()
        if n.GetSymbol() == "C"
    ]

        if not carbonyl_carbons:
            continue

        c_idx = carbonyl_carbons[0]
        c_pos = vec_add(donor_pos, vec_scale(direction, 0.55))
        coords[c_idx] = c_pos

        carbon_atom = ligand_mol.GetAtomWithIdx(c_idx)

        other_oxygens = [
        n.GetIdx()
        for n in carbon_atom.GetNeighbors()
        if n.GetSymbol() == "O" and n.GetIdx() != donor_idx
    ]

        for i, o_idx in enumerate(other_oxygens):
            sign = 1 if i % 2 == 0 else -1
            coords[o_idx] = vec_add(
            c_pos,
            vec_add(
                vec_scale(direction, 0.35),
                vec_scale(side, 0.25 * sign),
            ),
        )

    donor_set= set(donor_indices)
#place other heavy atoms
    for atom in ligand_mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in coords:
            continue
        if atom.GetSymbol()== "H":
            continue
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        connected_donors = [n for n in neighbors if n in donor_set]

        if connected_donors:
            donor_idx = connected_donors[0]
            donor_pos = coords[donor_idx]
            outward = unit(donor_pos)

            coords[idx] = vec_add(
                donor_pos,
                vec_scale(outward, 1.0),
            )
        else:
            old_pos = ligand_conf.GetAtomPosition(idx)
            old_pos = (old_pos.x, old_pos.y, old_pos.z)

            outward = unit(old_pos)
            coords[idx] = vec_scale(outward, max(norm(old_pos), 1.5))
#place hydrogens
    for atom in ligand_mol.GetAtoms():
        if atom.GetSymbol() != "H":
            continue
        h_idx=atom.GetIdx()
        neighbors = atom.GetNeighbors()
        if len(neighbors)==0:
            continue
        neighbor_idx=neighbors[0].GetIdx()
        if neighbor_idx not in coords:
            old_pos = ligand_conf.GetAtomPosition(h_idx)
            coords[h_idx] = (old_pos.x, old_pos.y, old_pos.z)
            continue
        base_pos=coords[neighbor_idx]
        outward=unit(base_pos)

        coords[h_idx] = vec_add(
            base_pos,
            vec_scale(outward, 0.75),
            )

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
    occupied=[False]*len(sites)
    tridentate_count=0
    def next_free_site():
        for i, is_occupied in enumerate(occupied):
            if not is_occupied:
                occupied[i] = True
                return i
        return None


    def next_free_pair(ligand_symbol: str, denticity: int):
        pairs = bidentate_site_pair_indices(
            geometry,
            len(sites),
            short_bidentate=is_short_bidentate_ligand(
                ligand_symbol,
                denticity,
            ),
        )

        for a, b in pairs:
            if not occupied[a] and not occupied[b]:
                occupied[a] = True
                occupied[b] = True
                return a, b

        return None


    def next_free_triplet():
        if "octahedral" in geometry.lower() and len(sites) >= 6:
            triplets = [(0, 2, 1), (3, 4, 5)]
        else:
            triplets = [(i, i + 1, i + 2) for i in range(0, len(sites) - 2, 3)]

        for a, b, c in triplets:
            if not occupied[a] and not occupied[b] and not occupied[c]:
                occupied[a] = True
                occupied[b] = True
                occupied[c] = True
                return a, b, c

        return None
    def next_free_quadruplet():
        if "octahedral" in geometry.lower() and len(sites)>=6:
            quadruplets=[(0,1,2,3)]
        else:
            quadruplets=[(i, i + 1, i + 2, i + 3)
            for i in range(0, len(sites) - 3, 4)]
        for a, b, c, d in quadruplets:
            if not occupied[a] and not occupied[b] and not occupied[c] and not occupied[d]:
                occupied[a] = True
                occupied[b] = True
                occupied[c] = True
                occupied[d] = True
                return a, b, c, d

        return None
    def next_free_sextuplet():
        if "octahedral" in geometry.lower() and len(sites)>=6:
            sextuplets=[(0,1,2,3,4,5)]
        else:
            sextuplets=[(i, i + 1, i + 2, i + 3,i+4,i+5)
            for i in range(0, len(sites) - 5, 6)]
        for a, b, c, d,e,f in sextuplets:
            if not occupied[a] and not occupied[b] and not occupied[c] and not occupied[d] and not occupied[e] and not occupied[f]:
                occupied[a] = True
                occupied[b] = True
                occupied[c] = True
                occupied[d] = True
                occupied[e] = True
                occupied[f] = True
                return a, b, c, d, e, f

        return None
    ligand_items = sorted(
    parsed.ligands.items(),
    key=lambda item: parsed.ligand_denticity.get(item[0], 1),
    reverse=True)
    
    for ligand_symbol, count in ligand_items:
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
        donor_overrides = donor_index_overrides_for_ligand(
            ligand_symbol,
            donor_symbol,
        )
        denticity = parsed.ligand_denticity.get(ligand_symbol, 1)

        for _ in range(count):
            if site_index >= len(sites):
                break

            try:
                ligand_mol = build_ligand_3d(smiles)
            except ValueError:
                continue
            #bidentate case
            if denticity == 2 and len(donor_overrides) >= 2:
                bidentate_index = site_index // 2
                
                pair=next_free_pair(ligand_symbol, denticity)
                if pair is None:
                    break
                target_sites = (sites[pair[0]],sites[pair[1]])
                donor_indices = donor_overrides[:2]

                if ligand_symbol=="en":
                    local_coords=ethylenediamine_ligand_positions(ligand_mol, donor_indices,target_sites)

                elif ligand_symbol=="bipy":
                    local_coords=bipyridine_ligand_positions(ligand_mol, donor_indices,target_sites, bidentate_index)
                elif ligand_symbol=="ox":
                    local_coords=oxalate_ligand_positions(ligand_mol, donor_indices,target_sites)
                else:
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
            #tridentate case
            elif denticity==3 and len(donor_overrides)>=3:
                triplet = next_free_triplet()
                if triplet is None:
                    break       
                target_sites = (sites[triplet[0]], sites[triplet[1]], sites[triplet[2]])
                tridentate_index = 0
                donor_indices = donor_overrides[:3]
                if ligand_symbol in ("terpy"):
                    local_coords = terpyridine_ligand_positions(
        ligand_mol,
        donor_overrides[:3],
        tridentate_count,
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
            #tetradentate case
            elif denticity==4 and len(donor_overrides)>=4:
                quadruplet=next_free_quadruplet()
                if quadruplet is None:
                    break       
                target_sites = (sites[quadruplet[0]], sites[quadruplet[1]],sites[quadruplet[2]],sites[quadruplet[3]])
                tridentate_index = 0
                donor_indices=donor_overrides[:4]

                local_coords=translate_tetradentate_ligand(ligand_mol,donor_indices,target_sites)
                offset=rw.GetNumAtoms()
                for atom in ligand_mol.GetAtoms():
                    global_idx = rw.AddAtom(atom)
                    coords[global_idx] = local_coords[atom.GetIdx()]

                for bond in ligand_mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx() + offset
                    a2 = bond.GetEndAtomIdx() + offset
                    rw.AddBond(a1, a2, bond.GetBondType())

                for donor_idx in donor_indices:
                    rw.AddBond(metal_idx,donor_idx + offset,Chem.BondType.DATIVE)

                site_index += 4
                continue

            #pentadentate case
            elif ligand_symbol == "Cp" or denticity==5:
                    site=next_free_site()
                    if site is None:
                        break
                    target=sites[site]
                    if ligand_symbol == "Cp":

                        local_coords=cyclopentadienyl_ligand_positions(ligand_mol,target,site)
                    else:
                        ligand_conf=ligand_mol.GetConformer()
                        local_coords={}
                        for atom in ligand_mol.GetAtoms():
                            idx=atom.GetIdx()
                            old_pos=ligand_conf.GetAtomPosition(idx)
                            local_coords[idx]=(old_pos.x+target[0],old_pos.y+target[1],old_pos.z+target[2])
                    offset = rw.GetNumAtoms()

                    for atom in ligand_mol.GetAtoms():
                        global_idx = rw.AddAtom(atom)
                        coords[global_idx] = local_coords[atom.GetIdx()]

                    for bond in ligand_mol.GetBonds():
                        a1 = bond.GetBeginAtomIdx() + offset
                        a2 = bond.GetEndAtomIdx() + offset
                        rw.AddBond(a1, a2, bond.GetBondType())
                    if ligand_symbol == "Cp":
                        #c_indices = [atom.GetIdx()for atom in ligand_mol.GetAtoms() if atom.GetSymbol() == "C" ]
                        dummy = Chem.Atom("He")
                        dummy_idx = rw.AddAtom(dummy)

                        coords[dummy_idx] = target

                        rw.AddBond( metal_idx,dummy_idx,Chem.BondType.DATIVE)

                        #for c_idx in c_indices[:5]:
                            #rw.AddBond(metal_idx,  c_idx + offset,  Chem.BondType.DATIVE )

                    site_index+=1
                    continue
            #else:
               # donor_idx_local = find_donor_atom(ligand_mol, donor_symbol)
               # rw.AddBond(
           # metal_idx, donor_idx_local + offset, Chem.BondType.DATIVE )
               # site_index+= 1
               # continue 
            elif denticity==6 and len(donor_overrides)>=6:

                if len(sites)<6:
                    break
                sextuplet = next_free_sextuplet()
                if sextuplet is None:
                    break

                
                donor_indices = donor_overrides[:6]
                if ligand_symbol in ("EDTA", "edta"):
                   target_sites = (
    sites[0],  
    sites[2], 
    sites[4], 

    sites[1],  
    sites[3], 
    sites[5],  
)

                   local_coords=edta_ligand_positions(ligand_mol,donor_indices, target_sites)
                else:
                    target_sites = (sites[sextuplet[0]],sites[sextuplet[1]],sites[sextuplet[2]],sites[sextuplet[3]],sites[sextuplet[4]], sites[sextuplet[5]])
                    local_coords = translate_hexadentate_ligand(ligand_mol,donor_indices,target_sites)

                offset = rw.GetNumAtoms()

                for atom in ligand_mol.GetAtoms():
                    global_idx = rw.AddAtom(atom)
                    coords[global_idx] = local_coords[atom.GetIdx()]

                for bond in ligand_mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx() + offset
                    a2 = bond.GetEndAtomIdx() + offset
                    rw.AddBond(a1, a2, bond.GetBondType())

                for donor_idx in donor_indices:
                    rw.AddBond(metal_idx,donor_idx + offset,Chem.BondType.DATIVE)

                site_index += 6
                continue

            primary_override = donor_overrides[0] if donor_overrides else None

            if primary_override is not None and primary_override < ligand_mol.GetNumAtoms():
                    donor_idx_local = primary_override
            else:
                try:
                    donor_idx_local = find_donor_atom(ligand_mol, donor_symbol)
                except ValueError:
                    continue
        
            #monodentate case
            site = next_free_site()
            if site is None:
                break
            target = sites[site]

            if ligand_symbol=="NH3":
                local_coords=ammonia_ligand_positions(ligand_mol, donor_idx_local,target)
            elif ligand_symbol=="py":
                local_coords=pyridine_ligand_positions(ligand_mol, donor_idx_local,target, site)    
            elif ligand_symbol=="CH3":
                local_coords=methyl_ligand_positions(ligand_mol, donor_idx_local,target, site)    
            elif ligand_symbol=="PMe3":
                local_coords=trimethyl_ligand_positions(ligand_mol, donor_idx_local,target)
            elif ligand_symbol=="dmso":
                local_coords=dmso_ligand_positions(ligand_mol, donor_idx_local,target,site)
            elif ligand_symbol=="PEt3":
                local_coords=triethyl_ligand_positions(ligand_mol,donor_idx_local,target)
            elif ligand_symbol in ("NO2","ONO"):
                local_coords=nitrito_ligand_positions(ligand_mol, donor_idx_local,target,ligand_symbol)
            elif ligand_symbol=="H2O":
                local_coords=water_ligand_positions(ligand_mol, donor_idx_local,target)
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

   



#functions to display diatomic ligands correctly
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

def parse_complex_input(complex_input: str | ParsedComplex) -> ParsedComplex:
    """Accept a formula string, compound name, or already parsed complex."""
    if isinstance(complex_input, ParsedComplex):
        return complex_input
    if isinstance(complex_input, str):
        try:
            return parse_formula(complex_input)
        except Exception:
            return parse_name(complex_input)
    parsed = getattr(complex_input, "parsed", None)
    if isinstance(parsed, ParsedComplex):
        return parsed
    raise TypeError(
        "Expected a ParsedComplex, Complex, or formula/name string"
    )


def _to_parsed(complex_or_formula) -> ParsedComplex:
    """Best-effort coercion to a ``ParsedComplex``."""
    return parse_complex_input(complex_or_formula)

#function to display the molecule on a notebook

def view_complex_3d(
    complex_or_formula,
    width: int = 400,
    height: int = 400,
    distance: float = 2.0,
    geometry: str | None = None,
):
    import py3Dmol 

    parsed = _to_parsed(complex_or_formula)
    mol = build_complex_3d(parsed, distance=distance, geometry=geometry)

    block = Chem.MolToMolBlock(mol, kekulize=False) #pour ne pas forcer la conversion des cycles aromatiques
    view = py3Dmol.view(width=width, height=height)
    view.removeAllModels()
    view.addModel(block, "sdf")
    view.setStyle({}, {"stick": {}, "sphere": {"scale": 0.25}})
    view.setStyle({}, {"stick": {"radius": 0.12}, "sphere": {"scale": 0.18}})
    view.setStyle({"elem": "He"}, {"sphere": {"scale": 0.01, "color":"white"}, "stick": {"radius": 0.03,"color":"grey"}})



    view.setBackgroundColor("white")
    view.zoomTo()
    return view

#useful function to display the complex on an app such as streamlit
def complex_3d_html(
    complex_or_formula,
    width: int = 400,
    height: int = 400,
    distance: float = 2.0,
    geometry: str | None = None,
) -> str:
    view = view_complex_3d(
        complex_or_formula,
        width=width,
        height=height,
        distance=distance,
        geometry=geometry,
    )
    return view._make_html()
