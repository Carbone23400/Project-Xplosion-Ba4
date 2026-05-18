"""
coordchem/viz/layout_2d.py
--------------------------
2D coordination-site layouts for coordination complex depictions.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, sin

from ..parser import ParsedComplex
from .ligand_data import is_short_bidentate_ligand


@dataclass(frozen=True)
class Site:
    """One coordination site around the metal."""

    x: float
    y: float
    style: str  # plain, wedge, dash


def _polar_site(angle_degrees: float, radius: float, style: str) -> Site:
    """Build one site from polar coordinates."""
    angle = angle_degrees * pi / 180
    return Site(radius * cos(angle), radius * sin(angle), style)


def _cn4_depth_sites(n_sites: int) -> list[Site]:
    """Return tetrahedral CN=4 sites with front/back vertical depth cues."""
    return [
        Site(-3.0, 0.0, "plain"),
        Site(0.0, 2.85, "dash"),
        Site(3.0, 0.0, "plain"),
        Site(0.0, -2.85, "wedge"),
    ][:n_sites]


def _square_planar_sites(n_sites: int) -> list[Site]:
    """Return flat square-planar CN=4 sites as a cross."""
    return [
        Site(0.0, 3.2, "plain"),
        Site(-3.2, 0.0, "plain"),
        Site(3.2, 0.0, "plain"),
        Site(0.0, -3.2, "plain"),
    ][:n_sites]


def _cn6_depth_sites(n_sites: int) -> list[Site]:
    """Return six sites with 2 plain, 2 wedge, and 2 dash bonds."""
    return [
        _polar_site(90, 3.4, "plain"),
        _polar_site(30, 3.4, "dash"),
        _polar_site(-30, 3.4, "wedge"),
        _polar_site(-90, 3.4, "plain"),
        _polar_site(-150, 3.4, "wedge"),
        _polar_site(150, 3.4, "dash"),
    ][:n_sites]


def _square_antiprismatic_sites(n_sites: int) -> list[Site]:
    """Return a staggered two-square projection for square-antiprismatic CN=8."""
    return [
        # Upper square: two thin bonds, one front wedge, one back dash.
        Site(-2.35, 2.55, "plain"),
        Site(0.55, 3.15, "dash"),
        Site(2.35, 2.25, "plain"),
        Site(-0.55, 1.65, "wedge"),
        # Lower square, shifted in quincunx: two dashed back bonds and two wedges.
        Site(1.25, -1.25, "dash"),
        Site(-1.25, -1.25, "dash"),
        Site(-2.65, -2.45, "wedge"),
        Site(0.85, -2.45, "wedge"),
    ][:n_sites]


def _pentagonal_bipyramidal_sites(n_sites: int) -> list[Site]:
    """Return a styled pentagonal-bipyramidal projection for CN=7."""
    return [
        Site(0.0, 4.2, "plain"),
        Site(0.0, -4.2, "plain"),
        _polar_site(0, 3.35, "plain"),
        _polar_site(72, 2.75, "dash"),
        _polar_site(144, 3.35, "dash"),
        _polar_site(216, 3.35, "wedge"),
        _polar_site(288, 2.75, "wedge"),
    ][:n_sites]


def _capped_octahedral_sites(n_sites: int) -> list[Site]:
    """Return a flat heptagonal projection for capped-octahedral CN=7."""
    return regular_polygon_sites(n_sites, radius=2.9)


def coordination_sites(geometry: str, n_sites: int) -> list[Site]:
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
        return _square_planar_sites(n_sites)

    if g == "tetrahedral":
        return _cn4_depth_sites(n_sites)

    if g == "octahedral":
        return _cn6_depth_sites(n_sites)

    if "square planar" in g:
        return coordination_sites("square planar", n_sites)

    if "tetrahedral" in g:
        return coordination_sites("tetrahedral", n_sites)

    if "trigonal bipyramidal" in g:
        return [
            Site(0.0, -3.2, "plain"),
            Site(0.0, 3.2, "plain"),
            Site(-3.2, 0.0, "plain"),
            Site(2.6, 1.8, "dash"),
            Site(2.6, -1.8, "wedge"),
        ][:n_sites]

    if "square pyramidal" in g:
        return [
            _polar_site(90, 3.4, "plain"),
            _polar_site(30, 3.4, "dash"),
            _polar_site(-30, 3.4, "wedge"),
            _polar_site(-150, 3.4, "wedge"),
            _polar_site(150, 3.4, "dash"),
        ][:n_sites]

    if "pentagonal bipyramidal" in g:
        return _pentagonal_bipyramidal_sites(n_sites)

    if "capped octahedral" in g:
        return _capped_octahedral_sites(n_sites)

    if "square antiprismatic" in g:
        return _square_antiprismatic_sites(n_sites)

    return regular_polygon_sites(n_sites)


def chelate_octahedral_sites(n_ligands: int) -> list[Site]:
    """
    Return spread-out site pairs for tris-bidentate octahedral complexes.

    Designed so that:
    - top and bottom bonds are plain (thin vertical bonds),
    - upper side bonds are dashed,
    - lower side bonds are wedged.
    This matches the common textbook-style projection for [M(en)3]-type complexes.
    """
    if n_ligands != 3:
        return coordination_sites("octahedral", n_ligands * 2)

    return [
        Site(-3.0, 1.7, "dash"),
        Site(0.0, 3.5, "plain"),
        Site(3.0, 1.7, "dash"),
        Site(3.0, -1.7, "wedge"),
        Site(0.0, -3.5, "plain"),
        Site(-3.0, -1.7, "wedge"),
    ]


def tridentate_octahedral_sites(n_ligands: int) -> list[Site]:
    """Return site triplets for bis-tridentate octahedral complexes."""
    if n_ligands != 2:
        return coordination_sites("octahedral", n_ligands * 3)

    radius = 3.4
    angles = [150, 90, 30, -30, -90, -150]
    styles = ["dash", "plain", "dash", "wedge", "plain", "wedge"]

    return [
        Site(
            radius * cos(angle * pi / 180),
            radius * sin(angle * pi / 180),
            style,
        )
        for angle, style in zip(angles, styles)
    ]


def edta_octahedral_sites() -> list[Site]:
    """Return a compact octahedral projection for EDTA-like hexadentate binding."""
    radius = 2.15
    return [
        _polar_site(150, radius, "dash"),
        _polar_site(90, radius, "plain"),
        _polar_site(30, radius, "dash"),
        _polar_site(-150, radius, "wedge"),
        _polar_site(-90, radius, "plain"),
        _polar_site(-30, radius, "wedge"),
    ]


def should_use_chelate_layout(
    parsed: ParsedComplex,
    ligand_items: list[str],
    geometry: str,
) -> bool:
    """Return True when a tris-bidentate octahedral layout will be clearer."""
    if geometry.lower().strip() != "octahedral":
        return False

    if len(ligand_items) != 3:
        return False

    return all(parsed.ligand_denticity.get(ligand, 1) == 2 for ligand in ligand_items)


def should_use_mixed_polydentate_layout(
    parsed: ParsedComplex,
    ligand_items: list[str],
    geometry: str,
) -> bool:
    """Return True when polydentate ligands should reserve sites first."""
    if not ligand_items:
        return False

    denticities = [parsed.ligand_denticity.get(ligand, 1) for ligand in ligand_items]
    has_polydentate = any(denticity > 1 for denticity in denticities)
    has_monodentate = any(denticity == 1 for denticity in denticities)

    if not (has_polydentate and has_monodentate):
        return False

    return True


def mixed_polydentate_site_groups(
    parsed: ParsedComplex,
    ligand_items: list[str],
    geometry: str,
) -> list[tuple[Site, ...]]:
    """Assign polydentate ligand sites first, then fill monodentate sites."""
    denticities = [parsed.ligand_denticity.get(ligand, 1) for ligand in ligand_items]
    n_sites = sum(denticities)
    normalized_geometry = geometry.lower().strip()

    if (
        normalized_geometry == "octahedral"
        and n_sites == 6
        and any(denticity == 3 for denticity in denticities)
    ):
        sites = tridentate_octahedral_sites(2)
    elif (
        normalized_geometry == "octahedral"
        and n_sites == 6
        and any(denticity == 2 for denticity in denticities)
    ):
        sites = chelate_octahedral_sites(3)
    elif (
        "pentagonal bipyramidal" in normalized_geometry
        and n_sites == 7
        and all(denticity in {1, 2} for denticity in denticities)
        and any(denticity == 2 for denticity in denticities)
    ):
        return pentagonal_bipyramidal_mixed_site_groups(ligand_items, denticities)
    elif (
        "capped octahedral" in normalized_geometry
        and n_sites == 7
        and all(denticity in {1, 2} for denticity in denticities)
        and any(denticity == 2 for denticity in denticities)
    ):
        return capped_octahedral_mixed_site_groups(ligand_items, denticities)
    else:
        sites = coordination_sites(geometry, n_sites)

    groups: list[tuple[Site, ...] | None] = [None] * len(ligand_items)
    site_cursor = 0

    ligand_order = [
        index
        for index, denticity in enumerate(denticities)
        if denticity > 1
    ] + [
        index
        for index, denticity in enumerate(denticities)
        if denticity == 1
    ]

    for index in ligand_order:
        denticity = denticities[index]
        groups[index] = tuple(sites[site_cursor: site_cursor + denticity])
        site_cursor += denticity

    return [group for group in groups if group is not None]


def pentagonal_bipyramidal_mixed_site_groups(
    ligand_items: list[str],
    denticities: list[int],
) -> list[tuple[Site, ...]]:
    """Assign CN=7 bidentates to compact styled site pairs."""
    sites = coordination_sites("pentagonal bipyramidal", 7)
    compact_pairs = [
        (Site(1.15, 2.1, "dash"), Site(-1.15, 2.1, "dash")),
        (Site(-1.15, -2.1, "wedge"), Site(1.15, -2.1, "wedge")),
        (Site(3.0, 0.6, "plain"), Site(1.55, 1.95, "dash")),
    ]
    compact_pair_site_indices = [(3, 4), (5, 6), (2, 3)]
    monodentate_sites = [0, 1, 2, 3, 4, 5, 6]
    used: set[int] = set()
    groups: list[tuple[Site, ...] | None] = [None] * len(ligand_items)
    pair_cursor = 0

    for index, (ligand, denticity) in enumerate(zip(ligand_items, denticities)):
        if not is_short_bidentate_ligand(ligand, denticity):
            continue

        if pair_cursor < len(compact_pairs):
            groups[index] = compact_pairs[pair_cursor]
            used.update(compact_pair_site_indices[pair_cursor])
            pair_cursor += 1

    for index, denticity in enumerate(denticities):
        if groups[index] is not None:
            continue

        available = [
            site_index for site_index in monodentate_sites
            if site_index not in used
        ]
        chosen = available[:denticity]
        used.update(chosen)
        groups[index] = tuple(sites[site_index] for site_index in chosen)

    return [group for group in groups if group is not None]


def capped_octahedral_mixed_site_groups(
    ligand_items: list[str],
    denticities: list[int],
) -> list[tuple[Site, ...]]:
    """Assign capped-octahedral bidentates to compact cis-like site pairs."""
    sites = coordination_sites("capped octahedral", 7)
    compact_pairs = [(1, 2), (4, 5), (6, 0)]
    monodentate_sites = [0, 1, 2, 3, 4, 5, 6]
    used: set[int] = set()
    groups: list[tuple[Site, ...] | None] = [None] * len(ligand_items)
    pair_cursor = 0

    for index, (ligand, denticity) in enumerate(zip(ligand_items, denticities)):
        if not is_short_bidentate_ligand(ligand, denticity):
            continue

        while pair_cursor < len(compact_pairs):
            pair = compact_pairs[pair_cursor]
            pair_cursor += 1
            if not any(site_index in used for site_index in pair):
                used.update(pair)
                groups[index] = tuple(sites[site_index] for site_index in pair)
                break

    for index, denticity in enumerate(denticities):
        if groups[index] is not None:
            continue

        available = [
            site_index for site_index in monodentate_sites
            if site_index not in used
        ]
        chosen = available[:denticity]
        used.update(chosen)
        groups[index] = tuple(sites[site_index] for site_index in chosen)

    return [group for group in groups if group is not None]


def should_use_tridentate_layout(
    parsed: ParsedComplex,
    ligand_items: list[str],
    geometry: str,
) -> bool:
    """Return True when a bis-tridentate octahedral layout will be clearer."""
    if geometry.lower().strip() != "octahedral":
        return False

    if len(ligand_items) != 2:
        return False

    return all(parsed.ligand_denticity.get(ligand, 1) == 3 for ligand in ligand_items)


def should_use_edta_layout(
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


def regular_polygon_sites(n_sites: int, radius: float = 3.2) -> list[Site]:
    """Fallback placement on a regular polygon."""
    if n_sites <= 0:
        return []

    sites: list[Site] = []
    for i in range(n_sites):
        angle = -pi / 2 + 2 * pi * i / n_sites
        sites.append(Site(radius * cos(angle), radius * sin(angle), "plain"))
    return sites


__all__ = [
    "Site",
    "coordination_sites",
    "chelate_octahedral_sites",
    "tridentate_octahedral_sites",
    "edta_octahedral_sites",
    "should_use_chelate_layout",
    "should_use_mixed_polydentate_layout",
    "mixed_polydentate_site_groups",
    "should_use_tridentate_layout",
    "should_use_edta_layout",
    "regular_polygon_sites",
]
