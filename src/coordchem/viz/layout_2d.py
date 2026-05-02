"""
coordchem/viz/layout_2d.py
--------------------------
2D coordination-site layouts for coordination complex depictions.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, sin

from ..parser import ParsedComplex


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
    "should_use_tridentate_layout",
    "should_use_edta_layout",
    "regular_polygon_sites",
]
