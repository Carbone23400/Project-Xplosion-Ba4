"""
coordchem/viz/transform_2d.py
-----------------------------
Geometric helpers and ligand-placement transforms for 2D coordination diagrams.
"""

from __future__ import annotations

from math import atan2, cos, sin, sqrt

from .layout_2d import Site


def dist(x: float, y: float) -> float:
    return sqrt(x * x + y * y)


def angle(x: float, y: float) -> float:
    return atan2(y, x)


def max_distance_from_point(
    coords: dict[int, tuple[float, float]],
    center: tuple[float, float],
) -> float:
    """Return the maximum distance of any coordinate from center."""
    if not coords:
        return 1.0
    return max(dist(x - center[0], y - center[1]) for x, y in coords.values())


def rotate_point(x: float, y: float, theta: float) -> tuple[float, float]:
    c = cos(theta)
    s = sin(theta)
    return (c * x - s * y, s * x + c * y)


def reflect_point_across_line(
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


def centroid(
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


def transform_monodentate(
    coords: dict[int, tuple[float, float]],
    donor_idx: int,
    anchor: Site,
) -> dict[int, tuple[float, float]]:
    """Place a monodentate ligand with donor atom on the coordination site."""
    donor_x, donor_y = coords[donor_idx]
    cx, cy = centroid(coords, exclude={donor_idx})

    current_vx = cx - donor_x
    current_vy = cy - donor_y
    if dist(current_vx, current_vy) < 1e-8:
        current_vx, current_vy = 1.0, 0.0

    theta = angle(anchor.x, anchor.y) - angle(current_vx, current_vy)

    new_coords: dict[int, tuple[float, float]] = {}
    for idx, (x, y) in coords.items():
        rx, ry = rotate_point(x - donor_x, y - donor_y, theta)
        new_coords[idx] = (rx + anchor.x, ry + anchor.y)

    return new_coords


def transform_polydentate(
    coords: dict[int, tuple[float, float]],
    donor_indices: tuple[int, ...],
    anchors: tuple[Site, ...],
    *,
    shrink_large_ligand: bool = True,
) -> dict[int, tuple[float, float]]:
    """Place a polydentate ligand on multiple coordination sites."""
    if len(donor_indices) != len(anchors):
        raise ValueError("Number of donor indices and anchors must match.")

    if len(donor_indices) == 1:
        return transform_monodentate(coords, donor_indices[0], anchors[0])

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
    cur_len = dist(cur_vx, cur_vy) or 1.0

    tar_mid_x = (a1.x + a2.x) / 2
    tar_mid_y = (a1.y + a2.y) / 2
    tar_vx = a2.x - a1.x
    tar_vy = a2.y - a1.y
    tar_len = dist(tar_vx, tar_vy) or 1.0

    theta = angle(tar_vx, tar_vy) - angle(cur_vx, cur_vy)
    scale = tar_len / cur_len

    new_coords: dict[int, tuple[float, float]] = {}
    for idx, (x, y) in coords.items():
        x0 = (x - cur_mid_x) * scale
        y0 = (y - cur_mid_y) * scale
        xr, yr = rotate_point(x0, y0, theta)
        new_coords[idx] = (xr + tar_mid_x, yr + tar_mid_y)

    ligand_radius = max_distance_from_point(new_coords, (tar_mid_x, tar_mid_y))
    target_radius = 3.2 if len(donor_indices) >= 3 else 2.8
    if shrink_large_ligand and len(donor_indices) < 3 and ligand_radius > target_radius:
        shrink = target_radius / ligand_radius
        new_coords = {
            idx: (
                tar_mid_x + (x - tar_mid_x) * shrink,
                tar_mid_y + (y - tar_mid_y) * shrink,
            )
            for idx, (x, y) in new_coords.items()
        }

    lig_centroid = centroid(new_coords, exclude=set(donor_indices))
    outward_x = tar_mid_x
    outward_y = tar_mid_y
    ligand_side_x = lig_centroid[0] - tar_mid_x
    ligand_side_y = lig_centroid[1] - tar_mid_y

    if ligand_side_x * outward_x + ligand_side_y * outward_y < 0:
        line_a = (a1.x, a1.y)
        line_b = (a2.x, a2.y)
        new_coords = {
            idx: reflect_point_across_line(point, line_a, line_b)
            for idx, point in new_coords.items()
        }

    return new_coords


def transform_acac(
    coords: dict[int, tuple[float, float]],
    donor_indices: tuple[int, ...],
    anchors: tuple[Site, ...],
) -> dict[int, tuple[float, float]]:
    """Place acac while keeping donor sites exact and terminal methyls separated."""
    new_coords = transform_polydentate(
        coords,
        donor_indices,
        anchors,
        shrink_large_ligand=False,
    )

    if len(donor_indices) != 2:
        return new_coords

    a1, a2 = anchors
    has_top_plain = any(anchor.style == "plain" and anchor.y > 0 for anchor in anchors)
    has_bottom_plain = any(anchor.style == "plain" and anchor.y < 0 for anchor in anchors)
    is_right_side = a1.x > 0 and a2.x > 0

    adjusted = dict(new_coords)

    if has_top_plain and 6 in adjusted:
        x, y = adjusted[6]
        adjusted[6] = (x + 0.10, y + 0.45)

    if has_bottom_plain and 0 in adjusted:
        x, y = adjusted[0]
        adjusted[0] = (x + 0.10, y - 0.45)

    if is_right_side:
        if 0 in adjusted:
            x, y = adjusted[0]
            adjusted[0] = (x + 0.30, y - 0.45)
        if 6 in adjusted:
            x, y = adjusted[6]
            adjusted[6] = (x + 0.30, y + 0.45)

    return adjusted


def transform_oxalate(
    coords: dict[int, tuple[float, float]],
    donor_indices: tuple[int, ...],
    anchors: tuple[Site, ...],
) -> dict[int, tuple[float, float]]:
    """Place oxalate with both carbonyl bonds pointing away from the metal."""
    if len(donor_indices) != 2 or len(anchors) != 2:
        return transform_polydentate(coords, donor_indices, anchors)

    a1, a2 = anchors
    mid_x = (a1.x + a2.x) / 2
    mid_y = (a1.y + a2.y) / 2
    outward_len = dist(mid_x, mid_y)
    if outward_len < 1e-8:
        outward_x, outward_y = 0.0, 1.0
    else:
        outward_x = mid_x / outward_len
        outward_y = mid_y / outward_len

    chord_x = a2.x - a1.x
    chord_y = a2.y - a1.y
    chord_len = dist(chord_x, chord_y) or 1.0
    along_x = chord_x / chord_len
    along_y = chord_y / chord_len

    carbon_out = 1.34
    carbon_inset = max(0.0, (chord_len - 1.85) / 2)
    carbonyl_out = 0.94

    c1 = (
        a1.x + outward_x * carbon_out + along_x * carbon_inset,
        a1.y + outward_y * carbon_out + along_y * carbon_inset,
    )
    c2 = (
        a2.x + outward_x * carbon_out - along_x * carbon_inset,
        a2.y + outward_y * carbon_out - along_y * carbon_inset,
    )
    a1_out_len = dist(a1.x, a1.y) or 1.0
    a2_out_len = dist(a2.x, a2.y) or 1.0
    a1_out_x = a1.x / a1_out_len
    a1_out_y = a1.y / a1_out_len
    a2_out_x = a2.x / a2_out_len
    a2_out_y = a2.y / a2_out_len

    o1 = (c1[0] + a1_out_x * carbonyl_out, c1[1] + a1_out_y * carbonyl_out)
    o2 = (c2[0] + a2_out_x * carbonyl_out, c2[1] + a2_out_y * carbonyl_out)

    d1, d2 = donor_indices
    adjusted = dict(coords)
    adjusted[d1] = (a1.x, a1.y)
    adjusted[1] = c1
    adjusted[2] = o1
    adjusted[3] = c2
    adjusted[4] = o2
    adjusted[d2] = (a2.x, a2.y)
    return adjusted


def transform_edta(anchors: tuple[Site, ...]) -> dict[int, tuple[float, float]]:
    """Place EDTA in a wrapped hexadentate drawing around the metal."""
    if len(anchors) != 6:
        raise ValueError("EDTA layout requires six coordination sites.")

    n_top_left, o_top, o_top_right, n_bottom_left, o_bottom, o_bottom_right = anchors

    return {
        0: (n_top_left.x, n_top_left.y),
        1: (-1.45, 1.90),
        2: (-0.55, 2.35),
        3: (-0.65, 3.05),
        4: (o_top.x, o_top.y),
        5: (-0.95, 1.65),
        6: (0.55, 1.65),
        7: (0.95, 2.28),
        8: (o_top_right.x, o_top_right.y),
        9: (-2.45, 0.50),
        10: (-2.45, -0.50),
        11: (n_bottom_left.x, n_bottom_left.y),
        12: (-1.45, -1.90),
        13: (-0.55, -2.35),
        14: (-0.65, -3.05),
        15: (o_bottom.x, o_bottom.y),
        16: (-0.95, -1.65),
        17: (0.55, -1.65),
        18: (0.95, -2.28),
        19: (o_bottom_right.x, o_bottom_right.y),
    }
