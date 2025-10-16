#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: ComputeQuantities.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  One-pass batch computation over .xyz polygonal loops:
    - Alpha(H1) -> I1 (finite-bar integral), H1_max (max persistence)
    - Alpha(H2) -> I2, H2_max
    - Geometry: curvature, torsion, diameter, circumscribing-sphere volume,
                convex-hull volume (optional), radius of gyration
  Output: ONE Excel sheet (or CSV fallback) with columns:
    [I1, H1_max, I2, H2_max, curvature, torsion, diameter,
     Volume_of_sphere, Volume_of_convex_hull, radius_of_gyration]

Usage:
  # Example (Windows path must be quoted and keep backslashes):
  python ComputeQuantities.py --input-dir "~\\output_coordinate"
  # Optional:
  #   --glob "*.xyz"
  #   --pts-per-edge 10
  #   --excel "C:\\path\\to\\results.xlsx"
  #   --sheet results

Requirements:
  - Python 3.9+
  - numpy
  - gudhi
  - pandas (optional; for Excel writing)
  - openpyxl (optional; pandas Excel engine)
  - scipy (optional; for ConvexHull volume)

Notes:
  - If pandas/openpyxl are missing, a CSV fallback <results.csv> is written.
  - --pts-per-edge densifies each segment by linear interpolation before AlphaComplex.
"""

from __future__ import annotations

import os
import glob
import math
import argparse
from typing import List
import numpy as np

# pandas is optional; fallback to CSV if unavailable
try:
    import pandas as pd  # type: ignore
except Exception:
    pd = None  # type: ignore

import gudhi as gd

# Optional ConvexHull support
try:
    from scipy.spatial import ConvexHull  # type: ignore
    _has_convex_hull = True
    _convex_import_err = None
except Exception as e:
    _has_convex_hull = False
    _convex_import_err = e


# ---------- I/O ----------
def load_xyz(path: str) -> np.ndarray:
    """Parse a simple XYZ-like file (with a count header)."""
    coords: List[tuple[float, float, float]] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    start = 2 if len(lines) >= 2 else 0  # skip header lines if present
    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) >= 4:
            x, y, z = map(float, parts[-3:])
            coords.append((x, y, z))
    P = np.asarray(coords, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"Failed to parse coordinates from {path}")
    return P


def densify_polyline_vertices_plus_m(P0: np.ndarray, m: int = 10) -> np.ndarray:
    """Linearly densify each edge into m segments (excluding the last endpoint)."""
    P0 = np.asarray(P0, float)
    if m <= 1:
        return P0.copy()
    out = []
    for i in range(len(P0)):
        a = P0[i]
        b = P0[(i + 1) % len(P0)]
        ts = np.arange(0, m) / float(m)
        seg = (1 - ts)[:, None] * a[None, :] + ts[:, None] * b[None, :]
        out.append(seg)
    return np.vstack(out)


# ---------- PH helpers ----------
def betti_integral(H: np.ndarray) -> float:
    """Sum of finite bar lengths."""
    if H.size == 0:
        return 0.0
    finite = np.isfinite(H[:, 1])
    return float(np.sum(H[finite, 1] - H[finite, 0]))


def max_persistence(H: np.ndarray) -> float:
    """Max finite bar length."""
    if H.size == 0:
        return 0.0
    finite = np.isfinite(H[:, 1])
    if not np.any(finite):
        return 0.0
    return float(np.max(H[finite, 1] - H[finite, 0]))


# ---------- Geometric helpers ----------
def _safe_acos(x: float) -> float:
    return math.acos(max(-1.0, min(1.0, float(x))))


def compute_curvature(points_list: np.ndarray) -> float:
    """Discrete total curvature normalized by 2π."""
    P = np.asarray(points_list, float)
    N = len(P)
    total = 0.0
    for i in range(N):
        a = P[i] - P[i - 1]
        b = P[(i + 1) % N] - P[i]
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        if na == 0.0 or nb == 0.0:
            continue
        cosang = np.dot(a, b) / (na * nb)
        ang = _safe_acos(cosang)
        total += ang
    return float(total / (2 * math.pi))


def compute_torsion(points_list: np.ndarray) -> float:
    """Discrete total torsion (angle between consecutive binormals), normalized by 2π."""
    P = np.asarray(points_list, float)
    N = len(P)
    total = 0.0
    for i in range(N):
        a = P[i] - P[i - 1]
        b = P[(i + 1) % N] - P[i]
        c = P[(i + 2) % N] - P[(i + 1) % N]
        na, nb, nc = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
        if min(na, nb, nc) == 0.0:
            continue
        b1 = np.cross(a, b)
        b2 = np.cross(b, c)
        nb1, nb2 = np.linalg.norm(b1), np.linalg.norm(b2)
        if nb1 == 0.0 or nb2 == 0.0:
            ang = 0.0
        else:
            b1 /= nb1
            b2 /= nb2
            ang = _safe_acos(np.clip(np.dot(b1, b2), -1.0, 1.0))
        total += ang
    return float(total / (2 * math.pi))


def compute_max_distance(P: np.ndarray) -> float:
    """Brute-force max pairwise distance (O(n^2/2))."""
    P = np.asarray(P, float)
    if len(P) == 0:
        return 0.0
    dmax = 0.0
    for i in range(len(P)):
        if i + 1 < len(P):
            D = np.linalg.norm(P[i + 1 :] - P[i], axis=1)
            if D.size:
                dmax = max(dmax, float(np.max(D)))
    return dmax


def circumscribing_sphere_volume_from_diameter(d: float) -> float:
    """Volume of sphere with diameter d."""
    r = 0.5 * float(d)
    return (4.0 / 3.0) * math.pi * (r**3)


def convex_hull_volume(points_list: np.ndarray) -> float:
    """3D convex hull volume (requires scipy)."""
    if not _has_convex_hull:
        raise RuntimeError(f"scipy not available for ConvexHull: {_convex_import_err}")
    P = np.asarray(points_list, float)
    if len(P) < 4:
        return 0.0
    hull = ConvexHull(P)  # type: ignore
    return float(hull.volume)


def radius_of_gyration(points_list: np.ndarray) -> float:
    """Root mean square distance to centroid."""
    P = np.asarray(points_list, float)
    if len(P) == 0:
        return 0.0
    c = np.mean(P, axis=0)
    sq = np.sum((P - c) ** 2, axis=1)
    return float(math.sqrt(np.mean(sq)))


# ---------- main ----------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="All-in-one: I1/H1_max, I2/H2_max, and geometric metrics -> ONE Excel"
    )
    ap.add_argument(
        "--input-dir",
        type=str,
        default=r"C:\Users\LENOVO\Desktop\UoR\Paper_Topoly\output_coordinate",
        help="Directory containing .xyz files (Windows path: quote and keep backslashes)",
    )
    ap.add_argument("--glob", type=str, default="*.xyz", help="Filename glob pattern")
    ap.add_argument("--pts-per-edge", type=int, default=10, help="Densification per edge (m >= 1)")
    ap.add_argument(
        "--excel",
        type=str,
        default=None,
        help="Output Excel path (.xlsx). Default: <input-dir>\\results.xlsx",
    )
    ap.add_argument("--sheet", type=str, default="results", help="Excel sheet name")
    args = ap.parse_args()

    in_dir = args.input_dir
    if not os.path.isdir(in_dir):
        raise SystemExit(f"Input dir not found: {in_dir}")
    files = sorted(glob.glob(os.path.join(in_dir, args.glob)))
    if not files:
        raise SystemExit(f"No .xyz matched: {os.path.join(in_dir, args.glob)}")

    out_xlsx = args.excel or os.path.join(in_dir, "results.xlsx")

    names: List[str] = []
    I1_vals: List[float] = []
    H1max_vals: List[float] = []
    I2_vals: List[float] = []
    H2max_vals: List[float] = []
    curvatures: List[float] = []
    torsions: List[float] = []
    diameters: List[float] = []
    V_spheres: List[float] = []
    V_convexes: List[float] = []
    r_gyrs: List[float] = []

    m = max(1, int(args.pts_per_edge))

    for i, path in enumerate(files, 1):
        name = os.path.basename(path)
        P0 = load_xyz(path)
        PK = densify_polyline_vertices_plus_m(P0, m=m)  # densified for AlphaComplex

        # Alpha complex on densified points
        ac = gd.AlphaComplex(points=PK)
        st = ac.create_simplex_tree()
        st.persistence(homology_coeff_field=2)

        def alpha2_to_alpha(Ha2: np.ndarray) -> np.ndarray:
            """Convert Alpha^2 birth/death to Alpha by sqrt (preserve inf)."""
            if Ha2.size == 0:
                return np.empty((0, 2))
            H = Ha2.astype(float, copy=True)
            finite = np.isfinite(H)
            H[finite] = np.sqrt(H[finite])  # α = sqrt(α^2)
            return H

        H1 = alpha2_to_alpha(np.array(st.persistence_intervals_in_dimension(1)))
        H2 = alpha2_to_alpha(np.array(st.persistence_intervals_in_dimension(2)))

        I1 = betti_integral(H1)
        H1max = max_persistence(H1)
        I2 = betti_integral(H2)
        H2max = max_persistence(H2)

        dmax = compute_max_distance(PK)
        try:
            V_conv = convex_hull_volume(PK)
        except RuntimeError as e:
            print(f"[WARN] {name}: {e}; set Volume_of_convex_hull=0.0")
            V_conv = 0.0

        names.append(name)
        I1_vals.append(I1)
        H1max_vals.append(H1max)
        I2_vals.append(I2)
        H2max_vals.append(H2max)
        curvatures.append(compute_curvature(P0))
        torsions.append(compute_torsion(P0))
        diameters.append(dmax)
        V_spheres.append(circumscribing_sphere_volume_from_diameter(dmax))
        V_convexes.append(V_conv)
        r_gyrs.append(radius_of_gyration(PK))

        print(f"[{i:03d}/{len(files)}] {name}  I1={I1:.6f}  I2={I2:.6f}  dmax={dmax:.6f}")

    # Fallback to CSV if pandas is missing
    if pd is None:  # type: ignore
        out_csv = os.path.splitext(out_xlsx)[0] + ".csv"
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            f.write("name,I1,H1_max,I2,H2_max,curvature,torsion,diameter,Volume_of_sphere,Volume_of_convex_hull,radius_of_gyration\n")
            for i, nm in enumerate(names):
                f.write(
                    f"{nm},{I1_vals[i]:.10f},{H1max_vals[i]:.10f},{I2_vals[i]:.10f},{H2max_vals[i]:.10f},"
                    f"{curvatures[i]:.10f},{torsions[i]:.10f},{diameters[i]:.10f},"
                    f"{V_spheres[i]:.10f},{V_convexes[i]:.10f},{r_gyrs[i]:.10f}\n"
                )
        print(f"[Done] pandas not found; CSV written -> {out_csv}")
        return

    # Write Excel with pandas
    try:
        import openpyxl  # type: ignore  # ensure engine available
    except Exception:
        pass  # if missing, pandas may still handle via default writer

    df = pd.DataFrame(  # type: ignore
        {
            "I1": I1_vals,
            "H1_max": H1max_vals,
            "I2": I2_vals,
            "H2_max": H2max_vals,
            "curvature": curvatures,
            "torsion": torsions,
            "diameter": diameters,
            "Volume_of_sphere": V_spheres,
            "Volume_of_convex_hull": V_convexes,
            "radius_of_gyration": r_gyrs,
        },
        index=names,
    )

    os.makedirs(os.path.dirname(out_xlsx), exist_ok=True)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl", mode="w") as writer:  # type: ignore
        df.to_excel(writer, sheet_name=args.sheet, index=True)

    print(f"[Done] Excel written -> {out_xlsx}  (sheet='{args.sheet}')")


if __name__ == "__main__":
    main()
