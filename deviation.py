#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: deviation.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Batch-compute the Alpha-based deviation metric δ_ε(K) from H1 barcodes for
  each knot type under a fixed N, save an Excel per type, and export β1(t) curves.

Usage (Windows CMD with line breaks):
  python deviation.py ^
    --root "~\\SpecificPLKnots" ^
    --N 50 ^
    --eps 0.02 ^
    --min-pers 0.05

Directory layout (fixed):
  <root>/specific_type/N=<N>/<type>/*.xyz
Outputs:
  Excel:  <root>/deviation/N=<N>/<type>/deviation.xlsx
  Plots:  <root>/deviation/N=<N>/<type>/beta1_curves/*.png

Requirements:
  - Python 3.9+
  - numpy, pandas, matplotlib
  - gudhi
  - openpyxl (pandas Excel writer)
"""

from __future__ import annotations

import os, glob, argparse
import numpy as np
import pandas as pd

# headless (set backend before importing pyplot)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gudhi as gd


# ---------- δ_ε helper ----------
def _integral_linear_cutoff(t1: float, t2: float, R: float, eps: float) -> float:
    """Integral of (1 - t/(R-eps)) on [t1, t2] assuming 0 <= t1 < t2 <= R-eps."""
    return (t2 - t1) - (t2**2 - t1**2) / (2.0 * (R - eps))


def delta_epsilon_from_barcodes(h1_intervals: np.ndarray, eps: float) -> float:
    """
    δ_ε(K) for a set of finite H1 bars in Alpha radius (NOT α^2).
    Implements the linear cutoff at (R - eps), where R = max death.
    Returns δ_ε(K) / R (normalized by R) as in your formula.
    """
    bars = np.array([bd for bd in h1_intervals if np.isfinite(bd[1])], float)
    if bars.size == 0:
        return 0.0
    R = float(np.max(bars[:, 1]))
    if not (0.0 < eps < R):
        return 0.0
    cutoff = R - eps

    # event points for sweep (including boundaries and all bar endpoints)
    pts = np.unique(np.clip(np.concatenate(([0.0, cutoff, R], bars.reshape(-1))), 0.0, R))
    if not np.any(np.isclose(pts, cutoff)):
        pts = np.sort(np.append(pts, cutoff))

    events = []
    for b, d in bars:
        b = max(0.0, min(b, R)); d = max(0.0, min(d, R))
        events.append((b, +1)); events.append((d, -1))
    events.sort(key=lambda x: (x[0], -x[1]))

    beta = 0; sweep = 0; total = 0.0
    for i in range(len(pts) - 1):
        t1, t2 = pts[i], pts[i + 1]
        if t1 >= R or t1 >= cutoff:
            break
        while sweep < len(events) and np.isclose(events[sweep][0], t1):
            beta += events[sweep][1]; sweep += 1
        while sweep < len(events) and (events[sweep][0] < t1):
            beta += events[sweep][1]; sweep += 1
        k = max(beta - 1, 0)  # only count when β1 >= 2
        if k > 0:
            a, bnd = t1, min(t2, cutoff)
            if bnd > a:
                total += k * _integral_linear_cutoff(a, bnd, R, eps)
    return total / R


# ---------- Alpha complex ----------
def load_points(path: str) -> np.ndarray:
    """Robust .xyz-like loader: skip header, allow (label x y z) or (x y z)."""
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    def _has_xyz_header(ls):
        try:
            n = int(ls[0]); return n > 0 and len(ls) >= n + 1
        except Exception:
            return False

    start = 2 if _has_xyz_header(lines) else 0
    for ln in lines[start:]:
        ln = ln.replace(",", " ")
        toks = [t for t in ln.split() if t]
        if toks and any(c.isalpha() for c in toks[0]):   # element label present
            toks = toks[1:]
        nums = []
        for t in toks:
            try:
                nums.append(float(t))
            except Exception:
                pass
        if len(nums) >= 3:
            rows.append(nums[:3])

    pts = np.array(rows, float)
    if pts.ndim != 2 or pts.shape[0] < 4 or pts.shape[1] != 3:
        raise ValueError(f"Bad file {path}: got shape {pts.shape}")
    return pts


def alpha_h1_intervals(points: np.ndarray, take_sqrt: bool = True) -> np.ndarray:
    """Return H1 intervals in Alpha radius (sqrt of α^2 if take_sqrt=True)."""
    ac = gd.AlphaComplex(points=points)
    st = ac.create_simplex_tree()
    st.persistence(homology_coeff_field=2, min_persistence=0.0)
    h1 = st.persistence_intervals_in_dimension(1) or []
    h1 = np.asarray(h1, float)
    if h1.size == 0:
        return np.empty((0, 2), float)
    if take_sqrt:
        h1 = np.where(np.isfinite(h1), np.sqrt(h1), h1)  # convert α^2 -> α
    return h1


def filter_by_persistence(intervals: np.ndarray, min_pers: float) -> np.ndarray:
    """Keep finite bars with length >= min_pers."""
    if intervals.size == 0:
        return intervals
    pers = intervals[:, 1] - intervals[:, 0]
    keep = np.isfinite(intervals[:, 1]) & (pers >= float(min_pers))
    return intervals[keep]


# ---------- β1(t) figure ----------
def build_beta1_steps(intervals: np.ndarray):
    """Build step function samples (x,y) for β1(t) with finite bars."""
    bars = np.array([bd for bd in intervals if np.isfinite(bd[1])], float)
    if bars.size == 0:
        return np.array([0.0, 1.0]), np.array([0, 0])
    times = np.unique(np.sort(bars.reshape(-1)))
    if times.size < 2:
        return np.array([times[0], times[0]]), np.array([0, 0])

    from collections import Counter
    births = Counter(bars[:, 0]); deaths = Counter(bars[:, 1])
    beta = 0; bvals = []
    for i in range(times.size - 1):
        t = times[i]
        beta += births.get(t, 0); bvals.append(beta); beta -= deaths.get(t, 0)

    x, y = [], []
    for i in range(times.size - 1):
        x.extend([times[i], times[i + 1]]); y.extend([bvals[i], bvals[i]])
    return np.array(x, float), np.array(y, float)


def save_beta1_curve(intervals: np.ndarray, out_png: str, title: str) -> None:
    """Save a step plot of β1(t)."""
    x, y = build_beta1_steps(intervals)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.figure()
    plt.step(x, y, where="post")
    plt.xlabel("scale t (alpha radius)")
    plt.ylabel(r"$\beta_1(t)$")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close()


# ---------- batch ----------
def main() -> None:
    ap = argparse.ArgumentParser(description="Batch Alpha δ_ε(K) -> Excel + β1 curves")
    ap.add_argument("--root", required=True, help="Folder that contains 'specific_type'")
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--types", nargs="*", default=None, help="Subset of types; default = all under N=<N>")
    ap.add_argument("--eps", type=float, required=True, help="Epsilon for δ_ε cutoff")
    ap.add_argument("--min-pers", type=float, default=0.0, help="Filter bars with persistence < min-pers")
    ap.add_argument("--no-sqrt", action="store_true", help="Do NOT sqrt α^2 -> α")
    args = ap.parse_args()

    # Expand user (~) and normalize
    root = os.path.abspath(os.path.expandvars(os.path.expanduser(args.root)))
    spec_dir = os.path.join(root, "specific_type", f"N={args.N}")
    out_root = os.path.join(root, "deviation", f"N={args.N}")
    os.makedirs(out_root, exist_ok=True)

    # Collect type folders
    if args.types:
        type_dirs = [os.path.join(spec_dir, t) for t in args.types]
    else:
        type_dirs = [
            os.path.join(spec_dir, d)
            for d in os.listdir(spec_dir)
            if os.path.isdir(os.path.join(spec_dir, d))
        ]

    take_sqrt = not args.no_sqrt

    for tdir in sorted(type_dirs):
        tname = os.path.basename(tdir)
        type_out = os.path.join(out_root, tname)
        curves_dir = os.path.join(type_out, "beta1_curves")
        os.makedirs(curves_dir, exist_ok=True)

        rows = []  # columns: label, delta_epsilon
        files = sorted(glob.glob(os.path.join(tdir, "*.xyz")))
        for fp in files:
            label = os.path.splitext(os.path.basename(fp))[0]
            try:
                pts = load_points(fp)
                h1 = alpha_h1_intervals(pts, take_sqrt=take_sqrt)
                h1 = filter_by_persistence(h1, args.min_pers)
                delta = 0.0 if h1.size == 0 else delta_epsilon_from_barcodes(h1, args.eps)

                # figure
                png_path = os.path.join(curves_dir, label + ".png")
                save_beta1_curve(h1, png_path, f"Beta1 (Alpha, min-pers={args.min_pers}) - {label}")

                rows.append({"label": label, "delta_epsilon": float(delta)})
                print(f"[OK] {tname}/{label}  delta={delta:.6g}")
            except Exception as e:
                # keep NaN row to preserve counting
                rows.append({"label": label, "delta_epsilon": np.nan})
                print(f"[ERR] {tname}/{label}: {e}")

        # write Excel per type
        df = pd.DataFrame(rows, columns=["label", "delta_epsilon"])
        excel_path = os.path.join(type_out, "deviation.xlsx")
        os.makedirs(type_out, exist_ok=True)
        df.to_excel(excel_path, index=False)
        print(f"[Saved] {excel_path}")
        print(f"[Saved] curves -> {curves_dir}")

    print("[Done]")


if __name__ == "__main__":
    main()
