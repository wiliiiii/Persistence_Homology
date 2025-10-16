#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: Compute_AvgBetti_Stats.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  For each chosen homology dimension k, compute across knot types and lengths (N):
    (a) averaged integrals of finite bars in β_k (sum(death-birth) per sample, then average),
    (b) average of per-sample max Betti curve (mean_s max_t β_k^(s)(t)).
  Save static PNGs, and optionally interactive Plotly HTMLs.

Usage:
  python Compute_AvgBetti_Stats.py --root . --base specific_type --dims 2 \
    --Ns 50 100 150 200 --types 0_1 3_1 4_1 5_1 6_1 6_2 6_3 \
    --outdir outputs --interactive --auto-open

Directory layout:
  <root>/<base>/N=<N>/<knot_type>/*.xyz
  (If your data are under "specific_type", pass --base specific_type.)

Requirements:
  - Python 3.9+
  - numpy
  - matplotlib
  - gudhi
  - plotly  (optional; for --interactive)

Notes:
  - AlphaComplex in GUDHI yields intervals in α^2; this script converts to α via sqrt.
  - PNGs are always written; HTML files are written only with --interactive.
"""

from __future__ import annotations

import os, glob, re, math, argparse
import numpy as np
import matplotlib.pyplot as plt
import gudhi as gd

# Optional Plotly (interactive)
try:
    import plotly.graph_objects as go
    from plotly.offline import plot as plot_offline
    _HAS_PLOTLY = True
except Exception:
    _HAS_PLOTLY = False


# ---------- I/O ----------
def _is_float(s: str) -> bool:
    try:
        float(s); return True
    except Exception:
        return False


def load_xyz(path: str) -> np.ndarray:
    """Read last 3 columns as (x,y,z); skip first two .xyz header lines if present."""
    coords = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    start = 2 if len(lines) >= 2 else 0
    for line in lines[start:]:
        parts = line.strip().split()
        if not parts:
            continue
        if len(parts) >= 4 and not _is_float(parts[0]):
            parts = parts[-3:]
        elif len(parts) >= 3:
            parts = parts[:3]
        else:
            continue
        x, y, z = map(float, parts)
        coords.append((x, y, z))
    P = np.asarray(coords, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"Failed to parse coordinates from {path}")
    return P


def densify_polyline_vertices_plus_m(P0: np.ndarray, m: int = 1) -> np.ndarray:
    """Insert (m-1) evenly spaced interior points per edge (m=1: no densify)."""
    if m <= 1:
        return P0.copy()
    P0 = np.asarray(P0, float)
    out = []
    for i in range(len(P0)):
        a = P0[i]
        b = P0[(i + 1) % len(P0)]
        ts = np.arange(0, m) / float(m)  # 0, 1/m, ..., (m-1)/m
        seg = (1 - ts)[:, None] * a[None, :] + ts[:, None] * b[None, :]
        out.append(seg[:-1])  # drop last to avoid duplication
    return np.vstack(out + [P0[-1:]])


# ---------- PH helpers ----------
def alpha_intervals(P: np.ndarray, dim: int) -> np.ndarray:
    """AlphaComplex intervals in α (not α^2)."""
    st = gd.AlphaComplex(points=P).create_simplex_tree()
    st.persistence()
    H_a2 = np.asarray(st.persistence_intervals_in_dimension(dim), dtype=float)
    if H_a2.size == 0:
        return H_a2
    H = H_a2.copy()
    finite = np.isfinite(H)
    H[finite] = np.sqrt(H[finite])  # convert α^2 -> α
    return H


def betti_integral(H: np.ndarray) -> float:
    """Sum of finite bar lengths."""
    if H.size == 0:
        return 0.0
    finite = np.isfinite(H[:, 1])
    return float(np.sum(H[finite, 1] - H[finite, 0]))


def betti_curve_on_grid(H: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """β(t) on a grid: count intervals covering each t (finite bars only)."""
    if H.size == 0:
        return np.zeros_like(grid, dtype=float)
    births = H[:, 0][:, None]
    deaths = H[:, 1][:, None]
    finite = np.isfinite(deaths[:, 0])
    births = births[finite]
    deaths = deaths[finite]
    if births.size == 0:
        return np.zeros_like(grid, dtype=float)
    t = grid[None, :]
    cover = (births <= t) & (t < deaths)
    return cover.sum(axis=0).astype(float)


# ---------- discovery ----------
def discover_Ns(base_root: str):
    """Find integer N from folders named 'N=<N>'."""
    Ns = []
    for name in os.listdir(base_root):
        m = re.fullmatch(r"N=(\d+)", name)
        if m and os.path.isdir(os.path.join(base_root, name)):
            Ns.append(int(m.group(1)))
    Ns.sort()
    return Ns


def collect_types(folder_N: str):
    """List subfolders (knot types) under N=<N>."""
    return sorted([d for d in os.listdir(folder_N) if os.path.isdir(os.path.join(folder_N, d))])


# ---------- core ----------
def compute_stats_for_N_type(
    folder_N: str, knot_type: str, dim: int, grid_steps: int, densify_m: int
):
    """
    Returns:
      avg_integral  : mean over samples of sum(death-birth)
      avg_max_betti : mean over samples of max_t β_k^(s)(t)
    """
    folder = os.path.join(folder_N, knot_type)
    xyz_list = sorted(glob.glob(os.path.join(folder, "*.xyz")))
    if not xyz_list:
        return None, None

    H_list, integrals = [], []
    deaths_all = []
    for p in xyz_list:
        try:
            P = load_xyz(p)
            P = densify_polyline_vertices_plus_m(P, densify_m)
            H = alpha_intervals(P, dim)
            H_list.append(H)
            integrals.append(betti_integral(H))
            if H.size:
                finite = np.isfinite(H[:, 1])
                deaths_all.extend(H[finite, 1].tolist())
        except Exception as e:
            print(f"[WARN] {p}: {e}")

    if not H_list:
        return None, None

    t_max = max(deaths_all) if deaths_all else 1.0
    if not math.isfinite(t_max) or t_max <= 0:
        t_max = 1.0
    grid = np.linspace(0.0, t_max, num=grid_steps, endpoint=False)

    avg_integral = float(np.mean(np.array(integrals))) if integrals else 0.0

    per_sample_maxima = []
    for H in H_list:
        curve = betti_curve_on_grid(H, grid)
        per_sample_maxima.append(float(np.max(curve)) if curve.size else 0.0)
    avg_max_betti = float(np.mean(np.array(per_sample_maxima))) if per_sample_maxima else 0.0

    return avg_integral, avg_max_betti


# ---------- interactive figures ----------
def plotly_lines(Ns, series_dict, title, y_label, out_html, auto_open=False):
    if not _HAS_PLOTLY:
        print("[WARN] Plotly not installed; skip interactive HTML. Try: pip install plotly")
        return
    fig = go.Figure()
    for tp, ys in series_dict.items():
        fig.add_trace(go.Scatter(
            x=Ns, y=ys, mode="lines+markers", name=tp,
            hovertemplate="N=%{x}<br>%{y:.4g}<extra>"+tp+"</extra>"
        ))
    fig.update_layout(
        title=title,
        xaxis_title="knot length (N)",
        yaxis_title=y_label,
        legend_title="type",
        margin=dict(l=10, r=10, t=50, b=10),
    )
    plot_offline(fig, filename=out_html, auto_open=auto_open, include_plotlyjs="cdn")
    print("[Saved]", out_html)


# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser(description="Compute averaged Betti integrals and average-of-max-Betti vs. N.")
    ap.add_argument("--root", type=str, default=".")
    ap.add_argument("--base", type=str, default="data",
                    help="subfolder under <root> containing N=<N>/ (e.g., 'data' or 'specific_type')")
    ap.add_argument("--dims", type=int, nargs="+", default=[1], help="homology dimensions, e.g., 1 2")
    ap.add_argument("--Ns", type=int, nargs="*", default=None,
                    help="explicit list of N; if omitted, auto-discover under <root>/<base>")
    ap.add_argument("--types", type=str, nargs="*", default=None,
                    help="knot types to include (e.g., 0_1 3_1 5_1); if omitted, use all available")
    ap.add_argument("--grid-steps", type=int, default=1200, help="grid size to sample Betti curves")
    ap.add_argument("--densify-m", type=int, default=1, help="per-edge densification factor (m>=1)")
    ap.add_argument("--outdir", type=str, default="outputs")
    ap.add_argument("--interactive", action="store_true", help="also save Plotly HTML figures")
    ap.add_argument("--auto-open", action="store_true", help="open HTML in browser after saving")
    args = ap.parse_args()

    base_root = os.path.join(args.root, args.base)
    if not os.path.isdir(base_root):
        raise SystemExit(f"Base folder not found: {base_root}")

    Ns = args.Ns if args.Ns else discover_Ns(base_root)
    if not Ns:
        raise SystemExit(f"No N=... folders found under {base_root}")
    Ns = sorted(Ns)

    os.makedirs(args.outdir, exist_ok=True)

    for dim in args.dims:
        # union of types across Ns
        type_set = set()
        for N in Ns:
            folder_N = os.path.join(base_root, f"N={N}")
            if os.path.isdir(folder_N):
                type_set.update(collect_types(folder_N))

        if args.types:
            types = [tp for tp in args.types if tp in type_set]
            missing = [tp for tp in args.types if tp not in type_set]
            if missing:
                print(f"[WARN] requested types not found (ignored): {missing}")
        else:
            types = sorted(type_set)

        if not types:
            print(f"[WARN] No matching knot types for dim={dim}")
            continue

        y_integral = {tp: [] for tp in types}
        y_avgmax  = {tp: [] for tp in types}

        for N in Ns:
            folder_N = os.path.join(base_root, f"N={N}")
            for tp in types:
                if os.path.isdir(os.path.join(folder_N, tp)):
                    ai, am = compute_stats_for_N_type(folder_N, tp, dim, args.grid_steps, args.densify_m)
                else:
                    ai, am = None, None
                if ai is None:
                    y_integral[tp].append(np.nan)
                    y_avgmax[tp].append(np.nan)
                else:
                    y_integral[tp].append(ai)
                    y_avgmax[tp].append(am)

        # ---- Matplotlib PNGs ----
        plt.figure()
        for tp in types:
            plt.plot(Ns, y_integral[tp], marker='o', label=tp)
        plt.xlabel("knot length (N)")
        plt.ylabel(f"averaged integrals (dim={dim})")
        plt.title(f"Averaged integrals vs. knot length (β_{dim})")
        plt.legend()
        plt.tight_layout()
        f1 = os.path.join(args.outdir, f"avg_integral_vs_length_dim{dim}.png")
        plt.savefig(f1, dpi=220)
        print("[Saved]", f1)

        plt.figure()
        for tp in types:
            plt.plot(Ns, y_avgmax[tp], marker='o', label=tp)
        plt.xlabel("knot length (N)")
        plt.ylabel(f"average of maxima Betti (dim={dim})")
        plt.title(f"Average of maxima Betti vs. knot length (β_{dim})")
        plt.legend()
        plt.tight_layout()
        f2 = os.path.join(args.outdir, f"avg_of_max_betti_vs_length_dim{dim}.png")
        plt.savefig(f2, dpi=220)
        print("[Saved]", f2)

        # ---- Plotly HTML (optional) ----
        if args.interactive:
            html1 = os.path.join(args.outdir, f"avg_integral_vs_length_dim{dim}.html")
            html2 = os.path.join(args.outdir, f"avg_of_max_betti_vs_length_dim{dim}.html")
            plotly_lines(Ns, y_integral,
                         title=f"Averaged integrals vs. knot length (β_{dim})",
                         y_label=f"averaged integrals (dim={dim})",
                         out_html=html1, auto_open=args.auto_open)
            plotly_lines(Ns, y_avgmax,
                         title=f"Average of maxima Betti vs. knot length (β_{dim})",
                         y_label=f"average of maxima Betti (dim={dim})",
                         out_html=html2, auto_open=args.auto_open)


if __name__ == "__main__":
    main()
