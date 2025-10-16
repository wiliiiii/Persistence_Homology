#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: GenRandomKnot.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Batch-generate random polygonal loops (knots) with TopoLy and save both
  coordinates (.xyz) and preview images (.png).

Usage:
  python GenRandomKnot.py --length 5 --count 500
  # Common flags:
  #   --bond-length 1.0          # edge length used by TopoLy
  #   --coord-dir output_coordinate
  #   --png-dir output_png
  #   --basename topoly_random
  #   --show                     # open interactive windows (matplotlib)
  #   --save-png                 # when --show, also save PNG

Requirements:
  - Python 3.9+
  - numpy
  - matplotlib
  - topoly  (https://pypi.org/project/topoly/)

Notes:
  - If --show is not set, a non-interactive backend is used and PNGs are saved automatically.
  - .xyz format here is a simple XYZ-like dump with a header line count.
"""

import os
import argparse
import numpy as np

# Parse minimal pre-args to choose matplotlib backend before importing pyplot
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("--show", action="store_true", help="Open interactive window(s)")
_pre.add_argument("--save-png", action="store_true", help="Also save PNG when --show is used")
_pre_args, _ = _pre.parse_known_args()

import matplotlib
if not _pre_args.show:
    matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt


def set_axes_equal(ax, P: np.ndarray) -> None:
    """Make 3D axes have equal scale for a nicer aspect."""
    x, y, z = P[:, 0], P[:, 1], P[:, 2]
    cx, cy, cz = (x.min() + x.max()) / 2, (y.min() + y.max()) / 2, (z.min() + z.max()) / 2
    r = 0.5 * max(x.max() - x.min(), y.max() - y.min(), z.max() - z.min())
    ax.set_xlim(cx - r, cx + r)
    ax.set_ylim(cy - r, cy + r)
    ax.set_zlim(cz - r, cz + r)


def plot_pl(P: np.ndarray, title: str, path_png: str | None, show: bool = False) -> None:
    """Plot a closed polygonal line and optionally save PNG."""
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111, projection="3d")
    Pc = np.vstack([P, P[:1]])  # close loop for visualization
    ax.plot(Pc[:, 0], Pc[:, 1], Pc[:, 2], lw=1.6)
    ax.scatter(P[:, 0], P[:, 1], P[:, 2], s=6, alpha=0.75)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    set_axes_equal(ax, P)
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    if show:
        plt.show()
    else:
        assert path_png is not None, "path_png must be provided when not showing interactively."
        plt.tight_layout()
        plt.savefig(path_png, dpi=160, bbox_inches="tight")
        plt.close(fig)


def write_xyz(P: np.ndarray, path_xyz: str, close_loop: bool = False) -> None:
    """Write coordinates to a simple XYZ-like text file."""
    Q = np.vstack([P, P[:1]]) if close_loop else P
    with open(path_xyz, "w", encoding="utf-8") as f:
        f.write(f"{len(Q)}\n\n")
        for x, y, z in Q:
            f.write(f"C {x:.6f} {y:.6f} {z:.6f}\n")


def main() -> None:
    # Import here to give a clearer error if missing
    try:
        import topoly  # type: ignore
    except Exception as e:
        raise ImportError(
            "Failed to import 'topoly'. Please install via 'pip install topoly'."
        ) from e

    ap = argparse.ArgumentParser(parents=[_pre])
    ap.add_argument(
        "--length", "--n-edges", dest="n_edges", type=int, required=True,
        help="Number of polygon edges/vertices (>=3)"
    )
    ap.add_argument("--count", type=int, default=10, help="How many samples to generate")
    ap.add_argument(
        "--bond-length", type=float, default=1.0,
        help="Edge length used by TopoLy (default 1.0)"
    )
    ap.add_argument(
        "--coord-dir", type=str, default="output_coordinate",
        help="Directory to save .xyz coordinate files"
    )
    ap.add_argument(
        "--png-dir", type=str, default="output_png",
        help="Directory to save .png figures"
    )
    ap.add_argument(
        "--basename", type=str, default="topoly_random",
        help="Base filename prefix (default: 'topoly_random')"
    )
    args = ap.parse_args()

    if args.n_edges < 3:
        raise ValueError("--length / --n-edges must be >= 3")

    os.makedirs(args.coord_dir, exist_ok=True)
    os.makedirs(args.png_dir, exist_ok=True)

    for i in range(1, args.count + 1):
        # Generate one loop with TopoLy
        loops = topoly.generate_loop(
            args.n_edges, 1, args.bond_length,
            print_with_index=False, output="list"
        )
        P = np.asarray(loops[0], dtype=float)

        base = f"{args.basename}_N{args.n_edges}_{i:03d}"
        path_xyz = os.path.join(args.coord_dir, base + ".xyz")
        path_png = os.path.join(args.png_dir, base + ".png")

        title = f"Random knot (TopoLy) | N={args.n_edges} | bond={args.bond_length} | #{i:03d}"

        # Plot (save or show)
        path_for_plot = path_png if (args.save_png or not args.show) else None
        plot_pl(P, title, path_png=path_for_plot, show=args.show)

        # Save coordinates
        write_xyz(P, path_xyz, close_loop=False)

        # Logs
        print(f"[{i:03d}/{args.count}] XYZ -> {path_xyz}")
        if args.show and not args.save_png:
            print("      (interactive window shown; PNG not saved)")
        else:
            print(f"      PNG -> {path_png}")

    print("[Done] All samples generated.")


if __name__ == "__main__":
    main()
