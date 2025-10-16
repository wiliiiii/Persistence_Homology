#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: ACN.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Compute Average Crossing Number (ACN) for batches of .xyz polylines using
  pyknotid, and append the results to existing per-N Excel files.

Usage:
  # Windows path must be quoted and keep backslashes
  python ACN.py \
    --src "C:\\Users\\LENOVO\\Desktop\\UoR\\Paper_Topoly\\PLknot" \
    --dst "C:\\Users\\LENOVO\\Desktop\\UoR\\Paper_Topoly\\result" \
    --Ns 10 20 30 40 50 60 70 80 90 100 \
    --samples 400

  # Optional flags:
  #   --xyz-glob "*.xyz"
  #   --xyz-subdir "output_coordinate"
  #   --sheet-name "results"
  #   --col-name "ACN"

Requirements:
  - Python 3.9+
  - numpy
  - pandas
  - openpyxl       (pandas Excel engine)
  - pyknotid       (pip install pyknotid)

Notes:
  - For each N, this script expects:
      <src>/N=<N>/<xyz-subdir>/*.xyz
    and writes into:
      <dst>/N=<N>.xlsx  (sheet=<sheet-name>)
  - The Excel index is assumed to be filenames; the script aligns by index.
"""

# ---- NumPy 1.24+ compatibility shims for deprecated aliases ----
import numpy as _np
if not hasattr(_np, "float"):   _np.float = float   # type: ignore[attr-defined]
if not hasattr(_np, "int"):     _np.int = int       # type: ignore[attr-defined]
if not hasattr(_np, "complex"): _np.complex = complex  # type: ignore[attr-defined]
if not hasattr(_np, "bool"):    _np.bool = bool     # type: ignore[attr-defined]

import os
import glob
import argparse
import numpy as np
import pandas as pd
from pyknotid.spacecurves import Knot, SpaceCurve


# ---------- I/O ----------
def load_xyz_xyzlike(path: str) -> np.ndarray:
    """Read last 3 columns as (x,y,z); skip first two .xyz header lines if present."""
    coords = []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    start = 2 if len(lines) >= 2 else 0
    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) >= 3:
            coords.append(tuple(map(float, parts[-3:])))
    P = np.asarray(coords, dtype=float)
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"Cannot parse 3D points from {path}")
    return P


def make_curve(P: np.ndarray, tol: float = 1e-9):
    """Return Knot if closed; otherwise SpaceCurve."""
    if P.shape[0] >= 2 and np.allclose(P[0], P[-1], atol=tol, rtol=0):
        return Knot(P)
    return SpaceCurve(P)


# ---------- Core ----------
def compute_acn_for_folder(folder: str, xyz_glob: str, samples: int) -> dict[str, float]:
    """Compute ACN for all .xyz files under folder."""
    rows: dict[str, float] = {}
    files = sorted(glob.glob(os.path.join(folder, xyz_glob)))
    if not files:
        print(f"[WARN] No xyz files found in: {folder}")
        return rows
    for i, path in enumerate(files, 1):
        name = os.path.basename(path)
        P = load_xyz_xyzlike(path)
        curve = make_curve(P)
        acn = float(curve.average_crossing_number(samples=samples))
        rows[name] = acn
        print(f"  [{i:03d}/{len(files)}] {name}  ACN={acn:.6f}")
    return rows


def read_results_excel(xlsx_path: str, sheet_name: str):
    """Read DataFrame from Excel; try sheet_name then fall back to first sheet."""
    if not os.path.isfile(xlsx_path):
        raise FileNotFoundError(f"Excel not found: {xlsx_path}")
    try:
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name, index_col=0)
        used = sheet_name
    except Exception:
        df = pd.read_excel(xlsx_path, sheet_name=0, index_col=0)
        used = 0
        print(f"[INFO] Sheet '{sheet_name}' not found in {os.path.basename(xlsx_path)}; fallback to the first sheet.")
    return df, used


def write_back_excel(xlsx_path: str, sheet_used, df: pd.DataFrame, col_name: str) -> None:
    """Write updated df back to Excel, placing col_name at the end."""
    cols = [c for c in df.columns if c != col_name] + [col_name]
    df = df[cols]
    with pd.ExcelWriter(xlsx_path, engine="openpyxl", mode="w") as w:
        sname = sheet_used if isinstance(sheet_used, str) else "results"
        df.to_excel(w, sheet_name=sname, index=True)
    print(f"[OK] Updated -> {xlsx_path}  (sheet='{sname}')")


# ---------- CLI ----------
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Compute ACN for .xyz batches and append to Excel.")
    ap.add_argument("--src", type=str, required=True, help="Root directory of source data (contains N=<N>/...)")
    ap.add_argument("--dst", type=str, required=True, help="Root directory of Excel outputs (contains N=<N>.xlsx)")
    ap.add_argument("--Ns", type=int, nargs="+", required=True, help="List of N values, e.g. 10 20 ... 100")
    ap.add_argument("--samples", type=int, default=400, help="Random projections for ACN estimation")
    ap.add_argument("--xyz-glob", type=str, default="*.xyz", help="Filename glob for .xyz")
    ap.add_argument("--xyz-subdir", type=str, default="output_coordinate", help="Subdir under N=<N> containing xyz files")
    ap.add_argument("--sheet-name", type=str, default="results", help="Excel sheet to read/write")
    ap.add_argument("--col-name", type=str, default="ACN", help="Column name for ACN")
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    for N in args.Ns:
        src_folder = os.path.join(args.src, f"N={N}", args.xyz_subdir)
        dst_excel = os.path.join(args.dst, f"N={N}.xlsx")

        print(f"\n=== N={N} ===")
        if not os.path.isdir(src_folder):
            print(f"[WARN] Source folder missing: {src_folder}")
            continue
        if not os.path.isfile(dst_excel):
            print(f"[WARN] Target Excel missing: {dst_excel}")
            continue

        # 1) compute ACN
        acn_map = compute_acn_for_folder(src_folder, args.xyz_glob, args.samples)
        acn_series = pd.Series(acn_map, name=args.col_name)

        # 2) read Excel
        df, used_sheet = read_results_excel(dst_excel, args.sheet_name)

        # 3) align by filename index and update column
        union_idx = df.index.union(acn_series.index)
        df = df.reindex(union_idx)
        df[args.col_name] = acn_series

        # 4) write back
        write_back_excel(dst_excel, used_sheet, df, args.col_name)


if __name__ == "__main__":
    main()
