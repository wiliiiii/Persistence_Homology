#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: Quantities.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Build a toy periodic signal, apply a fixed Takens (time-delay) embedding,
  compute Vietoris–Rips persistent homology up to H2, and report bar-length
  statistics (L2, L3, RMS, max, mean, count) for H1 and H2.

Usage:
  # reproduce the example in the paper/snippet
  python Quantities.py

  # customize embedding / signal / homology / parallelism
  python Quantities.py --dim 6 --tau 12 --stride 4 --x-max 50 --n-samples 1001 \
                       --coeff 2 --jobs 6 --hom-dims 0 1 2

Requirements:
  - Python 3.9+
  - numpy
  - giotto-tda   (SingleTakensEmbedding, VietorisRipsPersistence)

Notes:
  - Statistics ignore infinite bars (death = inf).
  - RMS is defined as L2 / sqrt(count), with count≥1 guarded.
"""

from __future__ import annotations

import argparse
import numpy as np
from gtda.time_series import SingleTakensEmbedding
from gtda.homology import VietorisRipsPersistence


# ---------------------- helpers ----------------------
def build_signal(x_max: float, n_samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Toy periodic 1D signal: cos(5x) + sin(7x)."""
    x = np.linspace(0.0, float(x_max), int(n_samples))
    y = np.cos(5.0 * x) + np.sin(7.0 * x)
    return x, y


def takens_embed(y: np.ndarray, d: int, tau: int, stride: int, n_jobs: int = 2) -> np.ndarray:
    """Fixed-parameter Takens embedding."""
    embedder = SingleTakensEmbedding(
        parameters_type="fixed",
        n_jobs=int(n_jobs),
        time_delay=int(tau),
        dimension=int(d),
        stride=int(stride),
    )
    Y = embedder.fit_transform(y)
    return Y


def stats_for_dim(diag_sample: np.ndarray, dim: int) -> tuple[float, float, float, float, int]:
    """
    Compute bar-length statistics for a given homology dimension.
    diag_sample: (n_points, 3) with columns (birth, death, hom_dim).
    Returns: (L2, L3, max, mean, count). Infinite bars are ignored.
    """
    D = diag_sample[diag_sample[:, 2] == dim, :2]  # (birth, death)
    D = D[np.isfinite(D[:, 1])]
    if D.size == 0:
        return 0.0, 0.0, 0.0, 0.0, 0
    pers = D[:, 1] - D[:, 0]
    l2 = float(np.linalg.norm(pers, ord=2))            # (sum pers^2)^(1/2)
    l3 = float(np.sum(pers ** 3) ** (1.0 / 3.0))       # (sum pers^3)^(1/3)
    pmax = float(np.max(pers))
    pmean = float(np.mean(pers))
    return l2, l3, pmax, pmean, int(pers.size)


# ------------------------ main -----------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Takens embedding + VR PH statistics (H1/H2).")
    # signal
    ap.add_argument("--x-max", type=float, default=50.0, help="Right-end of x domain")
    ap.add_argument("--n-samples", type=int, default=1001, help="Number of samples of the signal")
    # embedding
    ap.add_argument("--dim", type=int, default=6, help="Takens embedding dimension")
    ap.add_argument("--tau", type=int, default=12, help="Takens time delay")
    ap.add_argument("--stride", type=int, default=4, help="Stride for Takens embedding")
    ap.add_argument("--embed-jobs", type=int, default=2, help="Parallel jobs for embedding")
    # homology
    ap.add_argument("--hom-dims", type=int, nargs="+", default=[0, 1, 2], help="Homology dimensions")
    ap.add_argument("--coeff", type=int, default=2, help="Coefficient field for PH")
    ap.add_argument("--jobs", type=int, default=6, help="Parallel jobs for PH")
    args = ap.parse_args()

    # 1) signal
    _, y = build_signal(args.x_max, args.n_samples)

    # 2) fixed Takens embedding
    Y = takens_embed(y, d=args.dim, tau=args.tau, stride=args.stride, n_jobs=args.embed_jobs)

    # 3) Vietoris–Rips persistence
    batch = Y[None, :, :]  # shape: (1, T, d)
    vr = VietorisRipsPersistence(
        homology_dimensions=list(args.hom_dims),
        coeff=int(args.coeff),
        n_jobs=int(args.jobs),
    )
    diagrams = vr.fit_transform(batch)  # shape: (1, n_points, 3)
    diag0 = diagrams[0]

    # 4) statistics for H1 and H2 (if requested)
    for target_dim in [1, 2]:
        if target_dim not in args.hom_dims:
            print(f"H{target_dim} not computed (omit from --hom-dims).")
            continue
        l2, l3, pmax, pmean, cnt = stats_for_dim(diag0, dim=target_dim)
        rms = l2 / np.sqrt(max(cnt, 1))
        print(f"H{target_dim} ({target_dim}D classes):")
        print(f"  count = {cnt},  L2 = {l2:.6f},  L3 = {l3:.6f},  RMS = {rms:.6f},  "
              f"max = {pmax:.6f},  mean = {pmean:.6f}")

if __name__ == "__main__":
    main()
