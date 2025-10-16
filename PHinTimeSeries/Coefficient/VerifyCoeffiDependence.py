#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: VerifyCoeffiDependence.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Test coefficient-field dependence of Vietoris–Rips persistence diagrams.
  1. Build a synthetic periodic signal.
  2. Use Takens embedding (with automatic search for optimal dimension/time delay).
  3. Compute persistence diagrams (H1, H2) with field Z_p for various primes.
  4. Compare PDs across primes (e.g., p ≤ 13) to detect coefficient dependence.

Usage:
  python VerifyCoeffiDependence.py
  # optional arguments:
  #   --max-dim 30 --max-delay 30 --stride 3
  #   --coeff 7 --jobs 6
  #   --primes 2 3 5 7 11 13
  #   --tol 1e-9

Requirements:
  - Python 3.9+
  - numpy
  - giotto-tda (for SingleTakensEmbedding, VietorisRipsPersistence)

Notes:
  - The script prints summary statistics of H1/H2 bars and tests PD equality
    across coefficient primes p ≤ 13.
"""

import sys
import numpy as np
from gtda.time_series import SingleTakensEmbedding
from gtda.homology import VietorisRipsPersistence


# -------------------- signal & embedding --------------------
def build_signal(xmax: float = 50.0, n_samples: int = 1001) -> np.ndarray:
    """Create synthetic periodic signal y(x) = 0.6 cos(x) + 0.8 sin(pi x)."""
    x = np.linspace(0, xmax, n_samples)
    return 0.6 * np.cos(x) + 0.8 * np.sin(np.pi * x)


def fit_embedder(embedder: SingleTakensEmbedding, y: np.ndarray, verbose: bool = True) -> np.ndarray:
    """Fit Takens embedder and print optimal dimension and time delay."""
    y_embedded = embedder.fit_transform(y)
    if verbose:
        print(f"Shape of embedded time series: {y_embedded.shape}")
        print(f"Optimal embedding dimension = {embedder.dimension_}, time delay = {embedder.time_delay_}")
    return y_embedded


def takens_embedding(y: np.ndarray, max_dim: int = 30, max_delay: int = 30, stride: int = 3) -> np.ndarray:
    """Search optimal Takens parameters, then re-embed with fixed ones."""
    embedder_search = SingleTakensEmbedding(
        parameters_type="search",
        time_delay=max_delay,
        dimension=max_dim,
        stride=stride,
    )
    fit_embedder(embedder_search, y)
    embedder_fixed = SingleTakensEmbedding(
        parameters_type="fixed",
        n_jobs=2,
        time_delay=embedder_search.time_delay_,
        dimension=embedder_search.dimension_,
        stride=stride,
    )
    return embedder_fixed.fit_transform(y)


# -------------------- PH + statistics --------------------
def stats_for_dim(diag_sample: np.ndarray, dim: int):
    """
    Compute statistics for given homology dimension.
    Returns (L2, L3, max, mean, count, RMS).
    """
    D = diag_sample[diag_sample[:, 2] == dim, :2]
    D = D[np.isfinite(D[:, 1])]
    if D.size == 0:
        return 0.0, 0.0, 0.0, 0.0, 0, 0.0
    pers = D[:, 1] - D[:, 0]
    l2 = float(np.linalg.norm(pers, 2))
    l3 = float(np.sum(pers**3) ** (1.0 / 3.0))
    pmax = float(pers.max())
    pmean = float(pers.mean())
    n = pers.size
    rms = l2 / np.sqrt(n)
    return l2, l3, pmax, pmean, n, rms


def compute_diagram(points, coeff=7, homology_dimensions=(0, 1, 2), n_jobs=6):
    """Compute Vietoris–Rips persistence diagram for given coefficient field."""
    vr = VietorisRipsPersistence(
        homology_dimensions=list(homology_dimensions),
        coeff=coeff,
        n_jobs=n_jobs,
    )
    return vr.fit_transform(points[None, :, :])[0]  # shape: (n_pts, 3)


# -------------------- PD comparison tools --------------------
def PD_equal_H1_across_primes(points, primes=(2, 3, 5, 7, 11, 13), tol=1e-9, n_jobs=6):
    """
    Check if H1 persistence diagrams are identical across primes.
    Returns (True/False, first differing prime if any).
    """
    base = None
    X = points[None, :, :]
    for p in primes:
        vr = VietorisRipsPersistence(homology_dimensions=[1], coeff=p, n_jobs=n_jobs)
        D = vr.fit_transform(X)[0]
        D = D[np.isfinite(D[:, 1])]
        D = D[D[:, 2] == 1][:, :2]
        D = D[np.lexsort((D[:, 1], D[:, 0]))]
        if base is None:
            base = D
        else:
            same_shape = (D.shape == base.shape)
            same_vals = same_shape and (np.max(np.abs(D - base)) <= tol if D.size else True)
            if not same_vals:
                return False, p
    return True, None


def PD_H1_differences_across_primes(points, primes=(2, 3, 5, 7, 11, 13), tol=1e-9, n_jobs=6):
    """
    Compare H1 persistence diagrams across primes vs. base prime (first one).
    Returns (diff_primes, details, base_prime).
    details: list of (p, n_bars_p, n_bars_base, max_abs_diff_if_same_len_or_None)
    """
    X = points[None, :, :]

    def pd_H1_for_prime(p):
        vr = VietorisRipsPersistence(homology_dimensions=[1], coeff=p, n_jobs=n_jobs)
        D = vr.fit_transform(X)[0]
        D = D[np.isfinite(D[:, 1])]
        D = D[D[:, 2] == 1][:, :2]
        return D[np.lexsort((D[:, 1], D[:, 0]))] if D.size else D

    PDs = {p: pd_H1_for_prime(p) for p in primes}
    base_p = primes[0]
    base_D = PDs[base_p]

    diff_primes, details = [], []
    for p in primes[1:]:
        D = PDs[p]
        same_shape = (D.shape == base_D.shape)
        if same_shape:
            if D.size == 0 and base_D.size == 0:
                same_vals, max_diff = True, 0.0
            else:
                max_diff = float(np.max(np.abs(D - base_D))) if D.size else 0.0
                same_vals = (max_diff <= tol)
        else:
            same_vals, max_diff = False, None
        if not same_vals:
            diff_primes.append(p)
        details.append((p, D.shape[0], base_D.shape[0], max_diff))
    return diff_primes, details, base_p


# -------------------- main routine --------------------
def main():
    print("Using Python:", sys.executable)
    y = build_signal()
    Y = takens_embedding(y)

    # compute PH with coeff=7
    diag = compute_diagram(Y, coeff=7, homology_dimensions=(1, 2))
    l2_H1, l3_H1, max_H1, mean_H1, n1, rms_H1 = stats_for_dim(diag, dim=1)
    l2_H2, l3_H2, max_H2, mean_H2, n2, rms_H2 = stats_for_dim(diag, dim=2)

    print("H1 (1D classes):")
    print(f"  count={n1}, L2={l2_H1:.6f}, L3={l3_H1:.6f}, RMS={rms_H1:.6f}, "
          f"max={max_H1:.6f}, mean={mean_H1:.6f}")
    print("H2 (2D classes):")
    print(f"  count={n2}, L2={l2_H2:.6f}, L3={l3_H2:.6f}, RMS={rms_H2:.6f}, "
          f"max={max_H2:.6f}, mean={mean_H2:.6f}")

    same, badp = PD_equal_H1_across_primes(Y)
    print("H1 PD equal across primes ≤13?:", same, "| first differing prime:", badp)

    diff_ps, info, base_p = PD_H1_differences_across_primes(Y)
    if not diff_ps:
        print(f"H1 PD identical for all primes ≤13 (base p={base_p}).")
    else:
        print(f"H1 PD differ (vs p={base_p}) at primes:", diff_ps)
        for p, n_p, n_base, maxdiff in info:
            if p in diff_ps:
                if maxdiff is None:
                    print(f"  p={p}: bar count differs (n_p={n_p}, n_base={n_base})")
                else:
                    print(f"  p={p}: same bar count ({n_p}), max |Δ|={maxdiff:.3e}")


if __name__ == "__main__":
    main()
