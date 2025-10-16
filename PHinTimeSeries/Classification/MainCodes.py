#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: Maincodes.py
Author: Yuzhou He <ribosomehyz@gmail.com>

Purpose:
  Demo pipeline for time-delay embedding (Takens) with automatic parameter
  search (giotto-tda), 3D visualization (direct or PCA), and Vietoris–Rips
  persistent homology up to H2.

Usage:
  # minimal (defaults reproduce the example)
  python Maincodes.py

  # custom stride / search limits / coeff field / jobs
  python Maincodes.py --stride 4 --max-dim 30 --max-delay 30 --coeff 2 --jobs 6

  # change signal length and sampling
  python Maincodes.py --x-max 50 --n-samples 1001

Requirements:
  - Python 3.9+
  - numpy
  - scikit-learn
  - giotto-tda     (pip install giotto-tda)
  - plotly         (giotto-tda plotting backend; usually installed with giotto-tda)

Notes:
  - plot_point_cloud returns a Plotly Figure. Use .show() for an interactive view,
    or export with fig.write_html(...).
"""

from __future__ import annotations

import argparse
import numpy as np
from sklearn.decomposition import PCA
from gtda.time_series import SingleTakensEmbedding
from gtda.plotting import plot_point_cloud
from gtda.homology import VietorisRipsPersistence


# -------------------- core helpers --------------------
def build_signal(x_max: float, n_samples: int) -> tuple[np.ndarray, np.ndarray]:
    """Toy periodic signal y(x) = cos(5x) + cos(pi x)."""
    x = np.linspace(0.0, float(x_max), int(n_samples))
    y = np.cos(5.0 * x) + np.cos(np.pi * x)
    return x, y


def fit_embedder(embedder: SingleTakensEmbedding, y: np.ndarray, verbose: bool = True) -> np.ndarray:
    """Fit Takens embedder and print the discovered (dimension, time_delay)."""
    Y = embedder.fit_transform(y)
    if verbose:
        print(f"Shape of embedded time series: {Y.shape}")
        print(f"Optimal embedding dimension is {embedder.dimension_} and time delay is {embedder.time_delay_}")
    return Y


def make_embedding(y: np.ndarray, stride: int, max_dim: int, max_delay: int) -> tuple[np.ndarray, int, int]:
    """Search optimal (d, tau), then re-embed with fixed parameters."""
    # search
    searcher = SingleTakensEmbedding(
        parameters_type="search",
        time_delay=max_delay,
        dimension=max_dim,
        stride=stride,
    )
    _ = fit_embedder(searcher, y, verbose=True)

    # fixed with discovered params
    d_opt = int(searcher.dimension_)
    tau_opt = int(searcher.time_delay_)
    embedder = SingleTakensEmbedding(
        parameters_type="fixed",
        n_jobs=2,
        time_delay=tau_opt,
        dimension=d_opt,
        stride=stride,
    )
    Y = embedder.fit_transform(y)
    return Y, d_opt, tau_opt


# -------------------- CLI main --------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="Takens embedding + VR persistence demo (giotto-tda).")
    ap.add_argument("--x-max", type=float, default=50.0, help="Right endpoint of x domain")
    ap.add_argument("--n-samples", type=int, default=1001, help="Number of samples in the signal")
    ap.add_argument("--stride", type=int, default=4, help="Stride for Takens embedding")
    ap.add_argument("--max-dim", type=int, default=30, help="Max embedding dimension for search")
    ap.add_argument("--max-delay", type=int, default=30, help="Max time delay for search")
    ap.add_argument("--coeff", type=int, default=2, help="Coefficient field for homology")
    ap.add_argument("--jobs", type=int, default=6, help="Parallel jobs for VR persistence")
    args = ap.parse_args()

    # 1) build toy signal
    x, y = build_signal(args.x_max, args.n_samples)

    # 2) Takens embedding (search -> fixed)
    Y, d_opt, tau_opt = make_embedding(y, stride=args.stride, max_dim=args.max_dim, max_delay=args.max_delay)

    # 3) 3D visualization (direct if d=3 else PCA to 3D)
    if d_opt == 3:
        fig_embed = plot_point_cloud(Y)
    else:
        pca = PCA(n_components=3)
        Y3 = pca.fit_transform(Y)
        fig_embed = plot_point_cloud(Y3)

    print(f"Embedding ready: d*={d_opt}, tau*={tau_opt}. Opening 3D view...")
    fig_embed.show()

    # 4) Vietoris–Rips persistence (H0/H1/H2)
    batch = Y[None, :, :]  # shape (1, T, d)
    vr = VietorisRipsPersistence(
        homology_dimensions=[0, 1, 2],
        coeff=int(args.coeff),
        n_jobs=int(args.jobs),
    )
    print("Computing persistence diagram(s)...")
    vr.fit_transform_plot(batch)   # opens an interactive diagram figure


if __name__ == "__main__":
    main()
