#!/usr/bin/env python3
"""Render demo figures from a tail-evt JSON output.

Usage:
    python scripts/render_demo.py <tail.json> <out_dir>

Produces:
    <out_dir>/chi_matrix.png   — heatmap of pairwise χ tail coupling
    <out_dir>/gpd_fits.png     — bar plot of per-asset shape ξ and scale β

Kept dependency-light (numpy + matplotlib) so it doubles as a reference for
how the JSON emitted by evt-cli should be consumed downstream.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def render(payload: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    labels = payload["labels"]
    chi = np.asarray(payload["chi_matrix"], dtype=float)

    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    im = ax.imshow(chi, cmap="viridis", vmin=0.0, vmax=max(0.5, float(chi.max())))
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{chi[i, j]:.2f}", ha="center", va="center",
                    color="white" if chi[i, j] < 0.3 else "black", fontsize=9)
    fig.colorbar(im, ax=ax, label="χ(u)")
    ax.set_title(f"Pairwise tail coupling · u = q{int(payload['threshold_quantile']*100)}")
    fig.tight_layout()
    fig.savefig(out_dir / "chi_matrix.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    fits = payload["fits"]
    xs = np.arange(len(fits))
    xi = np.array([f["xi"] for f in fits])
    beta = np.array([f["beta"] for f in fits])
    ax.bar(xs - 0.18, xi, width=0.36, color="#b6ff4a", label="ξ (shape)")
    ax2 = ax.twinx()
    ax2.bar(xs + 0.18, beta, width=0.36, color="#4ec9e0", label="β (scale)")
    ax.set_xticks(xs)
    ax.set_xticklabels([f["asset"] for f in fits])
    ax.set_ylabel("ξ (shape)")
    ax2.set_ylabel("β (scale)")
    ax.set_title("Per-asset GPD fit (PWM estimator)")
    ax.spines[["top"]].set_visible(False)
    ax2.spines[["top"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_dir / "gpd_fits.png", dpi=160)
    plt.close(fig)


def main() -> None:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        sys.exit(2)
    payload = json.loads(Path(sys.argv[1]).read_text())
    render(payload, Path(sys.argv[2]))


if __name__ == "__main__":
    main()
