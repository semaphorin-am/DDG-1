"""Error metrics and convergence analysis (project Phase 6).

Compares numerical geodesic distances against the exact unit-sphere geodesic
``d(p, q) = arccos(p . q)`` and reports MAE, L-infinity and RMSE, plus an
empirical convergence rate fitted from a refinement sequence.
"""

from __future__ import annotations

from typing import Dict, Sequence

import numpy as np


def exact_sphere_distance(points: np.ndarray, source: np.ndarray) -> np.ndarray:
    r"""Exact geodesic distance on the unit sphere from ``source`` to every row
    of ``points``: ``d = arccos(p . q)`` (great-circle / central angle)."""
    p = points / np.linalg.norm(points, axis=1, keepdims=True)
    q = source / np.linalg.norm(source)
    dots = np.clip(p @ q, -1.0, 1.0)
    return np.arccos(dots)


def error_metrics(numeric: np.ndarray, exact: np.ndarray) -> Dict[str, float]:
    """Return MAE, L-infinity (max) and RMSE between two distance fields."""
    diff = numeric - exact
    return {
        "MAE": float(np.mean(np.abs(diff))),
        "Linf": float(np.max(np.abs(diff))),
        "RMSE": float(np.sqrt(np.mean(diff ** 2))),
    }


def convergence_rate(h: Sequence[float], err: Sequence[float]) -> float:
    """Fit ``err ~ C * h**p`` by least squares in log-log space; return ``p``.

    ``h`` is a representative mesh size (e.g. mean edge length) and ``err`` the
    corresponding error metric.  A slope ``p ~ 1`` indicates first-order and
    ``p ~ 2`` second-order convergence under refinement.
    """
    h = np.asarray(h, dtype=float)
    err = np.asarray(err, dtype=float)
    mask = (h > 0) & (err > 0)
    slope, _ = np.polyfit(np.log(h[mask]), np.log(err[mask]), 1)
    return float(slope)
