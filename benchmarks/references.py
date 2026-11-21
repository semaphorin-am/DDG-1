"""Independent analytic references for the "geometry beyond the sphere" benchmark.

This module deliberately imports nothing from ``ddg``: every function here is a
closed-form answer computed from coordinates alone, so it can be reviewed (and
re-derived) without reading the solver.  It was committed, together with
``fixtures.json`` and ``tolerances.json``, before any solver change.

Unit sphere
-----------
* geodesic distance      d(p, q) = atan2(|p x q|, p . q)   (great-circle angle)
* heat, harmonic mode    u(t) = exp(-l (l + 1) t) Y_l

Flat periodic torus  T = [0, Lx) x [0, Ly) with opposite edges identified
---------------------------------------------------------------------------
* geodesic distance      d = sqrt(wrap(dx, Lx)^2 + wrap(dy, Ly)^2),
                         wrap(d, L) = min(|d| mod L, L - |d| mod L)
* heat, Fourier mode     u(t) = exp(-lam t) cos(2 pi (k x / Lx + l y / Ly)),
                         lam = 4 pi^2 (k^2 / Lx^2 + l^2 / Ly^2); (0, 0) is constant
* heat, Gaussian bump    periodised Gaussian of variance s0^2 + 2 t
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- sphere
def sphere_distance(points: np.ndarray, source: np.ndarray) -> np.ndarray:
    """Great-circle distance on the unit sphere.

    Uses ``atan2(|p x q|, p . q)``, which is accurate near 0 and near the
    antipode (``arccos`` loses half the digits there), so it is independent of
    the ``arccos`` formula the original benchmark uses.
    """
    p = points / np.linalg.norm(points, axis=1, keepdims=True)
    q = np.asarray(source, dtype=float)
    q = q / np.linalg.norm(q)
    cross = np.linalg.norm(np.cross(p, q[None, :]), axis=1)
    return np.arctan2(cross, p @ q)


def sphere_harmonic(points: np.ndarray, l: int) -> np.ndarray:
    """Zonal harmonic of degree ``l`` (unnormalised Legendre P_l(z))."""
    z = points[:, 2] / np.linalg.norm(points, axis=1)
    if l == 0:
        return np.ones_like(z)
    if l == 1:
        return z
    if l == 2:
        return 0.5 * (3 * z ** 2 - 1)
    if l == 3:
        return 0.5 * (5 * z ** 3 - 3 * z)
    raise ValueError("only l <= 3 implemented")


def sphere_mode_rate(l: int) -> float:
    """Decay rate of a degree-``l`` harmonic under u_t = Laplace u."""
    return float(l * (l + 1))


# --------------------------------------------------------------------------- torus
def wrap(d: np.ndarray, period: float) -> np.ndarray:
    """Shortest absolute displacement on a circle of circumference ``period``."""
    d = np.abs(np.asarray(d, dtype=float)) % period
    return np.minimum(d, period - d)


def torus_distance(points: np.ndarray, source: np.ndarray, Lx: float, Ly: float) -> np.ndarray:
    """Exact geodesic distance on the flat torus (wraps through both seams)."""
    dx = wrap(points[:, 0] - source[0], Lx)
    dy = wrap(points[:, 1] - source[1], Ly)
    return np.hypot(dx, dy)


def torus_distance_images(points: np.ndarray, source: np.ndarray, Lx: float, Ly: float,
                          n_images: int = 2) -> np.ndarray:
    """Same distance by brute force: minimum over periodic images of the source.

    A second, structurally different derivation used to cross-check
    :func:`torus_distance` in the tests.
    """
    best = np.full(points.shape[0], np.inf)
    for a in range(-n_images, n_images + 1):
        for b in range(-n_images, n_images + 1):
            d = np.hypot(points[:, 0] - source[0] - a * Lx, points[:, 1] - source[1] - b * Ly)
            best = np.minimum(best, d)
    return best


def unwrapped_distance(points: np.ndarray, source: np.ndarray) -> np.ndarray:
    """KNOWN-WRONG answer for the torus: straight-line distance inside the
    rectangle, ignoring the identification of opposite edges.  Used only as a
    negative control to show the acceptance checks can fail."""
    return np.hypot(points[:, 0] - source[0], points[:, 1] - source[1])


def torus_diameter(Lx: float, Ly: float) -> float:
    """Largest possible geodesic distance on the flat torus."""
    return float(0.5 * np.hypot(Lx, Ly))


def torus_mode(points: np.ndarray, k: int, l: int, Lx: float, Ly: float) -> np.ndarray:
    return np.cos(2 * np.pi * (k * points[:, 0] / Lx + l * points[:, 1] / Ly))


def torus_mode_rate(k: int, l: int, Lx: float, Ly: float) -> float:
    return float(4 * np.pi ** 2 * (k ** 2 / Lx ** 2 + l ** 2 / Ly ** 2))


def torus_gaussian(points: np.ndarray, center, s0: float, t: float, Lx: float, Ly: float,
                   n_images: int = 6) -> np.ndarray:
    """Exact solution of u_t = Laplace u on the flat torus whose initial state
    is a unit-mass periodised Gaussian of standard deviation ``s0``.

    The heat flow of a Gaussian is a Gaussian with variance ``s0^2 + 2 t``;
    periodising commutes with the flow, so the solution is the image sum.
    """
    var = s0 ** 2 + 2.0 * t
    gx = np.zeros(points.shape[0])
    gy = np.zeros(points.shape[0])
    for a in range(-n_images, n_images + 1):
        gx += np.exp(-((points[:, 0] - center[0] - a * Lx) ** 2) / (2 * var))
        gy += np.exp(-((points[:, 1] - center[1] - a * Ly) ** 2) / (2 * var))
    return gx * gy / (2 * np.pi * var)
