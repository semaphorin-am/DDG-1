"""Spectral analysis of the Laplace-Beltrami operator (project Phase 7).

Solves the generalized eigenproblem

    L phi_i = lambda_i M phi_i

for the smallest eigenvalues.  On the unit sphere the continuous Laplace-Beltrami
spectrum is ``lambda_l = l (l + 1)`` for ``l = 0, 1, 2, ...`` with multiplicity
``2l + 1``; the eigenfunctions are the spherical harmonics ``Y_l^m``.  The
discrete eigenpairs should converge to these, which makes the sphere an ideal
validation case for spectral geometry.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def eigenpairs(
    L: sp.spmatrix, M: sp.spmatrix, k: int = 20
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the ``k`` smallest generalized eigenpairs of ``(L, M)``.

    Uses shift-invert Lanczos around ``sigma = 0`` (robust for the smallest
    eigenvalues of an SPD-definite pencil).  Eigenvalues are returned in
    ascending order; ``eigenvectors[:, i]`` is the ``i``-th eigenfunction.
    """
    L = L.tocsc()
    M = M.tocsc()
    # sigma slightly below 0 avoids the exactly-singular shifted matrix.
    vals, vecs = spla.eigsh(L, k=k, M=M, sigma=-1e-8, which="LM")
    order = np.argsort(vals)
    return vals[order], vecs[:, order]


def sphere_reference_spectrum(max_l: int) -> np.ndarray:
    """Analytical sphere eigenvalues ``l(l+1)`` repeated with multiplicity
    ``2l+1``, in ascending order, for ``l = 0 .. max_l``."""
    out = []
    for l in range(max_l + 1):
        out.extend([l * (l + 1)] * (2 * l + 1))
    return np.array(out, dtype=float)


def degree_for_index(i: int) -> int:
    """Band index ``l`` such that eigenvalue number ``i`` (0-based) belongs to
    the ``l``-th spherical-harmonic band."""
    l = 0
    count = 0
    while count + (2 * l + 1) <= i:
        count += 2 * l + 1
        l += 1
    return l
