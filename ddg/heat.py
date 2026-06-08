"""Heat diffusion on a surface (project Phase 4).

Solves the heat equation

    du/dt = Delta u

on the mesh using a backward (implicit) Euler step, which is unconditionally
stable.  With the positive-semidefinite stiffness ``L`` and the lumped mass
matrix ``M`` (see :mod:`ddg.geometry`), the implicit step is the linear system

    (M + dt * L) u^{n+1} = M u^n.

Because ``M + dt*L`` is symmetric positive-definite and constant in time, we
prefactor it once with a sparse Cholesky (via ``scipy``'s supernodal LU on the
SPD matrix) and reuse the factorization for every step.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class HeatSolver:
    """Reusable implicit-Euler heat solver for a fixed mesh / timestep."""

    def __init__(self, M: sp.spmatrix, L: sp.spmatrix, dt: float):
        self.M = M.tocsc()
        self.L = L.tocsc()
        self.dt = float(dt)
        A = (self.M + self.dt * self.L).tocsc()
        # SPD system -> prefactor once, reuse across timesteps.
        self._solve = spla.factorized(A)

    def step(self, u: np.ndarray) -> np.ndarray:
        """Advance the field ``u`` by one timestep ``dt``."""
        rhs = self.M @ u
        return self._solve(rhs)

    def run(
        self,
        u0: np.ndarray,
        n_steps: int,
        record_every: int = 1,
    ) -> List[np.ndarray]:
        """Integrate ``n_steps`` steps, recording snapshots.

        Returns a list of snapshots including the initial condition.
        """
        u = np.array(u0, dtype=np.float64)
        snaps = [u.copy()]
        for n in range(1, n_steps + 1):
            u = self.step(u)
            if n % record_every == 0:
                snaps.append(u.copy())
        return snaps


def point_source(n_vertices: int, source: int) -> np.ndarray:
    """Unit heat concentrated at a single vertex (a Kronecker delta)."""
    u0 = np.zeros(n_vertices)
    u0[source] = 1.0
    return u0


def total_heat(M: sp.spmatrix, u: np.ndarray) -> float:
    r"""Integrated heat ``\int u dA = 1^T M u`` -- conserved by the continuous
    flow on a closed surface; a useful diagnostic."""
    return float(np.asarray(M @ u).ravel().sum())
