"""Geodesic distance via the Heat Method (project Phase 5).

Implements Crane, Weischedel & Wardetzky, *"Geodesics in Heat"* (ACM TOG 2013).

The algorithm recovers the geodesic distance ``phi`` from a source set in three
linear-algebra steps:

1. **Heat flow.**  Integrate ``u_t = Delta u`` for a short time ``t`` from a
   delta heat source:  ``(M + t L) u = delta``.
2. **Normalized gradient.**  Form the per-face gradient ``grad u`` and reverse /
   normalize it:  ``X = -grad u / |grad u|``.  Where heat has barely reached,
   the level sets of ``u`` already look like distance contours, so the unit
   field ``X`` points along geodesics away from the source.
3. **Poisson recovery.**  Solve ``L phi = div X`` for the scalar field whose
   gradient best matches ``X``, then shift so the source has distance 0.

Only the right-hand side changes with the source, so both ``(M + t L)`` and the
Poisson matrix ``L`` are prefactored and reused.
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .mesh import TriMesh
from .ordering import nested_dissection
from .geometry import (
    face_areas,
    face_normals,
    interior_angles,
    cotangent_laplacian,
    mass_matrix,
)


# ---------------------------------------------------------------------------
# Differential operators used by the heat method
# ---------------------------------------------------------------------------
def face_gradient(
    mesh: TriMesh,
    u: np.ndarray,
    areas: Optional[np.ndarray] = None,
    normals: Optional[np.ndarray] = None,
) -> np.ndarray:
    r"""Per-face gradient of a vertex scalar field ``u``.

    For a triangle ``(i, j, k)`` with unit normal ``N`` and area ``A``,

        grad u = 1/(2A) * sum_m  u_m * ( N x e_m ),

    where ``e_m`` is the edge *opposite* the ``m``-th vertex.  Returns an
    ``(|F|, 3)`` array of gradient vectors lying in each triangle's plane.
    """
    c = mesh.corner_positions()
    f = mesh.faces
    if areas is None:
        areas = face_areas(mesh)
    if normals is None:
        normals = face_normals(mesh)

    grad = np.zeros((mesh.n_faces, 3))
    for m in range(3):
        # edge opposite local vertex m goes from vertex (m+1) to (m+2)
        e_opp = c[:, (m + 2) % 3] - c[:, (m + 1) % 3]
        grad += u[f[:, m]][:, None] * np.cross(normals, e_opp)
    grad /= (2.0 * areas)[:, None]
    return grad


def integrated_divergence(
    mesh: TriMesh,
    X: np.ndarray,
    cot: Optional[np.ndarray] = None,
) -> np.ndarray:
    r"""Integrated divergence of a per-face vector field at each vertex.

        (div X)_i = 1/2 * sum_{f ni i} [ cot(t1) (e1 . X_f) + cot(t2) (e2 . X_f) ],

    where ``e1, e2`` are the two edges of face ``f`` emanating from vertex ``i``
    and ``t1, t2`` are the angles opposite them.  Returns a ``(|V|,)`` array,
    the right-hand side of the Poisson system.
    """
    c = mesh.corner_positions()
    f = mesh.faces
    if cot is None:
        cot = 1.0 / np.tan(interior_angles(mesh))

    div = np.zeros(mesh.n_vertices)
    for m in range(3):
        i = f[:, m]
        e1 = c[:, (m + 1) % 3] - c[:, m]   # edge i->j, opposite angle at vertex k
        e2 = c[:, (m + 2) % 3] - c[:, m]   # edge i->k, opposite angle at vertex j
        cot_k = cot[:, (m + 2) % 3]
        cot_j = cot[:, (m + 1) % 3]
        contrib = 0.5 * (
            cot_k * np.einsum("ij,ij->i", e1, X)
            + cot_j * np.einsum("ij,ij->i", e2, X)
        )
        np.add.at(div, i, contrib)
    return div


# ---------------------------------------------------------------------------
# Pinned Poisson solve for the singular closed-surface Laplacian
# ---------------------------------------------------------------------------
def _pinned_solve(L: sp.csr_matrix, b: np.ndarray, pin: int = 0) -> np.ndarray:
    """Solve the consistent singular system ``L phi = b`` (``L @ 1 = 0``).

    The constant nullspace is removed by pinning ``phi[pin] = 0``: drop that row
    and column, solve the resulting symmetric positive-definite system, and
    reinsert the pinned value.  The global additive constant is irrelevant here
    because distances are shifted to zero at the source afterwards.
    """
    n = L.shape[0]
    keep = np.ones(n, dtype=bool)
    keep[pin] = False
    Lr = L[keep][:, keep].tocsc()
    br = b[keep]
    phi_r = spla.spsolve(Lr, br)
    phi = np.zeros(n)
    phi[keep] = phi_r
    return phi


# ---------------------------------------------------------------------------
# The heat method
# ---------------------------------------------------------------------------
class HeatMethodSolver:
    """Prefactored heat-method solver for geodesic distance on a fixed mesh.

    Parameters
    ----------
    mesh : TriMesh
    m : float
        Time-scale multiplier; the diffusion time is ``t = m * h_mean**2`` with
        ``h_mean`` the average edge length (Crane et al. recommend ``m ~ 1``).
    ordering : {"colamd", "nested"}
        Elimination order for the two sparse factorizations.  ``"colamd"`` is
        the original behaviour (SuperLU defaults).  ``"nested"`` reorders the
        vertices by nested dissection (:mod:`ddg.ordering`) and factors without
        pivoting; same answers to round-off, less fill.
    """

    def __init__(self, mesh: TriMesh, m: float = 1.0, ordering: str = "colamd"):
        if ordering not in ("colamd", "nested"):
            raise ValueError("ordering must be 'colamd' or 'nested'")
        self.mesh = mesh
        self.ordering = ordering
        self.areas = face_areas(mesh)
        self.normals = face_normals(mesh)
        self.angles = interior_angles(mesh)
        self.cot = 1.0 / np.tan(self.angles)

        self.L, _ = cotangent_laplacian(mesh)
        self.M = mass_matrix(mesh, self.areas)

        h = mesh.edge_lengths().mean()
        self.t = m * h * h

        A = (self.M + self.t * self.L).tocsc()
        # Poisson matrix with vertex 0 pinned (kept fixed for reuse).
        self._pin = 0
        n = mesh.n_vertices
        keep = np.ones(n, dtype=bool)
        keep[self._pin] = False
        self._keep = keep

        # Prefactor the two SPD-ish systems used on every query.
        if ordering == "colamd":
            heat_lu = spla.splu(A)
            pois_lu = spla.splu(self.L[keep][:, keep].tocsc())
            self._heat_solve = heat_lu.solve
            self._Lr_solve = pois_lu.solve
        else:
            perm = nested_dissection(mesh.vertices, mesh.edges)
            perm_r = perm[perm != self._pin]      # same order, pin removed
            opts = dict(permc_spec="NATURAL", diag_pivot_thresh=0.0,
                        options=dict(SymmetricMode=True))
            heat_lu = spla.splu(A[perm][:, perm].tocsc(), **opts)
            pois_lu = spla.splu(self.L[perm_r][:, perm_r].tocsc(), **opts)
            inv = np.empty(n, dtype=np.int64)
            inv[perm] = np.arange(n)
            # position of each kept vertex inside perm_r
            inv_r = np.empty(n, dtype=np.int64)
            inv_r[perm_r] = np.arange(n - 1)
            inv_r = inv_r[keep]
            keep_idx = np.flatnonzero(keep)
            perm_in_keep = np.searchsorted(keep_idx, perm_r)
            self._heat_solve = lambda b: heat_lu.solve(b[perm])[inv]
            self._Lr_solve = lambda b: pois_lu.solve(b[perm_in_keep])[inv_r]
        self.factor_nnz = int(heat_lu.L.nnz + heat_lu.U.nnz
                              + pois_lu.L.nnz + pois_lu.U.nnz)

    def _poisson(self, b: np.ndarray) -> np.ndarray:
        phi = np.zeros(self.mesh.n_vertices)
        phi[self._keep] = self._Lr_solve(b[self._keep])
        return phi

    def distance(self, source: int | Sequence[int]) -> np.ndarray:
        """Return geodesic distance from ``source`` (a vertex index or list)."""
        n = self.mesh.n_vertices
        u0 = np.zeros(n)
        u0[np.atleast_1d(source)] = 1.0

        # 1. short-time heat flow
        u = self._heat_solve(u0)

        # 2. normalized reverse gradient
        g = face_gradient(self.mesh, u, self.areas, self.normals)
        gnorm = np.linalg.norm(g, axis=1, keepdims=True)
        gnorm[gnorm == 0] = 1.0
        X = -g / gnorm

        # 3. Poisson recovery.  The paper solves the Laplace *operator* equation
        # Delta phi = div X; with our positive-semidefinite stiffness L we have
        # Delta = -L, so the system to solve is  L phi = -div X.
        b = integrated_divergence(self.mesh, X, self.cot)
        phi = self._poisson(-b)

        # shift so the source sits at distance 0
        phi -= phi[np.atleast_1d(source)].min()
        return phi


def heat_geodesic(
    mesh: TriMesh, source: int | Sequence[int], m: float = 1.0,
    ordering: str = "colamd",
) -> np.ndarray:
    """Convenience one-shot wrapper around :class:`HeatMethodSolver`."""
    return HeatMethodSolver(mesh, m=m, ordering=ordering).distance(source)
