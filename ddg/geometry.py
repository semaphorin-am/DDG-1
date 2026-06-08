"""Discrete geometry: areas, cotangent weights, mass matrix, Laplacian.

This module answers project question (2) on the *metric* side, and underpins
questions (3)-(4).  Given a :class:`~ddg.mesh.TriMesh` it builds the operators
of Discrete Exterior Calculus needed to solve PDEs on the surface:

* per-face areas and angles,
* per-vertex (lumped barycentric) areas -> diagonal **mass matrix** ``M``
  (the discrete Hodge star ``*_0`` on 0-forms),
* **cotangent weights** on edges (the Hodge star ``*_1`` on 1-forms),
* the **cotangent Laplacian** ``L`` -- a symmetric positive-semidefinite
  stiffness matrix with ``L @ 1 = 0``.

Sign convention
---------------
``L`` here is the *positive-semidefinite stiffness* matrix

    ``L_ii = sum_j w_ij``,   ``L_ij = -w_ij``   (``w_ij = (cot a + cot b)/2``).

The Laplace-Beltrami operator is then ``Delta = -M^{-1} L`` (negative
semidefinite), so the heat equation ``u_t = Delta u`` discretizes to the
implicit-Euler step ``(M + t L) u^{n+1} = M u^n``.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
import scipy.sparse as sp

from .mesh import TriMesh


# ---------------------------------------------------------------------------
# Per-face quantities
# ---------------------------------------------------------------------------
def face_areas(mesh: TriMesh) -> np.ndarray:
    """Return the ``(|F|,)`` array of triangle areas."""
    v = mesh.vertices
    f = mesh.faces
    e1 = v[f[:, 1]] - v[f[:, 0]]
    e2 = v[f[:, 2]] - v[f[:, 0]]
    return 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1)


def face_normals(mesh: TriMesh, normalize: bool = True) -> np.ndarray:
    """Return ``(|F|, 3)`` (unit) face normals from the stored winding."""
    v = mesh.vertices
    f = mesh.faces
    n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    if normalize:
        norm = np.linalg.norm(n, axis=1, keepdims=True)
        norm[norm == 0] = 1.0
        n = n / norm
    return n


def interior_angles(mesh: TriMesh) -> np.ndarray:
    """Return ``(|F|, 3)`` interior angles; column ``k`` is the angle at the
    ``k``-th vertex of each face."""
    v = mesh.vertices
    f = mesh.faces
    angles = np.empty((mesh.n_faces, 3))
    for k in range(3):
        a = f[:, k]
        b = f[:, (k + 1) % 3]
        c = f[:, (k + 2) % 3]
        u = v[b] - v[a]
        w = v[c] - v[a]
        u = u / np.linalg.norm(u, axis=1, keepdims=True)
        w = w / np.linalg.norm(w, axis=1, keepdims=True)
        cos = np.clip(np.einsum("ij,ij->i", u, w), -1.0, 1.0)
        angles[:, k] = np.arccos(cos)
    return angles


# ---------------------------------------------------------------------------
# Per-vertex quantities and the mass matrix
# ---------------------------------------------------------------------------
def vertex_areas(mesh: TriMesh, areas: np.ndarray | None = None) -> np.ndarray:
    """Lumped barycentric vertex areas: each vertex receives one third of the
    area of every incident triangle.  Their sum equals the total surface area.
    """
    if areas is None:
        areas = face_areas(mesh)
    out = np.zeros(mesh.n_vertices)
    np.add.at(out, mesh.faces[:, 0], areas / 3.0)
    np.add.at(out, mesh.faces[:, 1], areas / 3.0)
    np.add.at(out, mesh.faces[:, 2], areas / 3.0)
    return out


def mass_matrix(mesh: TriMesh, areas: np.ndarray | None = None) -> sp.dia_matrix:
    """Diagonal (lumped) mass matrix ``M`` = discrete Hodge star on 0-forms."""
    va = vertex_areas(mesh, areas)
    return sp.diags(va, format="dia")


# ---------------------------------------------------------------------------
# Cotangent Laplacian
# ---------------------------------------------------------------------------
def cotangent_laplacian(
    mesh: TriMesh,
) -> Tuple[sp.csr_matrix, sp.csr_matrix]:
    """Build the cotangent Laplacian.

    Returns
    -------
    L : (|V|, |V|) csr_matrix
        Symmetric positive-semidefinite stiffness matrix with ``L @ 1 = 0``.
    C : (|V|, |V|) csr_matrix
        The signed off-diagonal cotangent weight matrix (``C_ij = w_ij`` for an
        edge, 0 otherwise), useful for assembling other operators.

    For each triangle, the angle at a vertex contributes ``cot(angle)/2`` to the
    weight of the *opposite* edge.  Summing the (up to two) triangles incident
    to an interior edge gives ``w_ij = (cot a + cot b) / 2``.
    """
    v = mesh.vertices
    f = mesh.faces
    angles = interior_angles(mesh)
    cot = 1.0 / np.tan(angles)  # (|F|, 3)

    # For face vertex k, the opposite edge connects vertices (k+1, k+2).
    I, J, W = [], [], []
    for k in range(3):
        a = f[:, (k + 1) % 3]
        b = f[:, (k + 2) % 3]
        w = 0.5 * cot[:, k]
        I.append(a); J.append(b); W.append(w)
        I.append(b); J.append(a); W.append(w)  # symmetric

    I = np.concatenate(I)
    J = np.concatenate(J)
    W = np.concatenate(W)

    n = mesh.n_vertices
    # off-diagonal weight matrix W_ij = w_ij (summed over shared triangles)
    Wmat = sp.csr_matrix((W, (I, J)), shape=(n, n))
    # stiffness: L = diag(row sums of W) - W
    diag = np.asarray(Wmat.sum(axis=1)).ravel()
    L = sp.diags(diag) - Wmat
    return L.tocsr(), Wmat.tocsr()


def cotangent_weights_per_edge(mesh: TriMesh) -> np.ndarray:
    """Return the scalar cotangent weight ``w_ij`` for each edge in
    ``mesh.edges`` order (the Hodge star ``*_1`` on 1-forms)."""
    _, Wmat = cotangent_laplacian(mesh)
    return np.asarray(Wmat[mesh.edges[:, 0], mesh.edges[:, 1]]).ravel()
