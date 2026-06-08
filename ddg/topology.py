"""Discrete exterior derivatives / incidence matrices.

This module answers project question (2) on the *topological* side: the
boundary structure of the simplicial complex, encoded as sparse signed
incidence matrices.

Orientation conventions
------------------------
* Every undirected edge ``(i, j)`` stored in ``mesh.edges`` is oriented from
  its lower- to its higher-indexed endpoint: ``i -> j`` with ``i < j``.
* Every face ``(a, b, c)`` is oriented by its stored vertex winding, giving the
  three directed boundary edges ``a->b``, ``b->c``, ``c->a``.

``d0`` (the vertex->edge incidence) is the discrete exterior derivative on
0-cochains (functions on vertices); ``d1`` (edge->face incidence) is the
exterior derivative on 1-cochains.  The defining identity of any cochain
complex, ``d1 @ d0 == 0`` (discrete ``d^2 = 0``), is the key correctness check
and is exposed via :func:`verify_d_squared`.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .mesh import TriMesh


def build_d0(mesh: TriMesh) -> sp.csr_matrix:
    r"""Vertex-to-edge incidence ``d_0`` with shape ``(|E|, |V|)``.

    Row ``e`` for the oriented edge ``i -> j`` (``i < j``) has ``-1`` in column
    ``i`` and ``+1`` in column ``j``.  Applied to a function ``u`` on vertices,
    ``(d0 @ u)[e] = u[j] - u[i]`` -- the discrete differential of ``u``.
    """
    E = mesh.n_edges
    i = mesh.edges[:, 0]
    j = mesh.edges[:, 1]
    rows = np.concatenate([np.arange(E), np.arange(E)])
    cols = np.concatenate([i, j])
    data = np.concatenate([-np.ones(E), np.ones(E)])
    return sp.csr_matrix((data, (rows, cols)), shape=(E, mesh.n_vertices))


def build_d1(mesh: TriMesh) -> sp.csr_matrix:
    r"""Edge-to-face incidence ``d_1`` with shape ``(|F|, |E|)``.

    For each face, the three directed boundary edges contribute ``+1`` if their
    direction matches the stored edge orientation (``i -> j``, ``i < j``) and
    ``-1`` otherwise.
    """
    F = mesh.n_faces
    idx = mesh.edge_index
    rows = []
    cols = []
    data = []
    for f, (a, b, c) in enumerate(mesh.faces):
        for (u, v) in ((a, b), (b, c), (c, a)):
            u, v = int(u), int(v)
            if u < v:
                e = idx[(u, v)]
                s = 1.0
            else:
                e = idx[(v, u)]
                s = -1.0
            rows.append(f)
            cols.append(e)
            data.append(s)
    return sp.csr_matrix(
        (data, (rows, cols)), shape=(F, mesh.n_edges)
    )


def verify_d_squared(
    d0: sp.spmatrix, d1: sp.spmatrix, tol: float = 1e-12
) -> float:
    r"""Return ``max |d1 @ d0|``; should be ~0 (discrete ``d^2 = 0``).

    Raises ``AssertionError`` if the maximum entry exceeds ``tol``.
    """
    prod = (d1 @ d0).tocoo()
    err = float(np.abs(prod.data).max()) if prod.nnz else 0.0
    assert err <= tol, f"d1 @ d0 != 0 (max |entry| = {err})"
    return err
