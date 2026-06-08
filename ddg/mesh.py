"""Triangular mesh representation and icosphere generation.

This module answers project question (1): *how are simplicial complexes and
triangular meshes represented computationally?*

A triangle mesh is stored with three arrays:

* ``vertices``  -- ``(|V|, 3)`` float array of vertex coordinates.
* ``faces``     -- ``(|F|, 3)`` int array of vertex indices, one row per
  triangle, ordered counter-clockwise as seen from outside the surface
  (consistent outward orientation).
* ``edges``     -- ``(|E|, 2)`` int array of *unique* undirected edges,
  derived from the faces, with the convention ``edge[k] = (i, j)`` and
  ``i < j``.

The icosphere is built by recursively subdividing a regular icosahedron and
projecting every vertex back onto the unit sphere.  This produces meshes whose
triangles are nearly equilateral and uniform, which is exactly what we want for
a clean convergence study.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Base icosahedron
# ---------------------------------------------------------------------------
def _icosahedron() -> Tuple[np.ndarray, np.ndarray]:
    """Return the 12 vertices / 20 faces of a unit regular icosahedron.

    The faces are consistently oriented so that the vertex winding is
    counter-clockwise when viewed from outside the sphere (outward normals).
    """
    phi = (1.0 + np.sqrt(5.0)) / 2.0  # golden ratio
    verts = np.array(
        [
            (-1,  phi, 0), ( 1,  phi, 0), (-1, -phi, 0), ( 1, -phi, 0),
            ( 0, -1,  phi), ( 0,  1,  phi), ( 0, -1, -phi), ( 0,  1, -phi),
            ( phi, 0, -1), ( phi, 0,  1), (-phi, 0, -1), (-phi, 0,  1),
        ],
        dtype=np.float64,
    )
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)

    faces = np.array(
        [
            (0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
            (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
            (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
            (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1),
        ],
        dtype=np.int64,
    )
    return verts, faces


def _subdivide(
    verts: np.ndarray, faces: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """One step of 1-to-4 loop subdivision, projecting midpoints to the sphere.

    Each triangle ``(a, b, c)`` is split into four by inserting the three edge
    midpoints.  Midpoints are shared between adjacent triangles via a cache
    keyed on the (sorted) endpoint pair, so the result stays watertight.
    """
    verts_list = list(map(tuple, verts))
    midpoint_cache: Dict[Tuple[int, int], int] = {}
    new_faces = []

    def midpoint(i: int, j: int) -> int:
        key = (i, j) if i < j else (j, i)
        cached = midpoint_cache.get(key)
        if cached is not None:
            return cached
        m = (verts[i] + verts[j]) / 2.0
        m /= np.linalg.norm(m)  # project onto the unit sphere
        idx = len(verts_list)
        verts_list.append(tuple(m))
        midpoint_cache[key] = idx
        return idx

    for a, b, c in faces:
        ab = midpoint(a, b)
        bc = midpoint(b, c)
        ca = midpoint(c, a)
        new_faces.extend(
            [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        )

    return (
        np.array(verts_list, dtype=np.float64),
        np.array(new_faces, dtype=np.int64),
    )


# ---------------------------------------------------------------------------
# Mesh container
# ---------------------------------------------------------------------------
@dataclass
class TriMesh:
    """A triangulated surface mesh.

    Attributes
    ----------
    vertices : (|V|, 3) float ndarray
    faces    : (|F|, 3) int ndarray, consistently oriented
    edges    : (|E|, 2) int ndarray of unique undirected edges (i < j)
    edge_index : dict mapping (i, j) with i < j -> edge row
    """

    vertices: np.ndarray
    faces: np.ndarray
    edges: np.ndarray = field(default=None)
    edge_index: Dict[Tuple[int, int], int] = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.vertices = np.asarray(self.vertices, dtype=np.float64)
        self.faces = np.asarray(self.faces, dtype=np.int64)
        if self.edges is None:
            self._build_edges()

    # -- derived connectivity ------------------------------------------------
    def _build_edges(self) -> None:
        f = self.faces
        # all three directed edges of every face
        raw = np.vstack([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
        raw_sorted = np.sort(raw, axis=1)
        edges = np.unique(raw_sorted, axis=0)
        self.edges = edges
        self.edge_index = {(int(i), int(j)): k for k, (i, j) in enumerate(edges)}

    # -- sizes ---------------------------------------------------------------
    @property
    def n_vertices(self) -> int:
        return self.vertices.shape[0]

    @property
    def n_faces(self) -> int:
        return self.faces.shape[0]

    @property
    def n_edges(self) -> int:
        return self.edges.shape[0]

    # -- statistics ----------------------------------------------------------
    def edge_lengths(self) -> np.ndarray:
        v = self.vertices
        return np.linalg.norm(v[self.edges[:, 0]] - v[self.edges[:, 1]], axis=1)

    def euler_characteristic(self) -> int:
        r"""V - E + F.  Equal to 2 for any mesh homeomorphic to a sphere."""
        return self.n_vertices - self.n_edges + self.n_faces

    def statistics(self) -> Dict[str, float]:
        lengths = self.edge_lengths()
        return {
            "n_vertices": self.n_vertices,
            "n_edges": self.n_edges,
            "n_faces": self.n_faces,
            "euler_characteristic": self.euler_characteristic(),
            "mean_edge_length": float(lengths.mean()),
            "min_edge_length": float(lengths.min()),
            "max_edge_length": float(lengths.max()),
        }


# ---------------------------------------------------------------------------
# Public constructors
# ---------------------------------------------------------------------------
def icosphere(subdivisions: int = 3) -> TriMesh:
    """Build an icosphere of the given subdivision level.

    The vertex count follows ``|V| = 10 * 4**n + 2`` for ``n`` subdivisions::

        n   |V|       |F|
        0   12        20
        1   42        80
        2   162       320
        3   642       1280
        4   2562      5120
        5   10242     20480
        6   40962     81920
        7   163842    327680
        8   655362    1310720
    """
    if subdivisions < 0:
        raise ValueError("subdivisions must be >= 0")
    verts, faces = _icosahedron()
    for _ in range(subdivisions):
        verts, faces = _subdivide(verts, faces)
    # guard against tiny drift from repeated averaging
    verts /= np.linalg.norm(verts, axis=1, keepdims=True)
    return TriMesh(verts, faces)


# Map a *target* vertex count to the smallest subdivision level that reaches it.
def icosphere_for_target(target_vertices: int) -> TriMesh:
    """Return the icosphere whose vertex count is closest to ``target``."""
    counts = [(n, 10 * 4 ** n + 2) for n in range(0, 9)]
    best = min(counts, key=lambda nc: abs(nc[1] - target_vertices))
    return icosphere(best[0])
