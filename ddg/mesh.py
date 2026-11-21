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
    # Optional (|F|, 3, 3) corner positions, one copy per face.  When set, all
    # metric quantities (areas, angles, gradients, edge lengths) are computed
    # from these instead of ``vertices[faces]``.  This is what lets a mesh carry
    # an intrinsic metric that no single vertex embedding has, e.g. a flat
    # torus whose seam triangles are "unwrapped" across the identified edge.
    corners: np.ndarray = field(default=None, repr=False)
    # (Lx, Ly) for periodic flat meshes, else None.
    period: Tuple[float, float] = field(default=None)

    def __post_init__(self) -> None:
        self.vertices = np.asarray(self.vertices, dtype=np.float64)
        self.faces = np.asarray(self.faces, dtype=np.int64)
        if self.corners is not None:
            self.corners = np.asarray(self.corners, dtype=np.float64)
        if self.edges is None:
            self._build_edges()

    def corner_positions(self) -> np.ndarray:
        """``(|F|, 3, 3)`` positions of each face's three corners."""
        if self.corners is not None:
            return self.corners
        return self.vertices[self.faces]

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
        if self.corners is not None:
            # Intrinsic metric: read each edge's length off a face that owns it.
            c = self.corners
            f = self.faces
            n = self.n_vertices
            keys = np.concatenate([
                np.minimum(f[:, k], f[:, (k + 1) % 3]) * n
                + np.maximum(f[:, k], f[:, (k + 1) % 3]) for k in range(3)])
            lens = np.concatenate([
                np.linalg.norm(c[:, (k + 1) % 3] - c[:, k], axis=1) for k in range(3)])
            edge_keys = self.edges[:, 0] * n + self.edges[:, 1]
            order = np.argsort(keys)
            pos = np.searchsorted(keys[order], edge_keys)
            return lens[order][pos]
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


def _disk_samples(rng: np.random.Generator, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """``n`` points uniform in the unit disk, as (radius, angle)."""
    return np.sqrt(rng.random(n)), 2 * np.pi * rng.random(n)


def perturbed_icosphere(subdivisions: int, amplitude: float, seed: int = 0) -> TriMesh:
    """Icosphere with every vertex jittered *along the sphere*.

    Each vertex moves by a random tangential offset, uniform in a disk of
    radius ``amplitude * (mean length of its incident edges)``, and is then
    projected back to the unit sphere.  Connectivity is unchanged, and every
    vertex still lies exactly on the sphere, so great-circle distance remains
    the exact reference; only triangle shape degrades.
    """
    base = icosphere(subdivisions)
    if amplitude == 0.0:
        return base
    v = base.vertices
    e = base.edges
    lens = base.edge_lengths()
    n = base.n_vertices
    local = (np.bincount(e[:, 0], lens, n) + np.bincount(e[:, 1], lens, n)) / (
        np.bincount(e[:, 0], minlength=n) + np.bincount(e[:, 1], minlength=n))
    rng = np.random.default_rng([int(seed), int(subdivisions)])
    r, th = _disk_samples(rng, n)
    # orthonormal tangent frame at each vertex
    helper = np.where(np.abs(v[:, [2]]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    t1 = np.cross(v, helper)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(v, t1)
    disp = (amplitude * local * r)[:, None] * (
        np.cos(th)[:, None] * t1 + np.sin(th)[:, None] * t2)
    w = v + disp
    w /= np.linalg.norm(w, axis=1, keepdims=True)
    return TriMesh(w, base.faces.copy())


def _grid_positions(nx: int, ny: int, Lx: float, Ly: float, jitter: float,
                    seed: int) -> np.ndarray:
    hx, hy = Lx / nx, Ly / ny
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    P = np.stack([ii * hx, jj * hy], axis=-1).astype(np.float64)
    if jitter:
        rng = np.random.default_rng([int(seed), int(nx), int(ny)])
        r, th = _disk_samples(rng, nx * ny)
        d = jitter * min(hx, hy) * r
        P += np.stack([d * np.cos(th), d * np.sin(th)], axis=-1).reshape(nx, ny, 2)
    return P


def flat_torus(nx: int, ny: int, Lx: float = 1.0, Ly: float = 1.0,
               jitter: float = 0.0, seed: int = 0) -> TriMesh:
    """Flat periodic torus: the rectangle ``[0, Lx) x [0, Ly)`` with opposite
    edges identified (not an embedded doughnut).

    ``nx * ny`` vertices on a grid, two triangles per cell.  Column ``nx`` *is*
    column 0 and row ``ny`` *is* row 0, so there are no boundary vertices and no
    duplicated seam vertices.  ``vertices`` holds each vertex once, wrapped
    into the rectangle (z = 0); ``corners`` holds the unwrapped positions used
    by every triangle, so triangles that straddle a seam have their true flat
    shape.  ``jitter`` displaces vertices by up to that fraction of the grid
    spacing (same displacement for a vertex in all its triangles).
    """
    if nx < 3 or ny < 3:
        raise ValueError("need nx, ny >= 3 for a simplicial torus")
    P = _grid_positions(nx, ny, Lx, Ly, jitter, seed)
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    ii = ii.ravel(); jj = jj.ravel()

    def vid(i, j):
        return (i % nx) * ny + (j % ny)

    def pos(i, j):  # unwrapped position of grid node (i, j), i <= nx, j <= ny
        p = P[i % nx, j % ny].copy()
        p[:, 0] += (i // nx) * Lx
        p[:, 1] += (j // ny) * Ly
        return p

    quads = [(ii, jj), (ii + 1, jj), (ii + 1, jj + 1), (ii, jj + 1)]
    ids = [vid(i, j) for i, j in quads]
    xy = [pos(i, j) for i, j in quads]
    faces = np.concatenate([np.stack([ids[0], ids[1], ids[2]], 1),
                            np.stack([ids[0], ids[2], ids[3]], 1)])
    c2 = np.concatenate([np.stack([xy[0], xy[1], xy[2]], 1),
                         np.stack([xy[0], xy[2], xy[3]], 1)])
    corners = np.concatenate([c2, np.zeros(c2.shape[:2] + (1,))], axis=2)
    verts = np.zeros((nx * ny, 3))
    verts[:, 0] = P[..., 0].ravel() % Lx
    verts[:, 1] = P[..., 1].ravel() % Ly
    return TriMesh(verts, faces, corners=corners, period=(float(Lx), float(Ly)))


def open_rectangle(nx: int, ny: int, Lx: float = 1.0, Ly: float = 1.0) -> TriMesh:
    """The same grid as :func:`flat_torus` but with the seams left *open*
    (``(nx+1) * (ny+1)`` vertices, a disk with boundary).  Exists only as a
    negative control: distances on it cannot wrap around."""
    ii, jj = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), indexing="ij")
    verts = np.stack([ii.ravel() * Lx / nx, jj.ravel() * Ly / ny,
                      np.zeros(ii.size)], axis=1)
    i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij")
    i = i.ravel(); j = j.ravel()
    a = i * (ny + 1) + j; b = (i + 1) * (ny + 1) + j
    c = (i + 1) * (ny + 1) + j + 1; d = i * (ny + 1) + j + 1
    faces = np.concatenate([np.stack([a, b, c], 1), np.stack([a, c, d], 1)])
    return TriMesh(verts, faces)


# Map a *target* vertex count to the smallest subdivision level that reaches it.
def icosphere_for_target(target_vertices: int) -> TriMesh:
    """Return the icosphere whose vertex count is closest to ``target``."""
    counts = [(n, 10 * 4 ** n + 2) for n in range(0, 9)]
    best = min(counts, key=lambda nc: abs(nc[1] - target_vertices))
    return icosphere(best[0])
