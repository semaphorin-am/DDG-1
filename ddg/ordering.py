"""Fill-reducing vertex ordering for mesh Laplacians.

SciPy's SuperLU defaults to COLAMD, a column ordering designed for unsymmetric
matrices.  The matrices in this package are symmetric and come from surface
meshes, where *nested dissection* is the classical low-fill ordering: split
the vertices in two with a thin separator, order both halves first
(recursively) and the separator last.  Eliminating the halves never creates
fill between them, so the factors stay sparse.

The split used here is the simplest one that works on any mesh with vertex
coordinates: cut at the median along the longest axis of the bounding box, and
take as separator the vertices on one side that touch the other side.  The
separator is found from the edge graph, so periodic seams need no special
handling (a seam simply contributes more separator vertices).
"""

from __future__ import annotations

import sys

import numpy as np


def nested_dissection(coords: np.ndarray, edges: np.ndarray, leaf_size: int = 64) -> np.ndarray:
    """Return a permutation of ``range(n)``: the elimination order.

    Parameters
    ----------
    coords : (n, d) vertex coordinates (only used to choose the cuts).
    edges  : (m, 2) undirected edges.
    leaf_size : parts with at most this many vertices are not split further.
    """
    n = coords.shape[0]
    side = np.zeros(n, dtype=np.int8)

    def rec(idx: np.ndarray, e: np.ndarray):
        if idx.size <= leaf_size:
            return [idx]
        c = coords[idx]
        axis = int(np.argmax(c.max(axis=0) - c.min(axis=0)))
        x = c[:, axis]
        half = idx.size // 2
        low = x < np.partition(x, half)[half]
        if not low.any():
            return [idx]
        side[idx] = np.where(low, 1, 2)
        sa, sb = side[e[:, 0]], side[e[:, 1]]
        cut = sa != sb
        sep = np.unique(np.where(sa[cut] == 1, e[cut, 0], e[cut, 1]))
        side[sep] = 3
        sa, sb = side[e[:, 0]], side[e[:, 1]]
        s = side[idx]
        return (rec(idx[s == 1], e[(sa == 1) & (sb == 1)])
                + rec(idx[s == 2], e[(sa == 2) & (sb == 2)])
                + [sep])

    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, 10000))
    try:
        parts = rec(np.arange(n), np.asarray(edges))
    finally:
        sys.setrecursionlimit(old)
    return np.concatenate(parts)
