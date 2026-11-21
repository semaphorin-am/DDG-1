"""Triangle-quality diagnostics used to explain where the solvers degrade."""

from __future__ import annotations

from typing import Dict

import numpy as np

from .mesh import TriMesh
from .geometry import interior_angles, cotangent_laplacian, face_normals


def mesh_quality(mesh: TriMesh) -> Dict[str, float]:
    """Summary of triangle shape for one mesh.

    * ``min_angle_*``: statistics of each triangle's smallest angle (degrees);
      60 is equilateral, small values are slivers.
    * ``frac_obtuse``: share of triangles with an angle above 90 degrees.
    * ``frac_negative_weights``: share of edges whose cotangent weight is
      negative (the edge is not Delaunay).  Negative weights break the discrete
      maximum principle that the heat step relies on.
    * ``n_inverted``: faces whose orientation flipped (sphere: normal points
      inward; flat mesh: negative signed area).
    """
    ang = np.degrees(interior_angles(mesh))
    mn = ang.min(axis=1)
    _, W = cotangent_laplacian(mesh)
    w = np.asarray(W[mesh.edges[:, 0], mesh.edges[:, 1]]).ravel()
    normals = face_normals(mesh)
    if mesh.period is not None:
        inverted = int((normals[:, 2] < 0).sum())
    else:
        centroid = mesh.corner_positions().mean(axis=1)
        inverted = int((np.einsum("ij,ij->i", normals, centroid) < 0).sum())
    return {
        "min_angle_min": float(mn.min()),
        "min_angle_p01": float(np.percentile(mn, 1)),
        "min_angle_mean": float(mn.mean()),
        "max_angle_max": float(ang.max()),
        "frac_obtuse": float((ang.max(axis=1) > 90.0).mean()),
        "frac_negative_weights": float((w < -1e-12).mean()),
        "n_inverted": inverted,
    }


def worst_triangle_vertex(mesh: TriMesh) -> int:
    """Vertex sitting at the smallest angle of the worst-shaped triangle."""
    ang = interior_angles(mesh)
    f = int(np.argmin(ang.min(axis=1)))
    return int(mesh.faces[f, int(np.argmin(ang[f]))])
