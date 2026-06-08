"""Static visualization helpers (headless, matplotlib Agg backend).

The project suggests ``pyvista``/``polyscope`` for interactive rendering, but to
keep the benchmark suite reproducible and runnable on a headless machine we
render PNG figures with matplotlib.  Scalar fields are drawn on the sphere with
a 3-D triangulated surface coloured per face, plus an optional orthographic
"map" view (longitude/latitude) that makes contour structure easy to read.
"""

from __future__ import annotations

from typing import Optional, Sequence

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

from .mesh import TriMesh


def _face_values(mesh: TriMesh, vertex_values: np.ndarray) -> np.ndarray:
    return vertex_values[mesh.faces].mean(axis=1)


def plot_scalar_on_sphere(
    mesh: TriMesh,
    values: np.ndarray,
    title: str = "",
    out_path: Optional[str] = None,
    cmap: str = "viridis",
    contours: Optional[int] = None,
    elev: float = 25,
    azim: float = -60,
):
    """Render a per-vertex scalar field on the sphere as a coloured 3-D surface.

    If ``contours`` is given, also draw that many iso-contours on a
    latitude/longitude map projection in a second panel.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    v = mesh.vertices
    f = mesh.faces
    fv = _face_values(mesh, values)

    norm = plt.Normalize(fv.min(), fv.max())
    colors = plt.get_cmap(cmap)(norm(fv))

    has_map = contours is not None
    fig = plt.figure(figsize=(12, 5) if has_map else (6, 5.5))

    ax = fig.add_subplot(1, 2, 1, projection="3d") if has_map else fig.add_subplot(111, projection="3d")
    tris = v[f]
    coll = Poly3DCollection(tris, facecolors=colors, edgecolors="none", linewidths=0)
    ax.add_collection3d(coll)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(title)
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(fv)
    fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.02)

    if has_map:
        ax2 = fig.add_subplot(1, 2, 2)
        lon = np.arctan2(v[:, 1], v[:, 0])
        lat = np.arcsin(np.clip(v[:, 2], -1, 1))
        tri = Triangulation(lon, lat, triangles=f)
        # mask triangles that wrap around the +/-pi longitude seam
        wrap = (np.abs(lon[f[:, 0]] - lon[f[:, 1]]) > np.pi) | \
               (np.abs(lon[f[:, 1]] - lon[f[:, 2]]) > np.pi) | \
               (np.abs(lon[f[:, 2]] - lon[f[:, 0]]) > np.pi)
        tri.set_mask(wrap)
        tcf = ax2.tricontourf(tri, values, levels=contours, cmap=cmap)
        ax2.tricontour(tri, values, levels=contours, colors="k", linewidths=0.4)
        ax2.set_xlabel("longitude"); ax2.set_ylabel("latitude")
        ax2.set_title("map projection + iso-contours")
        fig.colorbar(tcf, ax=ax2, shrink=0.6)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return out_path


def plot_convergence(
    h: Sequence[float],
    series: dict,
    title: str = "Convergence",
    xlabel: str = "mean edge length h",
    out_path: Optional[str] = None,
    reference_slopes: Sequence[int] = (1, 2),
):
    """Log-log convergence plot of one or more error series vs mesh size ``h``."""
    h = np.asarray(h, dtype=float)
    fig, ax = plt.subplots(figsize=(6, 5))
    for label, err in series.items():
        ax.loglog(h, err, "o-", label=label)
    # reference triangles for slope 1 and 2
    h0 = h.max()
    e0 = max(np.max(list(series.values())[0]), 1e-12)
    for p in reference_slopes:
        ax.loglog(h, e0 * (h / h0) ** p, "--", color="gray", alpha=0.6,
                  label=f"O(h^{p})")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("error")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
    return out_path


def plot_eigenfunctions(
    mesh: TriMesh,
    vecs: np.ndarray,
    n: int = 20,
    out_path: Optional[str] = None,
    elev: float = 25,
    azim: float = -60,
):
    """Grid of the first ``n`` eigenfunctions rendered on the sphere."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    cols = 5
    rows = int(np.ceil(n / cols))
    fig = plt.figure(figsize=(3 * cols, 3 * rows))
    v = mesh.vertices
    f = mesh.faces
    tris = v[f]
    for i in range(n):
        ax = fig.add_subplot(rows, cols, i + 1, projection="3d")
        fv = vecs[:, i][f].mean(axis=1)
        norm = plt.Normalize(-np.abs(fv).max(), np.abs(fv).max())
        colors = plt.get_cmap("coolwarm")(norm(fv))
        coll = Poly3DCollection(tris, facecolors=colors, edgecolors="none")
        ax.add_collection3d(coll)
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(f"phi_{i}", fontsize=9)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    return out_path


def plot_heat_series(
    mesh: TriMesh,
    snapshots: Sequence[np.ndarray],
    times: Sequence[float],
    out_path: Optional[str] = None,
    elev: float = 25,
    azim: float = -60,
):
    """Render a row of heat-field snapshots over time."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    n = len(snapshots)
    fig = plt.figure(figsize=(3 * n, 3.2))
    v = mesh.vertices
    f = mesh.faces
    tris = v[f]
    for i, (u, t) in enumerate(zip(snapshots, times)):
        ax = fig.add_subplot(1, n, i + 1, projection="3d")
        fv = u[f].mean(axis=1)
        norm = plt.Normalize(fv.min(), fv.max())
        colors = plt.get_cmap("inferno")(norm(fv))
        coll = Poly3DCollection(tris, facecolors=colors, edgecolors="none")
        ax.add_collection3d(coll)
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=elev, azim=azim)
        ax.set_axis_off()
        ax.set_title(f"t = {t:.3g}", fontsize=9)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    return out_path
