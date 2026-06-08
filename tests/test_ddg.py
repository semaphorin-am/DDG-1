"""Validation tests for the DDG benchmark.

Run with:  python -m pytest tests/  (or `python tests/test_ddg.py`).

These check the mathematical invariants that make the pipeline trustworthy:
mesh topology (Euler characteristic), the cochain identity d^2 = 0, properties
of the cotangent Laplacian, exactness of the per-face gradient on linear
fields, heat conservation, and heat-method accuracy against arccos.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ddg import (
    icosphere,
    build_d0,
    build_d1,
    verify_d_squared,
    face_areas,
    mass_matrix,
    cotangent_laplacian,
    face_gradient,
    HeatSolver,
    point_source,
    total_heat,
    HeatMethodSolver,
    exact_sphere_distance,
    error_metrics,
)


def test_mesh_counts_and_euler():
    for n in range(0, 5):
        mesh = icosphere(n)
        assert mesh.n_vertices == 10 * 4 ** n + 2
        assert mesh.n_faces == 20 * 4 ** n
        # V - E + F = 2 for a sphere
        assert mesh.euler_characteristic() == 2
        # Euler also forces E = 3F/2 for a closed triangulation
        assert mesh.n_edges == 3 * mesh.n_faces // 2


def test_vertices_on_sphere():
    mesh = icosphere(4)
    r = np.linalg.norm(mesh.vertices, axis=1)
    assert np.allclose(r, 1.0, atol=1e-12)


def test_d_squared_zero():
    for n in range(0, 4):
        mesh = icosphere(n)
        d0 = build_d0(mesh)
        d1 = build_d1(mesh)
        assert d0.shape == (mesh.n_edges, mesh.n_vertices)
        assert d1.shape == (mesh.n_faces, mesh.n_edges)
        err = verify_d_squared(d0, d1)
        assert err == 0.0


def test_total_area_converges_to_sphere():
    # Surface area of unit sphere is 4*pi; icosphere underestimates but ->.
    areas = [face_areas(icosphere(n)).sum() for n in range(2, 6)]
    err = [abs(a - 4 * np.pi) for a in areas]
    # monotonically improving and close at the finest level
    assert all(err[i + 1] < err[i] for i in range(len(err) - 1))
    assert err[-1] < 1e-2


def test_laplacian_properties():
    mesh = icosphere(3)
    L, _ = cotangent_laplacian(mesh)
    # symmetry
    assert (abs(L - L.T) > 1e-10).nnz == 0
    # constant in the kernel: L @ 1 = 0
    ones = np.ones(mesh.n_vertices)
    assert np.allclose(L @ ones, 0.0, atol=1e-10)
    # positive semidefinite: x^T L x >= 0 for random x
    rng = np.random.default_rng(0)
    for _ in range(5):
        x = rng.standard_normal(mesh.n_vertices)
        assert x @ (L @ x) >= -1e-9
    # sparse: far fewer than dense entries
    assert L.nnz < mesh.n_vertices ** 2 / 10


def test_mass_matrix_total_area():
    mesh = icosphere(3)
    M = mass_matrix(mesh)
    assert np.isclose(M.diagonal().sum(), face_areas(mesh).sum())


def test_face_gradient_linear_field():
    # gradient of an ambient linear field a.x should be its in-plane projection.
    mesh = icosphere(3)
    from ddg.geometry import face_normals

    a = np.array([0.3, -0.7, 1.1])
    u = mesh.vertices @ a
    g = face_gradient(mesh, u)
    N = face_normals(mesh)
    proj = a[None, :] - (N @ a)[:, None] * N
    assert np.allclose(g, proj, atol=1e-9)


def test_heat_conserves_total():
    mesh = icosphere(3)
    L, _ = cotangent_laplacian(mesh)
    M = mass_matrix(mesh)
    h = mesh.edge_lengths().mean()
    solver = HeatSolver(M, L, dt=h * h)
    u = point_source(mesh.n_vertices, 0)
    q0 = total_heat(M, u)
    for _ in range(10):
        u = solver.step(u)
    # closed surface: integral of u is conserved by implicit Euler too
    assert np.isclose(total_heat(M, u), q0, rtol=1e-9)


def test_heat_method_accuracy():
    mesh = icosphere(5)  # ~10k vertices
    solver = HeatMethodSolver(mesh, m=1.0)
    source = 0
    phi = solver.distance(source)
    exact = exact_sphere_distance(mesh.vertices, mesh.vertices[source])
    m = error_metrics(phi, exact)
    # heat method on a 10k icosphere should be accurate to a few percent
    assert m["MAE"] < 0.02
    assert m["Linf"] < 0.1
    assert phi[source] == 0.0


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"FAIL  {fn.__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)
