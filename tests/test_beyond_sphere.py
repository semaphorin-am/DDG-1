"""Checks for the perturbed-sphere and flat-torus benchmark.

Every pass/fail decision goes through ``benchmarks.common`` and therefore
through ``benchmarks/tolerances.json`` -- the same code path as the benchmark.
Known failures are recorded as strict ``xfail`` with the reason, not hidden
and not loosened; if one starts passing, the suite fails until it is reviewed.

Run:  python -m pytest -q tests
"""

import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from ddg import (HeatMethodSolver, cotangent_laplacian, mass_matrix, face_areas, icosphere,
                 perturbed_icosphere, flat_torus, open_rectangle, nested_dissection)
from ddg.quality import mesh_quality
from benchmarks import references as R
from benchmarks.common import (FIX, TOL, RESULTS, build, diameter, exact_distance, source_vertices,
                               score_geodesic, score_mode, run_mode, structure_checks,
                               sphere_cases, torus_cases)

T = FIX["torus"]
LX, LY = T["Lx"], T["Ly"]
SMALL_TORUS = [r for r in T["resolutions"] if r[0] <= 96]
SMALL_SPHERE = [l for l in FIX["sphere"]["levels"] if l <= 4]


def _id(c):
    return f"{c['family']}-{c['label']}-{c['variant']}-{c['seed']}"


TORUS = torus_cases(resolutions=SMALL_TORUS)
SPHERE = [c for c in sphere_cases(levels=SMALL_SPHERE)]


# ------------------------------------------------------------------ references
def test_reference_torus_distance_matches_image_sum():
    rng = np.random.default_rng(0)
    p = rng.random((500, 2)) * [LX, LY]
    for s in ([0.0, 0.0], [LX * 0.97, LY * 0.37], [LX * 0.98, LY * 0.97]):
        a = R.torus_distance(p, np.array(s), LX, LY)
        b = R.torus_distance_images(p, np.array(s), LX, LY)
        assert np.abs(a - b).max() < 1e-14
        assert a.max() <= R.torus_diameter(LX, LY) + 1e-14


def test_reference_sphere_distance_matches_arccos():
    rng = np.random.default_rng(1)
    p = rng.normal(size=(500, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    s = p[0]
    assert np.abs(R.sphere_distance(p, s) - np.arccos(np.clip(p @ s, -1, 1))).max() < 1e-7


# ------------------------------------------------------------------- structure
@pytest.mark.parametrize("case", TORUS + SPHERE, ids=_id)
def test_structure(case):
    out = structure_checks(case, build(case))
    assert out["passed"], out


def test_torus_seam_connectivity():
    nx, ny = 12, 8
    mesh = flat_torus(nx, ny, LX, LY)
    assert mesh.n_vertices == nx * ny            # no duplicated seam vertices
    assert mesh.n_faces == 2 * nx * ny
    assert mesh.euler_characteristic() == 0
    edges = {tuple(e) for e in np.sort(mesh.edges, axis=1).tolist()}
    for j in range(ny):                           # x-seam: column nx-1 joins column 0
        assert tuple(sorted(((nx - 1) + nx * j, 0 + nx * j))) in edges or \
               tuple(sorted(((nx - 1) * ny + j, j))) in edges
    assert np.bincount(mesh.edges.ravel()).min() == 6 == np.bincount(mesh.edges.ravel()).max()
    assert abs(face_areas(mesh).sum() - LX * LY) < 1e-12
    # every edge, including seam-crossing ones, is short
    assert mesh.edge_lengths().max() < 1.01 * np.hypot(LX / nx, LY / ny)


def test_torus_translation_invariance():
    """On the regular torus every vertex is equivalent, so the error field
    must not depend on whether the source sits on a seam."""
    case = dict(family="torus", nx=24, ny=16, variant="regular", seed=0)
    mesh = build(case)
    solver = HeatMethodSolver(mesh)
    stats = []
    for src in source_vertices(case, mesh).values():
        err = solver.distance(src) - exact_distance(case, mesh, src)
        stats.append((np.abs(err).mean(), np.abs(err).max()))
    stats = np.array(stats)
    assert np.ptp(stats, axis=0).max() < 1e-10


# ------------------------------------------------------------------- geodesics
def _geo(case, m=None):
    mesh = build(case)
    h, D = float(mesh.edge_lengths().mean()), diameter(case)
    solver = HeatMethodSolver(mesh, m=FIX["heat_method_m"] if m is None else m)
    return {name: score_geodesic(solver.distance(s), exact_distance(case, mesh, s), h, D)
            for name, s in source_vertices(case, mesh).items()}


@pytest.mark.parametrize("case", [c for c in SPHERE if c["required"]], ids=_id)
def test_sphere_geodesic_within_tolerance(case):
    for name, sc in _geo(case).items():
        assert sc["passed"], (name, sc)


def test_sphere_baseline_preserved():
    with open(os.path.join(RESULTS, "baseline.json")) as fh:
        base = {r["subdivisions"]: r for r in json.load(fh)["rows"]}
    for level in SMALL_SPHERE:
        mesh = icosphere(level)
        assert np.array_equal(mesh.vertices, perturbed_icosphere(level, 0.0, 0).vertices)
        phi = HeatMethodSolver(mesh, m=1.0).distance(0)
        mae = np.abs(phi - R.sphere_distance(mesh.vertices, mesh.vertices[0])).mean()
        assert abs(mae - base[level]["MAE"]) < 1e-9


def test_sphere_strong_perturbation_is_worse_quality():
    q0 = mesh_quality(build(dict(family="sphere", level=3, variant="regular", seed=0)))
    q1 = mesh_quality(build(dict(family="sphere", level=3, variant="strong", seed=11)))
    assert q1["min_angle_min"] < q0["min_angle_min"]
    assert q1["frac_obtuse"] > q0["frac_obtuse"]


MAE_OK = [c for c in TORUS if not (c["variant"] == "jittered" and c["nx"] >= 96)]


@pytest.mark.parametrize("case", MAE_OK, ids=_id)
def test_torus_geodesic_mae_within_tolerance(case):
    """Mean error near both seams, the corner and the interior."""
    for name, sc in _geo(case).items():
        assert sc["pass_mae"], (name, sc)


@pytest.mark.xfail(strict=True, reason=(
    "Known failure, not relaxed: on the flat torus the heat method's worst-point error at "
    "t = h^2 exceeds the sphere-derived L-inf tolerance at every resolution, for every source "
    "including the interior one. See docs/beyond_sphere_progress.md section 4."))
@pytest.mark.parametrize("case", [c for c in TORUS if c["variant"] == "regular"], ids=_id)
def test_torus_geodesic_linf_within_tolerance(case):
    for name, sc in _geo(case).items():
        assert sc["pass_linf"], (name, sc)


@pytest.mark.xfail(strict=True, reason=(
    "Known failure, not relaxed: with the default t = h^2 the torus mean error stops "
    "improving near 0.6% of the diameter, so it misses the tolerance from 192x128 on."))
def test_torus_geodesic_jittered_96x64_mae():
    """Same plateau, reached one resolution earlier on jittered meshes
    (2 of 15 seed/source combinations miss the MAE tolerance)."""
    cases = [c for c in TORUS if c["variant"] == "jittered" and c["nx"] == 96]
    assert all(sc["pass_mae"] for c in cases for sc in _geo(c).values())


@pytest.mark.xfail(strict=True, reason=(
    "Known failure, not relaxed: with the default t = h^2 the torus mean error stops "
    "improving near 0.6% of the diameter, so it misses the tolerance from 192x128 on."))
def test_torus_geodesic_fine_grid_mae():
    case = dict(family="torus", nx=192, ny=128, variant="regular", seed=0)
    assert all(sc["pass_mae"] for sc in _geo(case).values())


# ------------------------------------------------------------------------ heat
@pytest.mark.parametrize("case", [c for c in TORUS if c["nx"] >= 48], ids=_id)
def test_torus_fourier_mode_decay(case):
    mesh = build(case)
    L, _ = cotangent_laplacian(mesh); M = mass_matrix(mesh)
    h = float(mesh.edge_lengths().mean())
    Tm = T["heat_modes"]["T"]
    for k, l in T["heat_modes"]["modes"]:
        lam = R.torus_mode_rate(k, l, LX, LY)
        a, res, drift, _ = run_mode(mesh, R.torus_mode(mesh.vertices, k, l, LX, LY), Tm, L, M)
        assert drift <= TOL["heat_mode"]["total_heat_rel_drift_max"]
        if (k, l) == (0, 0):
            assert abs(a - 1) + res <= TOL["heat_mode"]["constant_mode_max_dev"]
        else:
            sc = score_mode(a, res, lam, Tm, h)
            assert sc["passed"], ((k, l), sc)


def test_torus_constant_mode_pointwise():
    mesh = flat_torus(48, 32, LX, LY, jitter=0.25, seed=11)
    L, _ = cotangent_laplacian(mesh); M = mass_matrix(mesh)
    assert np.abs(L @ np.ones(mesh.n_vertices)).max() <= TOL["structure"]["laplacian_constant_residual_max"]
    from ddg import HeatSolver
    s = HeatSolver(M, L, dt=1e-3)
    u = np.full(mesh.n_vertices, 3.0)
    for _ in range(50):
        u = s.step(u)
    assert np.abs(u - 3.0).max() <= TOL["heat_mode"]["constant_mode_max_dev"]


def test_sphere_harmonic_decay():
    for case in [c for c in SPHERE if c["required"] and c["level"] >= 3]:
        mesh = build(case)
        h = float(mesh.edge_lengths().mean())
        for l in FIX["sphere"]["heat_modes"]["degrees"]:
            a, res, _, _ = run_mode(mesh, R.sphere_harmonic(mesh.vertices, l), FIX["sphere"]["heat_modes"]["T"])
            sc = score_mode(a, res, R.sphere_mode_rate(l), FIX["sphere"]["heat_modes"]["T"], h)
            assert sc["passed"], (_id(case), l, sc)


# ------------------------------------------------------------- test the tests
SEAM_SOURCES = [s for s in T["sources"] if s != "center"]


@pytest.mark.parametrize("res", SMALL_TORUS, ids=lambda r: f"{r[0]}x{r[1]}")
def test_negative_unwrapped_reference_is_rejected(res):
    """Scoring the straight-line-in-the-rectangle distance must fail."""
    case = dict(family="torus", nx=res[0], ny=res[1], variant="regular", seed=0)
    mesh = build(case)
    h, D = float(mesh.edge_lengths().mean()), diameter(case)
    srcs = source_vertices(case, mesh)
    for name in SEAM_SOURCES:
        p = mesh.vertices[srcs[name]]
        sc = score_geodesic(R.unwrapped_distance(mesh.vertices, p),
                            R.torus_distance(mesh.vertices, p, LX, LY), h, D)
        assert not sc["pass_mae"] and not sc["pass_linf"], (name, sc)
        assert sc["linf_rel"] > 5 * sc["tol_linf_rel"]


@pytest.mark.parametrize("res", SMALL_TORUS, ids=lambda r: f"{r[0]}x{r[1]}")
def test_negative_unstitched_mesh_is_rejected(res):
    """The real solver on a rectangle whose seam was left open must fail."""
    mesh = open_rectangle(res[0], res[1], LX, LY)
    assert mesh.euler_characteristic() == 1
    h, D = float(mesh.edge_lengths().mean()), R.torus_diameter(LX, LY)
    solver = HeatMethodSolver(mesh)
    for name in SEAM_SOURCES:
        f = T["sources"][name]["frac"]
        src = int(np.argmin(np.hypot(mesh.vertices[:, 0] - f[0] * LX, mesh.vertices[:, 1] - f[1] * LY)))
        sc = score_geodesic(solver.distance(src),
                            R.torus_distance(mesh.vertices, mesh.vertices[src], LX, LY), h, D)
        assert not sc["pass_mae"] and not sc["pass_linf"], (name, sc)


def test_negative_wrong_decay_rate_is_rejected():
    h = float(build(dict(family="torus", nx=96, ny=64, variant="regular", seed=0)).edge_lengths().mean())
    Tm = T["heat_modes"]["T"]
    for k, l in T["heat_modes"]["modes"]:
        if (k, l) == (0, 0):
            continue
        lam = R.torus_mode_rate(k, l, LX, LY)
        assert not score_mode(float(np.exp(-2 * lam * Tm)), 0.0, lam, Tm, h)["passed"]
        assert score_mode(float(np.exp(-lam * Tm)), 0.0, lam, Tm, h)["passed"]


def test_negative_constant_mode_leak_is_rejected():
    assert abs(0.999 - 1) > TOL["heat_mode"]["constant_mode_max_dev"]


# ---------------------------------------------------------------- optimisation
@pytest.mark.parametrize("case", [
    dict(family="sphere", level=4, variant="regular", seed=0),
    dict(family="sphere", level=4, variant="strong", seed=47),
    dict(family="torus", nx=48, ny=32, variant="regular", seed=0),
    dict(family="torus", nx=48, ny=32, variant="jittered", seed=23)], ids=["sphere-regular", "sphere-strong", "torus-regular", "torus-jittered"])
def test_nested_ordering_is_equivalent(case):
    mesh = build(case)
    perm = nested_dissection(mesh.vertices, mesh.edges)
    assert np.array_equal(np.sort(perm), np.arange(mesh.n_vertices))
    a = HeatMethodSolver(mesh, ordering="colamd")
    b = HeatMethodSolver(mesh, ordering="nested")
    for src in source_vertices(case, mesh).values():
        assert np.abs(a.distance(src) - b.distance(src)).max() <= 1e-8
