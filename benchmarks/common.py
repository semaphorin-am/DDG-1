"""Shared case construction and scoring for the beyond-the-sphere benchmark.

Everything that decides *pass/fail* lives here and reads its numbers from
``tolerances.json``; the tests and the benchmark runner call the same
functions so they cannot drift apart.
"""

from __future__ import annotations

import json
import os
from functools import lru_cache
from typing import Dict, List, Tuple

import numpy as np

from ddg import (TriMesh, icosphere, perturbed_icosphere, flat_torus,
                 cotangent_laplacian, mass_matrix, face_areas, HeatSolver,
                 eigenpairs)
from ddg.quality import worst_triangle_vertex
from . import references as R

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results", "beyond_sphere")
FIGURES = os.path.join(ROOT, "figures", "beyond_sphere")
REPORT = os.path.join(ROOT, "docs", "beyond_sphere_report.html")

with open(os.path.join(HERE, "fixtures.json")) as fh:
    FIX = json.load(fh)
with open(os.path.join(HERE, "tolerances.json")) as fh:
    TOL = json.load(fh)


# --------------------------------------------------------------------------- meshes
@lru_cache(maxsize=None)
def sphere_mesh(level: int, variant: str, seed: int) -> TriMesh:
    amp = FIX["sphere"]["variants"][variant]["amplitude"]
    return perturbed_icosphere(level, amp, seed)


@lru_cache(maxsize=None)
def torus_mesh(nx: int, ny: int, variant: str, seed: int) -> TriMesh:
    t = FIX["torus"]
    return flat_torus(nx, ny, t["Lx"], t["Ly"], t["variants"][variant]["jitter"], seed)


def sphere_cases(levels=None, include_extra=True) -> List[dict]:
    s = FIX["sphere"]
    out = []
    for level in (levels or s["levels"]):
        for variant, spec in s["variants"].items():
            seeds = [0] if spec["amplitude"] == 0 else s["seeds"]
            for seed in seeds:
                out.append(dict(family="sphere", level=level, variant=variant, seed=seed,
                                required=spec["required"], label=f"L{level}"))
    if include_extra and levels is None:
        for level in s.get("extra_levels_regular", []):
            out.append(dict(family="sphere", level=level, variant="regular", seed=0,
                            required=True, label=f"L{level}"))
    return out


def torus_cases(resolutions=None, include_extra=True) -> List[dict]:
    t = FIX["torus"]
    out = []
    for nx, ny in (resolutions or t["resolutions"]):
        for variant, spec in t["variants"].items():
            seeds = [0] if spec["jitter"] == 0 else t["seeds"]
            for seed in seeds:
                out.append(dict(family="torus", nx=nx, ny=ny, variant=variant, seed=seed,
                                required=spec["required"], label=f"{nx}x{ny}"))
    if include_extra and resolutions is None:
        for nx, ny in t.get("extra_resolutions_regular", []):
            out.append(dict(family="torus", nx=nx, ny=ny, variant="regular", seed=0,
                            required=True, label=f"{nx}x{ny}"))
    return out


def build(case: dict) -> TriMesh:
    if case["family"] == "sphere":
        return sphere_mesh(case["level"], case["variant"], case["seed"])
    return torus_mesh(case["nx"], case["ny"], case["variant"], case["seed"])


def diameter(case: dict) -> float:
    if case["family"] == "sphere":
        return float(np.pi)
    return R.torus_diameter(FIX["torus"]["Lx"], FIX["torus"]["Ly"])


def exact_distance(case: dict, mesh: TriMesh, src: int) -> np.ndarray:
    if case["family"] == "sphere":
        return R.sphere_distance(mesh.vertices, mesh.vertices[src])
    t = FIX["torus"]
    return R.torus_distance(mesh.vertices, mesh.vertices[src], t["Lx"], t["Ly"])


def source_vertices(case: dict, mesh: TriMesh) -> Dict[str, int]:
    """Resolve the fixture's source locations to vertex indices on this mesh."""
    out = {}
    if case["family"] == "sphere":
        for name, spec in FIX["sphere"]["sources"].items():
            if spec["kind"] == "vertex":
                out[name] = int(spec["index"])
            elif spec["kind"] == "direction":
                out[name] = int(np.argmin(R.sphere_distance(mesh.vertices, np.array(spec["xyz"]))))
            else:
                out[name] = worst_triangle_vertex(mesh)
    else:
        t = FIX["torus"]
        for name, spec in t["sources"].items():
            p = np.array([spec["frac"][0] * t["Lx"], spec["frac"][1] * t["Ly"]])
            out[name] = int(np.argmin(R.torus_distance(mesh.vertices, p, t["Lx"], t["Ly"])))
    return out


# --------------------------------------------------------------------------- scoring
def score_geodesic(numeric: np.ndarray, exact: np.ndarray, h: float, D: float) -> dict:
    g = TOL["geodesic"]
    err = numeric - exact
    i = int(np.argmax(np.abs(err)))
    mae, linf = float(np.mean(np.abs(err))), float(np.abs(err[i]))
    tol_mae = g["mae_rel_max"]["C"] * (h / D) ** g["mae_rel_max"]["p"]
    tol_linf = g["linf_rel_max"]["C"] * (h / D) ** g["linf_rel_max"]["p"]
    return dict(MAE=mae, Linf=linf, RMSE=float(np.sqrt(np.mean(err ** 2))),
                mean_signed=float(err.mean()),
                mae_rel=mae / D, linf_rel=linf / D,
                tol_mae_rel=tol_mae, tol_linf_rel=tol_linf,
                pass_mae=bool(mae / D <= tol_mae), pass_linf=bool(linf / D <= tol_linf),
                passed=bool(mae / D <= tol_mae and linf / D <= tol_linf),
                worst_vertex=i, exact_at_worst=float(exact[i]))


def score_mode(amplitude: float, residual: float, lam: float, T: float, h: float) -> dict:
    tol = lam * h * h
    if lam == 0:
        return dict(lam_eff=0.0, rate_rel_err=0.0, tol=0.0, residual=residual, passed=True,
                    vacuous=False)
    lam_eff = -np.log(amplitude) / T if amplitude > 0 else float("inf")
    rel = abs(lam_eff / lam - 1.0)
    return dict(lam_eff=float(lam_eff), rate_rel_err=float(rel), tol=float(tol),
                residual=float(residual),
                passed=bool(rel <= tol and residual <= tol),
                vacuous=bool(tol >= 0.5))


def run_mode(mesh: TriMesh, u0: np.ndarray, T: float, L=None, M=None) -> Tuple[float, float, float, int]:
    """Implicit Euler to time ``T`` with dt ~ h^2.  Returns the amplitude of
    ``u0`` left in the solution, the M-norm of everything else (relative), the
    relative drift of total heat, and the number of steps."""
    if L is None:
        L, _ = cotangent_laplacian(mesh)
        M = mass_matrix(mesh)
    h = mesh.edge_lengths().mean()
    n = max(1, int(np.ceil(T / (h * h))))
    solver = HeatSolver(M, L, dt=T / n)
    u = u0.copy()
    for _ in range(n):
        u = solver.step(u)
    w = M.diagonal()
    a = float((w * u * u0).sum() / (w * u0 * u0).sum())
    res = float(np.sqrt((w * (u - a * u0) ** 2).sum() / (w * u0 * u0).sum()))
    # Drift of the integral, relative to int |u0| dA.  Dividing by the integral
    # itself is ill-conditioned for oscillating modes (integral ~ 0).
    drift = abs(float((w * u).sum()) - float((w * u0).sum())) / float((w * np.abs(u0)).sum())
    return a, res, drift, n


def structure_checks(case: dict, mesh: TriMesh, eig: bool = True) -> dict:
    """Connectivity, area and Laplacian invariants for one mesh."""
    s = TOL["structure"]
    f = mesh.faces
    n = mesh.n_vertices
    half = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    und = np.sort(half, axis=1)
    _, counts = np.unique(und[:, 0] * n + und[:, 1], return_counts=True)
    directed_unique = np.unique(half[:, 0] * n + half[:, 1]).size == half.shape[0]
    valence = np.bincount(mesh.edges.ravel(), minlength=n)
    area = float(face_areas(mesh).sum())
    L, _ = cotangent_laplacian(mesh)
    M = mass_matrix(mesh)
    h = float(mesh.edge_lengths().mean())
    out = dict(n_vertices=n, n_edges=mesh.n_edges, n_faces=mesh.n_faces, h_mean=h,
               euler=mesh.euler_characteristic(),
               edge_face_min=int(counts.min()), edge_face_max=int(counts.max()),
               consistently_oriented=bool(directed_unique),
               valence_min=int(valence.min()), valence_max=int(valence.max()),
               area=area,
               laplacian_constant_residual=float(np.abs(L @ np.ones(n)).max()),
               laplacian_symmetry=float(np.abs((L - L.T)).max()))
    if case["family"] == "torus":
        t = FIX["torus"]
        exact_area = t["Lx"] * t["Ly"]
        out["area_rel_err"] = abs(area - exact_area) / exact_area
        out["area_tol"] = s["torus_area_rel_err_max"]
        out["euler_expected"] = 0
        lam1 = R.torus_mode_rate(1, 0, t["Lx"], t["Ly"])
    else:
        out["area_rel_err"] = abs(area - 4 * np.pi) / (4 * np.pi)
        out["area_tol"] = h * h
        out["euler_expected"] = 2
        lam1 = 2.0
    out["lam1_exact"] = lam1
    ok = (out["euler"] == out["euler_expected"] and out["edge_face_min"] == 2
          and out["edge_face_max"] == 2 and out["consistently_oriented"]
          and out["area_rel_err"] <= out["area_tol"]
          and out["laplacian_constant_residual"] <= s["laplacian_constant_residual_max"]
          and out["laplacian_symmetry"] <= s["laplacian_symmetry_max"])
    if eig:
        vals, _ = eigenpairs(L, M, k=3)
        out["eig0"], out["eig1"] = float(vals[0]), float(vals[1])
        out["n_null"] = int((np.abs(vals) < 1e-8 * lam1).sum())
        out["eig1_rel_err"] = abs(vals[1] - lam1) / lam1
        ok = ok and out["n_null"] == 1
    out["passed"] = bool(ok)
    return out
