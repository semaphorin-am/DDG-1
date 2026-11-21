#!/usr/bin/env python3
"""Accuracy sweep: geodesics, heat modes, structure, negative controls.

Writes CSV/JSON to ``results/beyond_sphere/`` plus ``fields.npz`` (the
per-vertex error fields the report displays).  Uses the *baseline* solver
configuration throughout; the optimisation is measured separately in
``benchmarks/perf.py``.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import scipy
from ddg import HeatMethodSolver, cotangent_laplacian, mass_matrix, open_rectangle, HeatSolver
from ddg.quality import mesh_quality
from ddg.geometry import interior_angles
from benchmarks import references as R
from benchmarks.common import (FIX, TOL, RESULTS, sphere_cases, torus_cases, build, diameter,
                               exact_distance, source_vertices, score_geodesic, score_mode,
                               run_mode, structure_checks)

# Meshes whose per-vertex fields are stored for the report (kept small).
VIEW_SPHERE_LEVELS = (2, 3, 4, 5)
VIEW_TORUS_RES = ((24, 16), (48, 32), (96, 64))
GAUSS_RES = ((24, 16), (48, 32), (96, 64))


def case_id(c):
    return f"{c['family']}|{c['label']}|{c['variant']}|{c['seed']}"


def meta(c):
    return dict(case=case_id(c), family=c["family"], resolution=c["label"],
                variant=c["variant"], seed=c["seed"], required=c["required"])


def geodesic_sweep(fields):
    rows, qrows, srows = [], [], []
    for c in sphere_cases() + torus_cases():
        mesh = build(c)
        h = float(mesh.edge_lengths().mean())
        D = diameter(c)
        q = mesh_quality(mesh)
        qrows.append(dict(**meta(c), n_vertices=mesh.n_vertices, h_mean=h, **q))
        srows.append(dict(**meta(c), **structure_checks(c, mesh, eig=mesh.n_vertices <= 30000)))
        srcs = source_vertices(c, mesh)
        ang = interior_angles(mesh)
        valence = np.bincount(mesh.edges.ravel(), minlength=mesh.n_vertices)
        view = ((c["family"] == "sphere" and c["level"] in VIEW_SPHERE_LEVELS)
                or (c["family"] == "torus" and (c["nx"], c["ny"]) in VIEW_TORUS_RES))
        if view:
            fields[f"V|{case_id(c)}"] = mesh.vertices.astype(np.float32)
            fields[f"F|{c['family']}|{c['label']}"] = mesh.faces.astype(np.int32)
        for m in FIX["heat_method_m_sweep"]:
            solver = HeatMethodSolver(mesh, m=m)
            for name, src in srcs.items():
                phi = solver.distance(src)
                exact = exact_distance(c, mesh, src)
                sc = score_geodesic(phi, exact, h, D)
                local = float(np.degrees(ang[(mesh.faces == src).any(axis=1)].min()))
                rows.append(dict(**{**meta(c), "required": bool(
                                     c["required"] and m == FIX["heat_method_m"])},
                                 n_vertices=mesh.n_vertices, h_mean=h, D=D, m=m,
                                 source=name, source_vertex=src,
                                 source_valence=int(valence[src]),
                                 source_local_min_angle=local, **sc))
                if view:
                    fields[f"E|{case_id(c)}|{m:g}|{name}"] = (phi - exact).astype(np.float32)
        print(f"  {case_id(c):34s} V={mesh.n_vertices:6d} done", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(qrows), pd.DataFrame(srows)


def heat_sweep(fields):
    rows = []
    t = FIX["torus"]
    Lx, Ly, T = t["Lx"], t["Ly"], t["heat_modes"]["T"]
    for c in torus_cases(include_extra=False):
        mesh = build(c)
        L, _ = cotangent_laplacian(mesh); M = mass_matrix(mesh)
        h = float(mesh.edge_lengths().mean())
        for k, l in t["heat_modes"]["modes"]:
            lam = R.torus_mode_rate(k, l, Lx, Ly)
            u0 = R.torus_mode(mesh.vertices, k, l, Lx, Ly)
            a, res, drift, n = run_mode(mesh, u0, T, L, M)
            sc = score_mode(a, res, lam, T, h)
            if (k, l) == (0, 0):
                # constant mode: must stay constant, pointwise
                dev = abs(a - 1.0) + res
                sc["passed"] = bool(dev <= TOL["heat_mode"]["constant_mode_max_dev"])
                sc["residual"] = dev
            sc["passed"] = bool(sc["passed"] and drift <= TOL["heat_mode"]["total_heat_rel_drift_max"])
            rows.append(dict(**meta(c), h_mean=h, mode=f"({k},{l})", lam=lam, T=T, n_steps=n,
                             amplitude=a, amplitude_exact=float(np.exp(-lam * T)),
                             total_heat_drift=drift, **sc))
    s = FIX["sphere"]
    for c in sphere_cases(include_extra=False):
        mesh = build(c)
        L, _ = cotangent_laplacian(mesh); M = mass_matrix(mesh)
        h = float(mesh.edge_lengths().mean())
        for l in s["heat_modes"]["degrees"]:
            lam = R.sphere_mode_rate(l); T = s["heat_modes"]["T"]
            u0 = R.sphere_harmonic(mesh.vertices, l)
            a, res, drift, n = run_mode(mesh, u0, T, L, M)
            sc = score_mode(a, res, lam, T, h)
            rows.append(dict(**meta(c), h_mean=h, mode=f"l={l}", lam=lam, T=T, n_steps=n,
                             amplitude=a, amplitude_exact=float(np.exp(-lam * T)),
                             total_heat_drift=drift, **sc))
    df = pd.DataFrame(rows)

    # Gaussian bump next to the corner: the heat has to leave through both seams.
    g = t["gaussian"]
    center = np.array([g["center_frac"][0] * Lx, g["center_frac"][1] * Ly])
    grows = []
    for c in torus_cases(resolutions=[list(r) for r in GAUSS_RES]):
        if c["seed"] not in (0, t["seeds"][0]):
            continue
        mesh = build(c)
        L, _ = cotangent_laplacian(mesh); M = mass_matrix(mesh)
        w = M.diagonal()
        h = float(mesh.edge_lengths().mean())
        u = R.torus_gaussian(mesh.vertices, center, g["s0"], 0.0, Lx, Ly)
        t_now, total0 = 0.0, float((w * u).sum())
        for it, tt in enumerate(g["times"]):
            if tt > t_now:
                n = max(1, int(np.ceil((tt - t_now) / (h * h))))
                solver = HeatSolver(M, L, dt=(tt - t_now) / n)
                for _ in range(n):
                    u = solver.step(u)
                t_now = tt
            exact = R.torus_gaussian(mesh.vertices, center, g["s0"], tt, Lx, Ly)
            err = u - exact
            grows.append(dict(**meta(c), h_mean=h, t=tt, peak_exact=float(exact.max()),
                              peak_numeric=float(u.max()), linf=float(np.abs(err).max()),
                              linf_rel_peak=float(np.abs(err).max() / exact.max()),
                              l2_rel=float(np.sqrt((w * err ** 2).sum() / (w * exact ** 2).sum())),
                              total_heat=float((w * u).sum()),
                              total_heat_rel_drift=abs(float((w * u).sum()) - total0) / total0,
                              min_numeric=float(u.min())))
            fields[f"G|{case_id(c)}|{it}"] = err.astype(np.float32)
    return df, pd.DataFrame(grows)


def negative_controls():
    """Known-wrong answers scored with the real scoring functions."""
    rows = []
    t = FIX["torus"]
    Lx, Ly = t["Lx"], t["Ly"]
    D = R.torus_diameter(Lx, Ly)
    for c in torus_cases(include_extra=False):
        if c["variant"] != "regular":
            continue
        mesh = build(c)
        h = float(mesh.edge_lengths().mean())
        srcs = source_vertices(c, mesh)
        openm = open_rectangle(c["nx"], c["ny"], Lx, Ly)
        open_solver = HeatMethodSolver(openm, m=FIX["heat_method_m"])
        h_open = float(openm.edge_lengths().mean())
        for name, src in srcs.items():
            p = mesh.vertices[src]
            exact = R.torus_distance(mesh.vertices, p, Lx, Ly)
            wrong = R.unwrapped_distance(mesh.vertices, p)
            sc = score_geodesic(wrong, exact, h, D)
            rows.append(dict(control="unwrapped_reference", resolution=c["label"], source=name,
                             metric="linf_rel", value=sc["linf_rel"], tolerance=sc["tol_linf_rel"],
                             mae_rel=sc["mae_rel"], tol_mae_rel=sc["tol_mae_rel"],
                             rejected=not sc["passed"],
                             expected_rejection=name != "center"))
            osrc = int(np.argmin(np.hypot(openm.vertices[:, 0] - p[0], openm.vertices[:, 1] - p[1])))
            phi = open_solver.distance(osrc)
            ex_open = R.torus_distance(openm.vertices, openm.vertices[osrc], Lx, Ly)
            sc = score_geodesic(phi, ex_open, h_open, D)
            rows.append(dict(control="unstitched_mesh", resolution=c["label"], source=name,
                             metric="linf_rel", value=sc["linf_rel"], tolerance=sc["tol_linf_rel"],
                             mae_rel=sc["mae_rel"], tol_mae_rel=sc["tol_mae_rel"],
                             rejected=not sc["passed"],
                             expected_rejection=name != "center"))
        T = t["heat_modes"]["T"]
        for k, l in t["heat_modes"]["modes"]:
            if (k, l) == (0, 0):
                continue
            lam = R.torus_mode_rate(k, l, Lx, Ly)
            sc = score_mode(float(np.exp(-2 * lam * T)), 0.0, lam, T, h)
            rows.append(dict(control="wrong_decay", resolution=c["label"], source=f"({k},{l})",
                             metric="rate_rel_err", value=sc["rate_rel_err"], tolerance=sc["tol"],
                             rejected=not sc["passed"],
                             expected_rejection=sc["tol"] < 1.0))
    return pd.DataFrame(rows)


def environment():
    def git(*a):
        try:
            return subprocess.check_output(["git", *a], cwd=os.path.dirname(RESULTS)).decode().strip()
        except Exception:
            return "unknown"
    return dict(commit=git("rev-parse", "HEAD"), dirty=bool(git("status", "--porcelain", "--", "ddg", "benchmarks")),
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                platform=platform.platform(), machine=platform.machine(),
                generated=time.strftime("%Y-%m-%d %H:%M:%S %Z"))


def main():
    os.makedirs(RESULTS, exist_ok=True)
    fields = {}
    t0 = time.perf_counter()
    print("[1/3] geodesic sweep + structure + quality")
    geo, qual, struct = geodesic_sweep(fields)
    print("[2/3] heat modes + Gaussian")
    heat, gauss = heat_sweep(fields)
    print("[3/3] negative controls")
    neg = negative_controls()
    geo.to_csv(os.path.join(RESULTS, "geodesic.csv"), index=False)
    qual.to_csv(os.path.join(RESULTS, "quality.csv"), index=False)
    struct.to_csv(os.path.join(RESULTS, "structure.csv"), index=False)
    heat.to_csv(os.path.join(RESULTS, "heat_modes.csv"), index=False)
    gauss.to_csv(os.path.join(RESULTS, "heat_gaussian.csv"), index=False)
    neg.to_csv(os.path.join(RESULTS, "negative_controls.csv"), index=False)
    np.savez_compressed(os.path.join(RESULTS, "fields.npz"), **fields)
    req = geo[geo.required]
    summary = dict(
        environment=environment(), runtime_s=time.perf_counter() - t0,
        geodesic=dict(n=len(geo), required=len(req), required_failed=int((~req.passed).sum()),
                      stress=int((~geo.required).sum()),
                      stress_failed=int((~geo[~geo.required].passed).sum())),
        heat_modes=dict(n=len(heat), failed=int((~heat.passed).sum()),
                        vacuous=int(heat.vacuous.sum())),
        structure=dict(n=len(struct), failed=int((~struct.passed).sum())),
        negative_controls=dict(n=len(neg), expected=int(neg.expected_rejection.sum()),
                               rejected_of_expected=int((neg.rejected & neg.expected_rejection).sum())),
    )
    with open(os.path.join(RESULTS, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
