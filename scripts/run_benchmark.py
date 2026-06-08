#!/usr/bin/env python3
"""End-to-end DDG benchmark driver.

Runs the full pipeline described in the project brief and writes results to
``results/`` (CSV/JSON) and ``figures/`` (PNG).  Phases can be selected
individually; by default everything runs.

Examples
--------
    python scripts/run_benchmark.py --all
    python scripts/run_benchmark.py --phase stats heat geodesic
    python scripts/run_benchmark.py --phase convergence --max-sub 7
    python scripts/run_benchmark.py --phase performance --max-sub 7

Phases
------
    stats         mesh family statistics + d^2 = 0 check
    heat          heat diffusion from a point source (snapshots)
    geodesic      heat-method distance field + error vs arccos
    convergence   error vs mesh resolution, empirical rate
    spectral      first eigenfunctions + spectrum vs l(l+1)
    performance   assembly / solve / eig timing + memory scaling
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import tracemalloc

import numpy as np
import pandas as pd
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
    HeatSolver,
    point_source,
    HeatMethodSolver,
    eigenpairs,
    sphere_reference_spectrum,
    degree_for_index,
    exact_sphere_distance,
    error_metrics,
    convergence_rate,
)
from ddg import visualize as viz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
RES = os.path.join(ROOT, "results")
os.makedirs(FIG, exist_ok=True)
os.makedirs(RES, exist_ok=True)


def _farthest_vertex(mesh, source):
    """Return the index of the vertex antipodal-most to ``source``."""
    d = mesh.vertices @ mesh.vertices[source]
    return int(np.argmin(d))


# ---------------------------------------------------------------------------
# Phase 1/2: mesh statistics and topology check
# ---------------------------------------------------------------------------
def phase_stats(subs):
    rows = []
    for n in subs:
        mesh = icosphere(n)
        s = mesh.statistics()
        s["subdivisions"] = n
        # d^2 = 0 correctness check
        d0 = build_d0(mesh)
        d1 = build_d1(mesh)
        s["d2_residual"] = verify_d_squared(d0, d1)
        s["total_area"] = float(face_areas(mesh).sum())
        s["area_error_vs_4pi"] = abs(s["total_area"] - 4 * np.pi)
        rows.append(s)
        print(f"  sub={n}: V={s['n_vertices']} E={s['n_edges']} "
              f"F={s['n_faces']} chi={s['euler_characteristic']} "
              f"h_mean={s['mean_edge_length']:.4f} d2={s['d2_residual']:.1e}")
    df = pd.DataFrame(rows)
    out = os.path.join(RES, "mesh_statistics.csv")
    df.to_csv(out, index=False)
    print(f"  -> {out}")
    return df


# ---------------------------------------------------------------------------
# Phase 4: heat diffusion
# ---------------------------------------------------------------------------
def phase_heat(sub=5):
    mesh = icosphere(sub)
    L, _ = cotangent_laplacian(mesh)
    M = mass_matrix(mesh)
    h = mesh.edge_lengths().mean()
    dt = h * h
    solver = HeatSolver(M, L, dt=dt)
    u = point_source(mesh.n_vertices, 0)

    record = [0, 1, 3, 10, 30, 100]
    snaps, times = [], []
    cur = u.copy()
    step = 0
    for target in record:
        while step < target:
            cur = solver.step(cur)
            step += 1
        snaps.append(cur.copy())
        times.append(step * dt)

    viz.plot_heat_series(
        mesh, snaps, times,
        out_path=os.path.join(FIG, "heat_diffusion_series.png"),
    )
    # final field with isotherms
    viz.plot_scalar_on_sphere(
        mesh, snaps[-1], title=f"heat field at t={times[-1]:.3g}",
        out_path=os.path.join(FIG, "heat_field_contours.png"),
        cmap="inferno", contours=12,
    )
    print(f"  V={mesh.n_vertices}, dt={dt:.2e}, "
          f"figures -> heat_diffusion_series.png, heat_field_contours.png")


# ---------------------------------------------------------------------------
# Phase 5/6: geodesic distance + error
# ---------------------------------------------------------------------------
def phase_geodesic(sub=6):
    mesh = icosphere(sub)
    source = 0
    solver = HeatMethodSolver(mesh, m=1.0)
    phi = solver.distance(source)
    exact = exact_sphere_distance(mesh.vertices, mesh.vertices[source])
    metrics = error_metrics(phi, exact)
    print(f"  V={mesh.n_vertices}  MAE={metrics['MAE']:.4e}  "
          f"Linf={metrics['Linf']:.4e}  RMSE={metrics['RMSE']:.4e}")

    viz.plot_scalar_on_sphere(
        mesh, phi, title="heat-method geodesic distance",
        out_path=os.path.join(FIG, "geodesic_distance.png"),
        cmap="viridis", contours=15,
    )
    viz.plot_scalar_on_sphere(
        mesh, np.abs(phi - exact), title="|numeric - exact| distance error",
        out_path=os.path.join(FIG, "geodesic_error.png"),
        cmap="magma", contours=12,
    )
    with open(os.path.join(RES, "geodesic_metrics.json"), "w") as fh:
        json.dump({"n_vertices": mesh.n_vertices, **metrics}, fh, indent=2)
    print("  figures -> geodesic_distance.png, geodesic_error.png")


# ---------------------------------------------------------------------------
# Phase 6: convergence study
# ---------------------------------------------------------------------------
def phase_convergence(subs):
    rows = []
    for n in subs:
        mesh = icosphere(n)
        solver = HeatMethodSolver(mesh, m=1.0)
        phi = solver.distance(0)
        exact = exact_sphere_distance(mesh.vertices, mesh.vertices[0])
        m = error_metrics(phi, exact)
        h = mesh.edge_lengths().mean()
        rows.append({"subdivisions": n, "n_vertices": mesh.n_vertices,
                     "h_mean": h, **m})
        print(f"  sub={n} V={mesh.n_vertices} h={h:.4f} "
              f"MAE={m['MAE']:.3e} Linf={m['Linf']:.3e} RMSE={m['RMSE']:.3e}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "convergence.csv"), index=False)

    h = df["h_mean"].values
    rates = {k: convergence_rate(h, df[k].values) for k in ("MAE", "Linf", "RMSE")}
    print(f"  empirical rates (error ~ h^p):  "
          + "  ".join(f"{k}: p={v:.2f}" for k, v in rates.items()))

    viz.plot_convergence(
        h, {k: df[k].values for k in ("MAE", "Linf", "RMSE")},
        title="Heat-method geodesic error vs mesh size",
        out_path=os.path.join(FIG, "convergence_vs_h.png"),
    )
    viz.plot_convergence(
        df["n_vertices"].values,
        {k: df[k].values for k in ("MAE", "Linf", "RMSE")},
        title="Heat-method geodesic error vs vertex count",
        xlabel="number of vertices",
        out_path=os.path.join(FIG, "convergence_vs_n.png"),
        reference_slopes=(),
    )
    with open(os.path.join(RES, "convergence_rates.json"), "w") as fh:
        json.dump(rates, fh, indent=2)
    print("  figures -> convergence_vs_h.png, convergence_vs_n.png")
    return df, rates


# ---------------------------------------------------------------------------
# Phase 7: spectral analysis
# ---------------------------------------------------------------------------
def phase_spectral(sub=5, k=20):
    mesh = icosphere(sub)
    L, _ = cotangent_laplacian(mesh)
    M = mass_matrix(mesh)
    vals, vecs = eigenpairs(L, M, k=k)

    # compare against analytical l(l+1) spectrum
    max_l = degree_for_index(k - 1)
    ref = sphere_reference_spectrum(max_l)[:k]
    rows = [{"index": i, "lambda_num": float(vals[i]),
             "degree_l": degree_for_index(i),
             "lambda_exact": float(ref[i]),
             "abs_error": float(abs(vals[i] - ref[i]))}
            for i in range(k)]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "spectrum.csv"), index=False)
    print(df.to_string(index=False))

    viz.plot_eigenfunctions(
        mesh, vecs, n=k,
        out_path=os.path.join(FIG, "eigenfunctions.png"),
    )

    # spectrum scatter
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(range(k), ref, "s--", label="exact  l(l+1)", color="gray")
    ax.plot(range(k), vals, "o-", label="numerical")
    ax.set_xlabel("eigenvalue index"); ax.set_ylabel("eigenvalue")
    ax.set_title("Laplace-Beltrami spectrum vs sphere reference")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "spectrum.png"), dpi=130)
    print("  figures -> eigenfunctions.png, spectrum.png")


# ---------------------------------------------------------------------------
# Phase 8: performance benchmarking
# ---------------------------------------------------------------------------
def phase_performance(subs, k_eig=10):
    rows = []
    for n in subs:
        tracemalloc.start()
        t0 = time.perf_counter()
        mesh = icosphere(n)
        t_mesh = time.perf_counter() - t0

        t0 = time.perf_counter()
        L, _ = cotangent_laplacian(mesh)
        M = mass_matrix(mesh)
        t_assembly = time.perf_counter() - t0

        h = mesh.edge_lengths().mean()
        t0 = time.perf_counter()
        solver = HeatSolver(M, L, dt=h * h)
        u = solver.step(point_source(mesh.n_vertices, 0))
        t_solve = time.perf_counter() - t0

        t0 = time.perf_counter()
        try:
            eigenpairs(L, M, k=min(k_eig, mesh.n_vertices - 2))
            t_eig = time.perf_counter() - t0
        except Exception:
            t_eig = float("nan")

        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        rows.append({
            "subdivisions": n, "n_vertices": mesh.n_vertices,
            "n_faces": mesh.n_faces, "nnz_L": int(L.nnz),
            "t_mesh_s": t_mesh, "t_assembly_s": t_assembly,
            "t_solve_s": t_solve, "t_eig_s": t_eig,
            "peak_mem_mb": peak / 1e6,
        })
        print(f"  sub={n} V={mesh.n_vertices}: assembly={t_assembly:.3f}s "
              f"solve={t_solve:.3f}s eig={t_eig:.3f}s mem={peak/1e6:.1f}MB")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(RES, "performance.csv"), index=False)

    import matplotlib.pyplot as plt
    nv = df["n_vertices"].values
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for col, lab in [("t_assembly_s", "assembly"), ("t_solve_s", "linear solve"),
                     ("t_eig_s", "eig (k modes)")]:
        axes[0].loglog(nv, df[col].values, "o-", label=lab)
    # O(n) reference
    axes[0].loglog(nv, df["t_assembly_s"].values[0] * nv / nv[0], "--",
                   color="gray", alpha=0.6, label="O(n)")
    axes[0].set_xlabel("vertices"); axes[0].set_ylabel("seconds")
    axes[0].set_title("Runtime scaling"); axes[0].legend(); axes[0].grid(alpha=0.3, which="both")

    axes[1].loglog(nv, df["peak_mem_mb"].values, "o-", color="C3", label="peak memory")
    axes[1].loglog(nv, df["peak_mem_mb"].values[0] * nv / nv[0], "--",
                   color="gray", alpha=0.6, label="O(n)")
    axes[1].set_xlabel("vertices"); axes[1].set_ylabel("MB")
    axes[1].set_title("Memory scaling"); axes[1].legend(); axes[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "performance_scaling.png"), dpi=130)
    print("  figure -> performance_scaling.png")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
PHASES = ["stats", "heat", "geodesic", "convergence", "spectral", "performance"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", nargs="+", choices=PHASES, default=None,
                    help="phases to run (default: all)")
    ap.add_argument("--all", action="store_true", help="run all phases")
    ap.add_argument("--max-sub", type=int, default=6,
                    help="finest icosphere subdivision level for sweeps")
    ap.add_argument("--min-sub", type=int, default=2,
                    help="coarsest icosphere subdivision level for sweeps")
    args = ap.parse_args()

    phases = PHASES if (args.all or args.phase is None) else args.phase
    subs = list(range(args.min_sub, args.max_sub + 1))

    print(f"Running phases: {phases}")
    print(f"Subdivision sweep: {subs} "
          f"(V = {[10*4**n+2 for n in subs]})\n")

    if "stats" in phases:
        print("[Phase 1/2] Mesh statistics & topology (d^2=0):")
        phase_stats(subs)
    if "heat" in phases:
        print("\n[Phase 4] Heat diffusion:")
        phase_heat(sub=min(5, args.max_sub))
    if "geodesic" in phases:
        print("\n[Phase 5/6] Geodesic distance (heat method):")
        phase_geodesic(sub=min(6, args.max_sub))
    if "convergence" in phases:
        print("\n[Phase 6] Convergence study:")
        phase_convergence(subs)
    if "spectral" in phases:
        print("\n[Phase 7] Spectral analysis:")
        phase_spectral(sub=min(5, args.max_sub))
    if "performance" in phases:
        print("\n[Phase 8] Performance benchmarking:")
        phase_performance(subs)

    print("\nDone.")


if __name__ == "__main__":
    main()
