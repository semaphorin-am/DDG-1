#!/usr/bin/env python3
"""Baseline (SuperLU COLAMD) versus optimised (nested dissection) cost.

Every measurement runs in a fresh single-threaded subprocess so that peak
memory (``ru_maxrss``) belongs to one configuration only.  Orderings are
interleaved across repeats.  Setup (two sparse factorizations) and per-query
(two back-substitutions + gradient/divergence) are timed separately.
Writes ``perf_runs.csv`` (every run) and ``perf_summary.csv`` (median/min/max
plus the verdict defined in ``tolerances.json``).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

ORDERINGS = ("colamd", "nested")
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")


def worker(spec: dict) -> None:
    import resource
    import time
    import numpy as np
    from ddg import HeatMethodSolver
    from benchmarks.common import FIX, build, diameter, exact_distance, source_vertices, score_geodesic

    def rss_mb():
        r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return r / (1024.0 * 1024.0) if sys.platform == "darwin" else r / 1024.0

    case = spec["case"]
    mesh = build(case)
    srcs = list(source_vertices(case, mesh).values())[: spec["queries"]]
    h = float(mesh.edge_lengths().mean())
    rss0 = rss_mb()
    t0 = time.perf_counter()
    solver = HeatMethodSolver(mesh, m=FIX["heat_method_m"], ordering=spec["ordering"])
    setup = time.perf_counter() - t0
    rss1 = rss_mb()
    q = []
    phis = []
    for s in srcs:
        t0 = time.perf_counter()
        phis.append(solver.distance(s))
        q.append(time.perf_counter() - t0)
    rss2 = rss_mb()
    sc = [score_geodesic(p, exact_distance(case, mesh, s), h, diameter(case)) for p, s in zip(phis, srcs)]
    if spec.get("save"):
        np.save(spec["save"], np.array(phis))
    print(json.dumps(dict(n_vertices=mesh.n_vertices, setup_s=setup, query_s=float(np.median(q)),
                          query_total_s=float(np.sum(q)), rss_before_mb=rss0, peak_rss_mb=rss2,
                          solver_peak_mb=rss2 - rss0, setup_peak_mb=rss1 - rss0,
                          factor_nnz=int(solver.factor_nnz),
                          mae_rel_mean=float(np.mean([s["mae_rel"] for s in sc])),
                          linf_rel_max=float(np.max([s["linf_rel"] for s in sc])),
                          n_pass=int(sum(s["passed"] for s in sc)))))


def main() -> None:
    import numpy as np
    import pandas as pd
    from benchmarks.common import FIX, TOL, RESULTS

    P = FIX["performance"]
    env = dict(os.environ, **{v: str(P["threads"]) for v in THREAD_VARS})
    rows = []
    tmp = tempfile.mkdtemp()
    for mi, mspec in enumerate(P["meshes"]):
        case = dict(mspec, required=True)
        if case["family"] == "sphere":
            case["label"] = f"L{case['level']}"
        else:
            case["nx"], case["ny"] = case["resolution"]
            case["label"] = f"{case['nx']}x{case['ny']}"
        name = f"{case['family']} {case['label']} {case['variant']}"
        saved = {}
        for rep in range(P["repeats"]):
            for o in (ORDERINGS if rep % 2 == 0 else ORDERINGS[::-1]):
                spec = dict(case=case, ordering=o, queries=P["queries_per_repeat"])
                if rep == 0:
                    spec["save"] = saved[o] = os.path.join(tmp, f"{mi}_{o}.npy")
                out = subprocess.check_output([sys.executable, os.path.abspath(__file__), "--worker",
                                               json.dumps(spec)], env=env, cwd=ROOT)
                rows.append(dict(mesh=name, family=case["family"], resolution=case["label"],
                                 variant=case["variant"], seed=case["seed"], ordering=o, repeat=rep,
                                 threads=P["threads"], **json.loads(out.decode().strip().splitlines()[-1])))
        a, b = np.load(saved["colamd"]), np.load(saved["nested"])
        diff = float(np.abs(a - b).max())
        for r in rows:
            if r["mesh"] == name:
                r["max_field_diff"] = diff
        print(f"  {name:28s} max|phi_nested - phi_colamd| = {diff:.2e}", flush=True)
    runs = pd.DataFrame(rows)
    runs.to_csv(os.path.join(RESULTS, "perf_runs.csv"), index=False)

    out = []
    for (mesh, metric), _ in [((m, k), 0) for m in runs.mesh.unique()
                              for k in ("setup_s", "query_s", "solver_peak_mb")]:
        g = runs[runs.mesh == mesh]
        base, opt = g[g.ordering == "colamd"][metric], g[g.ordering == "nested"][metric]
        ratio = float(opt.median() / base.median())
        separated = bool(opt.max() < base.min() or opt.min() > base.max())
        verdict = "no reliable difference"
        if separated and ratio <= 0.8:
            verdict = "win"
        elif separated and ratio > 1.05:
            verdict = "regression"
        elif separated and ratio < 1.0:
            verdict = "small improvement (<20%)"
        f = g.iloc[0]
        gb, go = g[g.ordering == "colamd"], g[g.ordering == "nested"]
        out.append(dict(mesh=mesh, n_vertices=int(f.n_vertices), metric=metric,
                        baseline_median=float(base.median()), baseline_min=float(base.min()),
                        baseline_max=float(base.max()), optimized_median=float(opt.median()),
                        optimized_min=float(opt.min()), optimized_max=float(opt.max()),
                        ratio=ratio, ranges_separated=separated, verdict=verdict,
                        max_field_diff=float(f.max_field_diff),
                        equal_accuracy=bool(f.max_field_diff <= 1e-8 and gb.n_pass.iloc[0] == go.n_pass.iloc[0]),
                        baseline_mae_rel=float(gb.mae_rel_mean.iloc[0]),
                        optimized_mae_rel=float(go.mae_rel_mean.iloc[0]),
                        baseline_factor_nnz=int(gb.factor_nnz.iloc[0]),
                        optimized_factor_nnz=int(go.factor_nnz.iloc[0]),
                        repeats=P["repeats"], threads=P["threads"]))
    summ = pd.DataFrame(out)
    summ.to_csv(os.path.join(RESULTS, "perf_summary.csv"), index=False)
    print(summ[["mesh", "metric", "baseline_median", "optimized_median", "ratio", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--worker":
        worker(json.loads(sys.argv[2]))
    else:
        main()
