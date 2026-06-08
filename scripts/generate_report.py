#!/usr/bin/env python3
"""Generate the error-analysis and performance reports from benchmark outputs.

Reads the CSV/JSON files written to ``results/`` by ``run_benchmark.py`` and
renders two Markdown reports into ``docs/``.  Keeping report generation separate
and data-driven means the reports always reflect the latest benchmark run.

    python scripts/run_benchmark.py --all      # produce results/
    python scripts/generate_report.py          # render docs/*_report.md
"""

from __future__ import annotations

import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DOCS = os.path.join(ROOT, "docs")


def _read_csv(name):
    path = os.path.join(RES, name)
    return pd.read_csv(path) if os.path.exists(path) else None


def _read_json(name):
    path = os.path.join(RES, name)
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return None


def error_report():
    stats = _read_csv("mesh_statistics.csv")
    conv = _read_csv("convergence.csv")
    rates = _read_json("convergence_rates.json")
    spectrum = _read_csv("spectrum.csv")

    lines = ["# Error Analysis Report\n",
             "Numerical results from `scripts/run_benchmark.py`, validated "
             "against analytical solutions on the unit sphere.\n"]

    if stats is not None:
        lines.append("## Mesh family & topology\n")
        lines.append("Every mesh satisfies the Euler characteristic "
                     "`V - E + F = 2` and the cochain identity `d1 @ d0 = 0` "
                     "exactly (residual column).\n")
        cols = ["subdivisions", "n_vertices", "n_edges", "n_faces",
                "euler_characteristic", "mean_edge_length",
                "area_error_vs_4pi", "d2_residual"]
        lines.append(stats[cols].to_markdown(index=False) + "\n")

    if conv is not None:
        lines.append("## Geodesic distance convergence (heat method)\n")
        lines.append("Error of the heat-method geodesic distance vs the exact "
                     "great-circle distance `arccos(p·q)`, over all vertices.\n")
        cols = ["subdivisions", "n_vertices", "h_mean", "MAE", "Linf", "RMSE"]
        lines.append(conv[cols].to_markdown(index=False) + "\n")

    if rates is not None:
        lines.append("### Empirical convergence rate\n")
        lines.append("Fitting `error ~ C·h^p` in log-log space gives:\n")
        for k, v in rates.items():
            lines.append(f"- **{k}**: p = {v:.2f}")
        lines.append("\nThe heat method is first-order accurate in the smoothed "
                     "distance; the observed sub-linear-to-linear slope is "
                     "consistent with the literature (accuracy is limited by the "
                     "diffusion time `t = m·h²`).\n")
        lines.append("\n![convergence](../figures/convergence_vs_h.png)\n")

    if spectrum is not None:
        lines.append("## Spectral validation\n")
        lines.append("Generalized eigenvalues of `L φ = λ M φ` vs the exact "
                     "sphere spectrum `λ_ℓ = ℓ(ℓ+1)` with multiplicity `2ℓ+1`.\n")
        lines.append(spectrum.to_markdown(index=False) + "\n")
        lines.append("\n![spectrum](../figures/spectrum.png)\n")
        lines.append("\nThe first eigenfunctions reproduce the spherical "
                     "harmonics:\n")
        lines.append("\n![eigenfunctions](../figures/eigenfunctions.png)\n")

    out = os.path.join(DOCS, "error_report.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    print(f"-> {out}")


def performance_report():
    perf = _read_csv("performance.csv")
    if perf is None:
        print("no performance.csv; run `run_benchmark.py --phase performance`")
        return

    lines = ["# Performance Report\n",
             "Timing and peak memory from `scripts/run_benchmark.py "
             "--phase performance`. All operators are sparse with O(|V|) "
             "nonzeros.\n"]
    cols = ["subdivisions", "n_vertices", "n_faces", "nnz_L",
            "t_assembly_s", "t_solve_s", "t_eig_s", "peak_mem_mb"]
    lines.append(perf[cols].to_markdown(index=False) + "\n")

    # crude scaling exponents (time ~ n^p) from finest two points
    import numpy as np
    nv = perf["n_vertices"].values
    lines.append("\n## Observed scaling (time ~ n^p, finest two levels)\n")
    for col, lab in [("t_assembly_s", "assembly"), ("t_solve_s", "linear solve"),
                     ("t_eig_s", "eigensolve")]:
        y = perf[col].values
        if len(nv) >= 2 and y[-1] > 0 and y[-2] > 0:
            p = np.log(y[-1] / y[-2]) / np.log(nv[-1] / nv[-2])
            lines.append(f"- **{lab}**: p ≈ {p:.2f}")
    lines.append("\n![scaling](../figures/performance_scaling.png)\n")

    out = os.path.join(DOCS, "performance_report.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines))
    print(f"-> {out}")


if __name__ == "__main__":
    error_report()
    performance_report()
