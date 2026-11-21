#!/usr/bin/env python3
"""Static figures for readers who do not open the HTML report.

Everything is drawn from ``results/beyond_sphere/*.csv`` and ``fields.npz``.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from benchmarks import references as R
from benchmarks.common import FIX, RESULTS, FIGURES

COL = dict(regular="#2456c8", mild="#19806f", strong="#c2571a", jittered="#c2571a")
rd = lambda n: pd.read_csv(os.path.join(RESULTS, n))


def save(fig, name):
    fig.savefig(os.path.join(FIGURES, name), dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("  figures/beyond_sphere/" + name)


def convergence(geo):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex="col")
    g = geo[geo.m == 1.0]
    for j, fam in enumerate(("sphere", "torus")):
        for i, (met, tol, lab) in enumerate((("mae_rel", "tol_mae_rel", "mean error / diameter"),
                                             ("linf_rel", "tol_linf_rel", "worst-point error / diameter"))):
            ax = axes[i, j]
            f = g[g.family == fam]
            for v in f.variant.unique():
                q = f[f.variant == v]
                ax.scatter(q.n_vertices, q[met], s=9, alpha=0.35, color=COL[v])
                mu = q.groupby("n_vertices")[met].mean()
                ax.plot(mu.index, mu.values, "-o", ms=3, color=COL[v], label=v)
            t = f[f.variant == "regular"].groupby("n_vertices")[tol].first()
            ax.plot(t.index, t.values, "--", color="#b3261e", label="tolerance")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(alpha=0.3, which="both")
            ax.set_title(f"{'flat torus' if fam == 'torus' else 'sphere'}: {lab}", fontsize=10)
            if i == 1:
                ax.set_xlabel("vertices")
            ax.legend(fontsize=8)
    fig.suptitle("Heat-method geodesic error vs resolution (m = 1; dots = seeds x sources)")
    save(fig, "convergence.png")


def quality(geo, qual):
    g = geo[(geo.family == "sphere") & (geo.m == 1.0)]
    case = g.groupby("case").mae_rel.mean()
    q = qual[qual.family == "sphere"].set_index("case").join(case)
    base = q[q.variant == "regular"].set_index("resolution").mae_rel
    q["ratio"] = q.mae_rel / q.resolution.map(base)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for v in ("regular", "mild", "strong"):
        s = q[q.variant == v]
        axes[0].scatter(100 * s.frac_obtuse, s.ratio, color=COL[v], label=v)
        axes[1].scatter(s.min_angle_min, s.ratio, color=COL[v], label=v)
    axes[0].set_xlabel("obtuse triangles (%)"); axes[1].set_xlabel("smallest angle in mesh (degrees)")
    for ax in axes:
        ax.set_ylabel("mean error / regular mesh at same level"); ax.grid(alpha=0.3); ax.legend()
    fig.suptitle("Sphere: triangle quality vs geodesic error (one dot per mesh)")
    save(fig, "quality_vs_error.png")


def seam(geo, z):
    t = FIX["torus"]; Lx, Ly = t["Lx"], t["Ly"]
    key = "torus|48x32|regular|0"
    V, F = z["V|" + key].astype(float), z["F|torus|48x32"]
    row = geo[(geo.case == key) & (geo.m == 1.0) & (geo.source == "corner")].iloc[0]
    src = V[int(row.source_vertex)]
    err = z[f"E|{key}|1|corner"].astype(float)
    exact = R.torus_distance(V, src, Lx, Ly)
    wrong = R.unwrapped_distance(V, src)
    nx, ny = 48, 32
    ext = [0, Lx, 0, Ly]
    grid = lambda a: a.reshape(ny, nx) if abs(V[1, 0] - V[0, 0]) > 0 else a.reshape(nx, ny).T
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.4))
    for ax, a, title, cm in ((axes[0], exact, "exact wrapped distance", "viridis"),
                             (axes[1], exact + err, "computed (heat method)", "viridis"),
                             (axes[2], wrong, "unwrapped distance (wrong; rejected)", "viridis")):
        im = ax.imshow(grid(a), origin="lower", extent=ext, cmap=cm, vmin=0, vmax=wrong.max(), aspect="equal")
        ax.set_title(title, fontsize=10); ax.plot(*src, "o", color="orange", mec="w")
        fig.colorbar(im, ax=ax, shrink=0.8)
    lim = np.abs(err).max()
    im = axes[3].imshow(grid(err), origin="lower", extent=ext, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="equal")
    axes[3].set_title("signed error (computed - exact)", fontsize=10); axes[3].plot(*src, "o", color="orange", mec="w")
    fig.colorbar(im, ax=axes[3], shrink=0.8)
    fig.suptitle("Flat torus 48x32, source next to the corner: distance wraps through both seams")
    save(fig, "torus_seam.png")


def heat(hm, gauss):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    f = hm[(hm.family == "torus") & (hm.lam > 0)]
    for (mode, v), q in f.groupby(["mode", "variant"]):
        mu = q.groupby("h_mean").rate_rel_err.mean()
        axes[0].loglog(mu.index, mu.values, "-o" if v == "regular" else "--s", ms=3, label=f"{mode} {v}")
    h = np.array(sorted(f.h_mean.unique()))
    axes[0].loglog(h, 0.5 * f.rate_rel_err.min() * (h / h[0]) ** 2, "k:", label="slope 2")
    axes[0].set_xlabel("mean edge length h"); axes[0].set_ylabel("relative decay-rate error")
    axes[0].set_title("Torus Fourier modes: decay rate vs exact", fontsize=10); axes[0].legend(fontsize=6, ncol=2); axes[0].grid(alpha=0.3, which="both")
    for case, q in gauss.groupby("case"):
        axes[1].semilogy(q.t, q.linf_rel_peak.clip(lower=1e-6), "-o", ms=3, label=case.replace("torus|", "").replace("|", " "))
    axes[1].set_xlabel("diffusion time t"); axes[1].set_ylabel("worst error / exact peak")
    axes[1].set_title("Gaussian next to the corner vs periodised exact solution", fontsize=10); axes[1].legend(fontsize=7); axes[1].grid(alpha=0.3, which="both")
    save(fig, "heat_accuracy.png")


def perf(ps, runs):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    meshes = list(ps.mesh.unique())
    for ax, (met, lab) in zip(axes, (("setup_s", "setup (s)"), ("query_s", "per query (s)"), ("solver_peak_mb", "peak solver memory (MB)"))):
        for k, (o, c, name) in enumerate((("colamd", "#2456c8", "baseline"), ("nested", "#c2571a", "nested dissection"))):
            med = [runs[(runs.mesh == m) & (runs.ordering == o)][met].median() for m in meshes]
            ax.barh(np.arange(len(meshes)) + (0.2 - 0.4 * k), med, 0.38, color=c, label=name, alpha=0.85)
            for i, m in enumerate(meshes):
                r = runs[(runs.mesh == m) & (runs.ordering == o)][met]
                ax.plot(r, np.full(len(r), i + (0.2 - 0.4 * k)), "|", color="k", ms=8)
        ax.set_yticks(range(len(meshes))); ax.set_yticklabels(meshes if ax is axes[0] else []); ax.invert_yaxis()
        ax.set_xscale("log"); ax.set_xlabel(lab); ax.grid(alpha=0.3, axis="x", which="both")
        for i, m in enumerate(meshes):
            r = ps[(ps.mesh == m) & (ps.metric == met)].iloc[0]
            ax.text(ax.get_xlim()[1], i, f" {r.ratio:.2f}x {r.verdict}", va="center", fontsize=7)
    axes[0].legend(fontsize=8)
    fig.suptitle("Baseline vs optimised at identical accuracy (median of repeats, ticks = each run, 1 thread)")
    save(fig, "cost.png")


def main():
    os.makedirs(FIGURES, exist_ok=True)
    geo = rd("geodesic.csv")
    convergence(geo); quality(geo, rd("quality.csv")); seam(geo, np.load(os.path.join(RESULTS, "fields.npz")))
    heat(rd("heat_modes.csv"), rd("heat_gaussian.csv")); perf(rd("perf_summary.csv"), rd("perf_runs.csv"))


if __name__ == "__main__":
    main()
