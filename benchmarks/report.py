#!/usr/bin/env python3
"""Build the self-contained HTML report from the raw benchmark outputs.

Reads only ``results/beyond_sphere/*`` plus fixtures/tolerances; every number
and picture in the report comes from those files.  The output is one HTML file
with all data, script and style inline: no network, no server.
"""

from __future__ import annotations

import base64
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from benchmarks.common import FIX, TOL, RESULTS, REPORT, HERE

FIELD_M = {"sphere": (1.0, 4.0, 16.0), "torus": (1.0, 4.0)}
FIELD_RES = {"sphere": ("L2", "L3", "L4"), "torus": ("24x16", "48x32", "96x64")}


def b64(a: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(a).tobytes()).decode()


def quant(a: np.ndarray) -> dict:
    """16-bit quantisation of an error field (relative precision 3e-5)."""
    s = float(np.abs(a).max()) or 1.0
    return dict(s=s, d=b64(np.round(a / s * 32767).astype("<i2")))


def records(df: pd.DataFrame, cols) -> dict:
    df = df[list(cols)]
    return dict(cols=list(cols), rows=json.loads(df.to_json(orient="values", double_precision=15)))


def main() -> None:
    rd = lambda n: pd.read_csv(os.path.join(RESULTS, n))
    geo, qual, struct = rd("geodesic.csv"), rd("quality.csv"), rd("structure.csv")
    heat, gauss, neg = rd("heat_modes.csv"), rd("heat_gaussian.csv"), rd("negative_controls.csv")
    perf, runs = rd("perf_summary.csv"), rd("perf_runs.csv")
    with open(os.path.join(RESULTS, "summary.json")) as fh:
        summary = json.load(fh)
    with open(os.path.join(RESULTS, "baseline.json")) as fh:
        baseline = json.load(fh)
    z = np.load(os.path.join(RESULTS, "fields.npz"))

    meshes, faces, fields, gfields = {}, {}, {}, {}
    for k in z.files:
        p = k.split("|")
        if p[0] == "F" and p[2] in FIELD_RES[p[1]]:
            faces[f"{p[1]}|{p[2]}"] = b64(z[k].astype("<u2"))
        elif p[0] == "V" and p[2] in FIELD_RES[p[1]]:
            meshes["|".join(p[1:])] = b64(z[k].astype("<f4"))
        elif p[0] == "E" and p[2] in FIELD_RES[p[1]] and float(p[5]) in FIELD_M[p[1]]:
            fields["|".join(p[1:])] = quant(z[k])
        elif p[0] == "G":
            gfields["|".join(p[1:])] = quant(z[k])

    # Worst viewable case per family, relative to its own tolerance (annotation targets).
    view = geo[(geo.m == 1.0) & geo.apply(lambda r: r.resolution in FIELD_RES[r.family], axis=1)].copy()
    view["ratio"] = np.maximum(view.mae_rel / view.tol_mae_rel, view.linf_rel / view.tol_linf_rel)
    ann = {}
    for fam in ("sphere", "torus"):
        v = view[view.family == fam]
        if fam == "sphere":
            v = v[v.resolution == "L4"]
        else:
            v = v[v.resolution == "96x64"]
        r = v.sort_values("ratio").iloc[-1]
        ann[fam] = dict(case=r.case, source=r.source, m=float(r.m), ratio=float(r.ratio))

    gcols = ["case", "family", "resolution", "variant", "seed", "m", "source", "source_vertex",
             "n_vertices", "h_mean", "D", "MAE", "Linf", "mae_rel", "linf_rel", "tol_mae_rel",
             "tol_linf_rel", "pass_mae", "pass_linf", "passed", "required", "worst_vertex",
             "exact_at_worst", "mean_signed", "source_valence", "source_local_min_angle"]
    data = dict(
        summary=summary, baseline=baseline, fixtures=FIX, tolerances=TOL, annotations=ann,
        geo=records(geo, gcols),
        quality=records(qual, [c for c in qual.columns]),
        structure=records(struct, [c for c in struct.columns]),
        heat=records(heat, ["case", "family", "resolution", "variant", "seed", "h_mean", "mode", "lam",
                            "T", "n_steps", "amplitude", "amplitude_exact", "rate_rel_err", "tol",
                            "residual", "total_heat_drift", "passed", "vacuous"]),
        gauss=records(gauss, [c for c in gauss.columns]),
        neg=records(neg, [c for c in neg.columns]),
        perf=records(perf, [c for c in perf.columns]),
        perf_earlier=(lambda q: records(q, list(q.columns)))(rd("perf_summary_run1.csv"))
        if os.path.exists(os.path.join(RESULTS, "perf_summary_run1.csv")) else None,
        runs=records(runs, ["mesh", "ordering", "repeat", "setup_s", "query_s", "solver_peak_mb"]),
        meshes=meshes, faces=faces, fields=fields, gfields=gfields,
    )
    with open(os.path.join(HERE, "report_template.html")) as fh:
        tpl = fh.read()
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html = tpl.replace("/*__DATA__*/null", payload)
    assert "http://" not in tpl.replace("http://www.w3.org", "") and "https://" not in tpl, "no remote resources"
    with open(REPORT, "w") as fh:
        fh.write(html)
    print(f"wrote {os.path.relpath(REPORT)} ({len(html) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
