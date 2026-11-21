#!/usr/bin/env bash
# One command from a fresh checkout: install, run every check, regenerate the
# raw results, static figures and the standalone HTML report.
#   ./reproduce.sh            full run (about 15-25 minutes, needs ~4 GB RAM)
#   SKIP_PERF=1 ./reproduce.sh   reuse the committed performance measurements
set -euo pipefail
cd "$(dirname "$0")"

PY=${PYTHON:-python3}
if [ -z "${NO_VENV:-}" ]; then
  [ -d .venv ] || "$PY" -m venv .venv
  # shellcheck disable=SC1091
  . .venv/bin/activate
  PY=python
fi
"$PY" -m pip install -q -r requirements.txt pytest

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1

echo "== 1/5 tests (original sphere suite + beyond-sphere suite)"
"$PY" -m pytest -q tests
echo "== 2/5 accuracy sweep -> results/beyond_sphere/"
"$PY" benchmarks/run.py
echo "== 3/5 performance (baseline vs nested dissection)"
if [ -z "${SKIP_PERF:-}" ]; then "$PY" benchmarks/perf.py; else echo "skipped (SKIP_PERF set)"; fi
echo "== 4/5 static figures -> figures/beyond_sphere/"
"$PY" benchmarks/figures.py
echo "== 5/5 report -> docs/beyond_sphere_report.html"
"$PY" benchmarks/report.py
echo "Done. Open docs/beyond_sphere_report.html in a browser (no server needed)."
