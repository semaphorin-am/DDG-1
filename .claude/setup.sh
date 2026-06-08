#!/usr/bin/env bash
# SessionStart hook: ensure the DDG benchmark can run in a fresh web session.
# Installs the Python dependencies if the core stack is missing. Idempotent and
# quiet on the happy path.
set -euo pipefail
cd "$(dirname "$0")/.."

if python3 -c "import numpy, scipy, matplotlib, pandas, trimesh, tabulate" \
    >/dev/null 2>&1; then
  echo "DDG benchmark: dependencies already present."
else
  echo "DDG benchmark: installing dependencies from requirements.txt ..."
  pip install --quiet -r requirements.txt
  echo "DDG benchmark: dependencies installed."
fi
