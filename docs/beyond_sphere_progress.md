# Progress log: geometry beyond the sphere

Short, append-only log. Newest entries at the bottom.

## 1. Starting point

- Default branch: `claude/ddg-sphere-benchmark-zZ5Qf`
- Starting commit: `018fdb4a00d155006ed660c8da8720e4c37f74d1`
- Baseline tests: `python -m pytest -q tests` -> 9 passed.
- Baseline accuracy/timing on this machine: `results/beyond_sphere/baseline.json`
  (regular icosphere, source = vertex 0, m = 1). MAE matches the committed
  `results/convergence.csv` to 11 digits.

## 2. References, fixtures, tolerances (committed before any solver change)

- `benchmarks/references.py`: closed-form answers, no `ddg` import.
- `benchmarks/fixtures.json`: resolutions, seeds, sources, modes, negative controls,
  performance mesh list.
- `benchmarks/tolerances.json`: acceptance thresholds.
- How the tolerances were chosen: geodesic thresholds are the baseline
  icosphere fit times 1.5; heat thresholds are the textbook implicit-Euler +
  second-order-space estimate. No perturbed-sphere or torus result existed
  when they were written.
- The "strong" sphere perturbation is declared a stress case up front: it is
  scored against the same tolerances but is not required to pass.
