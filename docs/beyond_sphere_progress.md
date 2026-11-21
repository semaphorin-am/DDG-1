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

## 3. Meshes and intrinsic geometry

- `perturbed_icosphere`, `flat_torus` (one vertex per periodic class, per-face
  unwrapped corners), `open_rectangle` (negative control). Geometry operators
  now read per-face corners, so seam triangles get their true flat shape.
- Original 9 tests still pass unchanged.

## 4. First full accuracy sweep: checks that fail, and why

Tolerances were not touched. Outcomes of the first sweep (`benchmarks/run.py`):

- **Regular torus, geodesic, default m = 1: fails.** L-inf fails at every
  resolution; MAE stops improving at about 0.6 % of the diameter from 96x64
  on, so it also fails MAE from 192x128. The result is identical for all five
  sources (translation invariance holds to rounding), so this is not a seam
  bug. It is the heat method's known smoothing bias at t = h^2, larger on the
  flat torus because the cut locus is a pair of lines rather than one point.
  m = 4 passes MAE at every resolution but still not L-inf. Reported as a
  finding; the sphere-derived tolerance stays.
- **Regular sphere L7 and mild L6: partial fail** for the same reason (MAE
  plateau), first visible one level beyond the original benchmark.
- **Heat-mode "total heat drift" failed on jittered tori - this check was
  wrong.** My first implementation divided the drift by the mode's own
  integral, which is ~0 for an oscillating mode on a jittered mesh, turning
  rounding (1e-17 absolute) into 1e-8 "relative". Fixed the normalisation to
  int |u0| dA. The threshold (1e-10) is unchanged; decay-rate and residual
  checks were passing throughout.

## 5. Bottleneck and the two approaches tried

- Measured bottleneck: sparse LU factorisation in `HeatMethodSolver` setup
  (two factorizations; >90 % of setup time at 40k+ vertices).
- Approach 1: SuperLU built-in orderings / symmetric mode / RCM. No reliable
  gain; symmetric mode was ~100x slower on the sphere, RCM increased fill.
- Approach 2: geometric nested-dissection ordering (`ddg/ordering.py`),
  opt-in via `ordering="nested"`. Fields agree with the default to ~1e-13.
  Kept. No third approach attempted.
