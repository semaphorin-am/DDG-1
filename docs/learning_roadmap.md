# Codebase Review & Learning Roadmap

A pedagogical review of this repository: what it already teaches well, what is
missing, and a ladder of projects for learning **programming** and **discrete
differential geometry (DDG)** from an applied, computational perspective.

---

## 1. Assessment: what this codebase is

This is a small (~1,500 line) but unusually clean implementation of the core
DDG/DEC pipeline, validated end-to-end against closed-form answers on the unit
sphere. As a learning artifact it has several real strengths:

- **Theory ↔ code traceability.** Every module docstring states the
  mathematical object it implements, and `docs/theory.md` maps each concept to
  the function that realizes it. This is the single most valuable pedagogical
  property of the repo — preserve it in everything added.
- **Validation-first culture.** Each stage is checked against an invariant or
  an analytic truth: `χ = 2`, `d² = 0` to machine precision, `tr(M) → 4π`,
  `L𝟙 = 0`, exactness of the face gradient on linear fields, heat
  conservation, `λ_ℓ = ℓ(ℓ+1)`. The tests in `tests/test_ddg.py` read like a
  syllabus of "what must be true."
- **Honest empirics.** The convergence study reports `p ≈ 0.7` rather than the
  textbook `p ≈ 1` and says so — leaving a genuine open question to
  investigate (see §5).
- **Right-sized scope.** One pipeline (mesh → operators → heat → geodesics →
  spectrum), one domain (the sphere), no framework overhead. A learner can
  read the whole thing in an afternoon.

The main pedagogical *gaps*, which the project ladder below addresses:

- **Only one topology and one mesh family.** Everything runs on near-perfect
  icospheres of a genus-0 surface. Most of DDG's interesting behavior (negative
  cotangent weights, harmonic forms, Betti numbers, mesh-quality sensitivity)
  is invisible here.
- **No curvature.** A DDG codebase without angle defects, Gauss–Bonnet, or the
  mean-curvature normal is missing the most geometrically vivid material.
- **DEC is asserted, not constructed.** `theory.md` states `L = d₀ᵀ ⋆₁ d₀`,
  but the code builds `L` directly from cotangents and never assembles it from
  the DEC operators it so carefully constructs in `topology.py`. The Hodge
  stars `⋆₁`, `⋆₂` and the codifferential never appear as matrices.
- **No software-engineering scaffolding.** No packaging, CI, type checking, or
  notebooks — all of which are themselves good learning projects (§4).

---

## 2. Project ladder — DDG concepts

Ordered roughly by difficulty. Each project states the concept, the task, and
the validation criterion (keeping the repo's "test against truth" ethos).

### Tier A — exercises inside the existing pipeline

**A1. Close the DEC loop: build `L` from `d₀` and `⋆₁`.**
Assemble `⋆₁ = diag(cotangent edge weights)` using
`cotangent_weights_per_edge`, form `d₀ᵀ ⋆₁ d₀`, and verify it equals the
output of `cotangent_laplacian` to machine precision. One of the highest
insight-per-line exercises available: it makes the "topology × metric"
factorization of the Laplacian concrete.
*Validate:* `‖d₀ᵀ⋆₁d₀ − L‖_max < 1e-12`.

**A2. Curvature module (`ddg/curvature.py`).**
(a) Gaussian curvature via angle defect `K_i = (2π − Σθ)/A_i`; (b) the
mean-curvature normal `H n = (L @ V) / (2 A)` ("the Laplacian of position is
curvature"); (c) discrete **Gauss–Bonnet**: `Σ angle defects = 2πχ` — note
this holds *exactly*, independent of mesh quality, just like `d² = 0`.
*Validate:* on the unit sphere `K → 1`, `H → 1`; Gauss–Bonnet exact at every
subdivision level.

**A3. Mass-matrix variants.**
Implement (a) the full Galerkin (linear FEM) mass matrix and (b) mixed
Voronoi vertex areas (Meyer et al. 2003), alongside the current lumped
barycentric one. Compare spectral accuracy and geodesic error across all three.
*Validate:* Galerkin `M` improves eigenvalue accuracy; all give `tr/sum → 4π`.

**A4. Time integrators for the heat equation.**
Add explicit Euler and Crank–Nicolson to `HeatSolver`. Find the explicit CFL
limit experimentally as a function of `h`; show Crank–Nicolson is second-order
in `dt` against the exact spherical heat kernel
`u(θ,t) = Σ_ℓ e^{−ℓ(ℓ+1)t} (2ℓ+1)/(4π) P_ℓ(cosθ)`.
*Validate:* measured orders ≈ 1 and 2; explicit blow-up threshold ∝ h².

**A5. Spectral fidelity beyond eigenvalues.**
The current spectral check compares eigen*values* only. Project each computed
eigenvector band onto the span of analytic spherical harmonics `Y_ℓ^m`
(sampled at vertices) and report the subspace alignment residual per band.
Teaches why individual eigenvectors in a multiplicity-`(2ℓ+1)` band are
arbitrary but the *subspace* is not.

**A6. Baseline geodesics: Dijkstra and graph-distance bias.**
Implement edge-graph Dijkstra (`scipy.sparse.csgraph`) and compare to the heat
method. Dijkstra over-estimates by a constant factor that does *not* vanish
under refinement — a memorable lesson that "finer mesh" ≠ "convergent
algorithm."
*Validate:* heat-method error → 0; Dijkstra error plateaus.

**A7. The `m`-parameter study.**
Sweep the diffusion-time multiplier `m` in `HeatMethodSolver` over
`[0.1, 100]` across subdivision levels; plot error vs `m`. Connects directly
to the open `p ≈ 0.7` question (§5).

### Tier B — new domains and meshes

**B1. The torus (`χ = 0`).**
Generate a torus mesh; verify `V − E + F = 0`; compute its spectrum and compare
to known flat-torus eigenvalues for the corresponding lattice. The heat method
now meets a **cut locus** — visualize where the distance field is least
accurate.

**B2. Mesh quality and negative cotangent weights.**
Randomly perturb icosphere vertices (tangentially, then re-project) to create
obtuse triangles. Show that cotangent weights go negative, the discrete maximum
principle fails (heat `u` dips below 0), and geodesic error degrades. Then
implement **intrinsic Delaunay edge flipping** (Bobenko–Springborn; Sharp &
Crane 2020) and show it repairs the weights without moving vertices. This is
the canonical "robustness" arc in modern DDG.

**B3. Real meshes via `trimesh`.**
Load a scanned mesh (e.g. the Stanford bunny). No analytic truth exists, so
introduce **self-convergence** (compare against the finest level) and
cross-validation against an exact polyhedral geodesic implementation. Teaches
how to do numerics when `arccos(p·q)` isn't available — the normal situation.

**B4. Surfaces with boundary.**
Cut the sphere into a hemisphere. Implement Dirichlet and Neumann conditions
for the heat equation and the Poisson solve; observe what happens to
`L𝟙 = 0` and heat conservation. The pinning trick in `geodesic._pinned_solve`
becomes unnecessary for the Dirichlet problem — understand why (the kernel
disappears with the boundary).

### Tier C — capstone projects

**C1. Full DEC and Hodge decomposition.**
Build `⋆₀, ⋆₁, ⋆₂` as diagonal matrices, the codifferential
`δ = ⋆⁻¹ dᵀ ⋆`, and the 1-form Laplacian `Δ₁ = dδ + δd`. Decompose a random
1-cochain into exact + coexact + harmonic parts. Compute Betti numbers from
`dim ker Δ₁`: 0 on the sphere, 2 on the torus — *cohomology, computed*.
*Validate:* the three parts are mutually `⋆`-orthogonal and sum to the input;
harmonic dimension matches `2g`.

**C2. Mean curvature flow and conformalized MCF.**
Semi-implicit MCF: `(M − dt·(−L)) V^{n+1} = M V^n` applied to vertex positions
(note the stiffness `L` is reassembled or frozen — compare both). Then
conformalized MCF (Kazhdan–Solomon–Ben-Chen 2012) to flow the bunny to a
sphere. Visually spectacular; numerically subtle.

**C3. The vector heat method / logarithmic map.**
Extend the heat method to parallel-transport *vectors* from a source (Sharp,
Soliman & Crane 2019): connection Laplacian, complex/2×2-block structure,
holonomy on the sphere. The natural sequel to `geodesic.py`.

**C4. Spectral shape analysis.**
Heat kernel signature and diffusion distance from the eigenpairs already
computed in `spectral.py`; verify Weyl's law `N(λ) ~ (A/4π)λ` on the sphere
(area recovered from eigenvalue counting!).

**C5. Reaction–diffusion (Turing patterns) on the sphere.**
Couple two heat solvers with a Gray–Scott reaction term. Everything needed
already exists (`M`, `L`, prefactored implicit steps); the payoff is striking
pattern formation on a curved surface, and a lesson in operator splitting.

---

## 3. Project ladder — learning to program

These use the repo as a vehicle for software craft; each is independently
mergeable.

**P1. Packaging.** Add `pyproject.toml`, make `pip install -e .` work, remove
the `sys.path.insert` hack in `tests/test_ddg.py`. Learn: modern Python
packaging, editable installs.

**P2. Real pytest.** Convert the hand-rolled `__main__` runner to idiomatic
pytest: `@pytest.mark.parametrize` over subdivision levels, fixtures for
shared meshes (they're expensive — learn `scope="module"`), `pytest.approx`.
Add **property-based tests** with `hypothesis` (e.g. `d² = 0` and Gauss–Bonnet
for randomly generated valid triangulations — they hold for *any* mesh).

**P3. Continuous integration.** GitHub Actions workflow running pytest + a
small benchmark smoke test on push. Learn: CI, caching dependencies, keeping
the suite fast enough to run on every commit.

**P4. Static quality.** Fix the typing bug in `TriMesh`
(`edges: np.ndarray = field(default=None)` — the annotation says it can never
be `None`); run `mypy` and `ruff` clean. Learn: how type checkers catch real
API ambiguities (`Optional` vs sentinel).

**P5. Vectorize `build_d1`.** It is the only Python-level per-face loop in the
operator stack. Vectorize it with NumPy (sort each face's vertex pairs, look
up edge ids via `np.searchsorted` on the sorted edge array), benchmark before
and after at subdivision 7–8. Learn: profiling (`cProfile`,
`line_profiler`), the cost model of NumPy vs interpreter loops.

**P6. Refactor: a cached geometry context.** `face_areas`, `face_normals`,
`interior_angles` are recomputed in several places, and
`geodesic._pinned_solve` duplicates the pinning logic inside
`HeatMethodSolver`. Design a small `Geometry` object that lazily computes and
caches derived quantities. Learn: API design, `functools.cached_property`,
when caching is worth its complexity.

**P7. Honest docstrings as a code-review exercise.** `heat.py` claims the
system is prefactored "with a sparse Cholesky" — `scipy.sparse.linalg.
factorized` actually performs SuperLU LU. Fix the docstring, then go further:
benchmark SuperLU vs CHOLMOD (`scikit-sparse`) vs conjugate gradients with an
AMG preconditioner (`pyamg`) on subdivision 7–8. Learn: reading library
internals, sparse-solver tradeoffs.

**P8. Notebooks and interactive visualization.** One Jupyter notebook per
pipeline phase, mirroring `theory.md`; add an optional `polyscope` viewer for
eigenfunctions and distance fields. Learn: literate programming, separating
library code from presentation.

**P9. Reproducibility hardening.** Pin dependency versions, seed all RNGs,
make `run_benchmark.py` emit a manifest (git SHA, versions, parameters) next
to its results. Learn: why "it ran on my machine in June" is not a result.

---

## 4. Open questions already in the data (mini-research projects)

These are real, currently-unexplained observations in `results/` — ideal for
teaching the *investigative* side of computational science:

1. **Why `p ≈ 0.7` and not `1.0`?** The heat method is nominally first-order
   in the smoothed distance. Hypotheses to test: (a) the error is dominated by
   a neighborhood of the source and/or the antipode (cut locus) — plot
   `|φ − d|` on the sphere per level; (b) the `t = m·h²` coupling is
   suboptimal — re-fit `p` for the best `m` per level (project A7); (c) the
   MAE is contaminated by a global shift — compare after removing the mean
   error. Each hypothesis is one afternoon and a plot.
2. **Eigenvalue error vs band.** From `results/spectrum.csv`, fit how the
   relative error of `λ` grows with `ℓ` at fixed `h`. The classical estimate
   is `O(λ h²)`-ish — does it hold? This motivates project A3 (consistent
   mass matrix).
3. **Where does the area deficit live?** `tr(M) − 4π` shrinks like `h²`.
   Decompose it per-face: is it uniform, or concentrated at the 12 original
   icosahedron vertices (the only valence-5 vertices)? Connects mesh
   combinatorics to metric error.

---

## 5. Suggested sequencing

For a learner doing both tracks simultaneously, a natural interleaving:

| Week | DDG | Programming |
|---|---|---|
| 1 | Read `theory.md` + code; rerun benchmark | P1 packaging, P2 pytest |
| 2 | A1 (DEC Laplacian), A2 (curvature, Gauss–Bonnet) | P4 typing/lint |
| 3 | A4 (integrators), A6 (Dijkstra baseline) | P3 CI, P5 vectorize `d₁` |
| 4 | A7 + open question 1 (`p ≈ 0.7`) | P6 geometry cache |
| 5 | B1 (torus) or B2 (mesh quality / intrinsic Delaunay) | P8 notebooks |
| 6+ | One capstone: C1 (Hodge), C2 (MCF), C3 (vector heat), or C5 (Turing) | P7, P9 |

---

## 6. Companion resources

- **Keenan Crane, *Discrete Differential Geometry: An Applied Introduction*** —
  free CMU course notes; this repo is essentially a Python lab manual for
  chapters on simplicial complexes, exterior calculus, the Laplacian, and
  geodesics. The course's coding assignments (written for JS/C++) port
  naturally to this codebase.
- Crane, Weischedel & Wardetzky, *Geodesics in Heat* (TOG 2013) and the
  follow-up survey (CACM 2017) — the paper behind `geodesic.py`.
- Botsch et al., *Polygon Mesh Processing* — mesh data structures and the
  engineering side (relevant to projects B2, B3).
- Sharp & Crane, *A Laplacian for Nonmanifold Triangle Meshes* /
  *navigating intrinsic triangulations* — background for B2.
- Desbrun, Hirani, Leok & Marsden, *Discrete Exterior Calculus* — the formal
  backbone for C1.
