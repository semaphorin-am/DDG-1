# Learning Roadmap

This codebase is a small, working scientific program (~1,500 lines). It builds
a triangle mesh of a sphere, assembles some matrices, solves equations on the
surface, and checks every answer against known exact results.

That makes it an ideal vehicle for learning **how to program and how to design
and build software**. It is small enough to read in full, real enough to have
genuine design problems, and — because the right answers are known — you can
always tell whether a change broke something.

The plan below is one ordered sequence, not a menu. Each step teaches one
fundamental skill, and each step builds on the previous one. The mathematics
stays fixed; the *software* is what improves.

---

## The principle

> Don't add features. Take a working program and make it *good*, one
> discipline at a time. Then add one feature using everything you learned.

---

## Step 1 — Read it and run it

**Skill: reading code you didn't write.**

- Run the tests (`python tests/test_ddg.py`) and the benchmark
  (`python scripts/run_benchmark.py --all`).
- Read the modules in pipeline order: `mesh.py` → `topology.py` →
  `geometry.py` → `heat.py` → `geodesic.py` → `spectral.py` → `analysis.py`.
  Read `docs/theory.md` alongside.
- Write down, in your own words, what each module's inputs and outputs are.
  If you can't, read again. You cannot improve code you can't summarize.

**Done when:** you can sketch the data flow (which arrays go where) from
memory.

---

## Step 2 — Make it a real package

**Skill: project structure and packaging.**

Right now the tests hack the import path
(`sys.path.insert(...)` at the top of `tests/test_ddg.py`). That's a smell:
the project isn't installable.

- Add a `pyproject.toml` so `pip install -e .` works.
- Delete the `sys.path` hack; the tests should just `import ddg`.

**Done when:** a fresh clone + `pip install -e .` + `pytest` passes with no
path tricks.

**Lesson:** every serious project starts with making it trivially runnable by
someone else.

---

## Step 3 — Real tests

**Skill: automated testing as a design tool.**

The tests exist but are run by a hand-rolled loop at the bottom of the file.
Convert them to idiomatic `pytest`:

- Use `@pytest.mark.parametrize` for the loops over subdivision levels.
- Meshes are expensive to build, and several tests rebuild the same one.
  Use a fixture with `scope="module"` to build each mesh once and share it.
- Use `pytest.approx` instead of manual tolerance comparisons.

**Done when:** `pytest` runs everything, the suite is noticeably faster than
before (because meshes are shared), and the hand-rolled runner is deleted.

**Lesson:** tests are code too — they deserve structure, and slow tests stop
getting run.

---

## Step 4 — Types and linting

**Skill: letting tools find your bugs.**

- Run `ruff` and `mypy` on the project and fix what they find.
- There is at least one real typing bug waiting: in `mesh.py`, `TriMesh`
  declares `edges: np.ndarray = field(default=None)`. The annotation says
  "always an array", the default says "might be None". A type checker flags
  this immediately; fix it properly (with `Optional` and a clear contract,
  or by restructuring so the field is never None).

**Done when:** `ruff` and `mypy` both pass cleanly.

**Lesson:** type annotations are claims about your code. Tools can check
claims; comments they can't.

---

## Step 5 — Continuous integration

**Skill: automation.**

Add a GitHub Actions workflow that runs the linter, type checker, and test
suite on every push.

**Done when:** a push to a branch shows a green check, and an intentionally
broken commit shows a red one.

**Lesson:** quality checks that depend on you remembering to run them will
eventually not be run.

---

## Step 6 — Refactor: remove duplication and recomputation

**Skill: software design — the heart of this roadmap.**

Now that tests and CI protect you, you can reshape the code safely. Two
concrete design problems already exist:

1. **Recomputation.** `face_areas`, `face_normals`, and `interior_angles` are
   computed from scratch in several places (look at how `HeatMethodSolver`
   passes them around by hand to avoid this).
2. **Duplication.** `geodesic.py` contains `_pinned_solve`, and
   `HeatMethodSolver` contains a second copy of the same pinning logic.

Design a small `Geometry` class that owns a mesh and lazily computes and
caches derived quantities (`functools.cached_property` is the natural tool).
Migrate the callers. Unify the pinned solve into one implementation.

**Done when:** each derived quantity is computed in exactly one place, the
duplication is gone, and all tests still pass.

**Lesson:** this is what design *is* — deciding what owns what, and making
each fact live in one place. The tests from Step 3 are what make the
refactor safe; that's why testing came first.

---

## Step 7 — Measure, then optimize

**Skill: performance work driven by evidence, not guesses.**

- Profile the benchmark (`cProfile`) at a high subdivision level. Find where
  the time actually goes.
- One known hotspot: `build_d1` in `topology.py` is the only per-face Python
  loop in the whole operator stack — everything else uses NumPy arrays.
  Vectorize it and measure the speedup.

**Done when:** you have before/after timings written down, and the tests
(including `d² = 0`, which would catch any sign mistake) still pass.

**Lesson:** never optimize without a profile, and never optimize without a
test that proves you didn't change the answer.

---

## Step 8 — Make the documentation true

**Skill: documentation hygiene.**

Docs are part of the software. There is a real example here: the docstring in
`heat.py` says the solver prefactors the matrix with "a sparse Cholesky", but
`scipy.sparse.linalg.factorized` actually performs an LU factorization
(SuperLU). Small, but false.

- Fix it, and audit the other docstrings against the code while you're there.

**Lesson:** wrong documentation is worse than none — readers trust it.

---

## Step 9 — Capstone: build one feature end-to-end

**Skill: putting it all together.**

Add **one** new module, `ddg/curvature.py`, using every discipline from the
previous steps. It's a natural fit (curvature is the one obvious thing this
geometry library lacks) and it's small:

- Gaussian curvature at each vertex via the angle defect
  (`2π − sum of incident angles`, divided by vertex area).
- A test against exact truth, in the existing style: on the unit sphere the
  curvature should be ≈ 1 everywhere, and the **sum of angle defects must
  equal `2π × Euler characteristic` exactly** — for any mesh, like the
  existing `d² = 0` check.

Work the full loop: design the function signatures first, write the test,
implement, type-check, lint, document, let CI verify.

**Done when:** the feature ships with tests, types, and docs — and you didn't
have to be reminded of any of those, because by now it's just how you work.

---

## Summary of what each step teaches

| Step | Fundamental |
|---|---|
| 1 | Reading and summarizing unfamiliar code |
| 2 | Project structure and packaging |
| 3 | Testing as a first-class activity |
| 4 | Types and static analysis |
| 5 | Automation (CI) |
| 6 | Design: ownership, caching, removing duplication |
| 7 | Evidence-driven performance work |
| 8 | Documentation that stays true |
| 9 | Shipping a feature with all of the above |

The mathematics (discrete differential geometry) stays in the background as
the thing the software *does* — and the existing exact checks (`χ = 2`,
`d² = 0`, spectrum `= ℓ(ℓ+1)`) act as an incorruptible safety net for every
refactor. If you later want to go deeper into the geometry itself, Keenan
Crane's free course notes *Discrete Differential Geometry: An Applied
Introduction* pair naturally with this code.
