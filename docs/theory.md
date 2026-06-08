# Discrete Differential Geometry on the Sphere — Theory Notes

These notes document the mathematics behind the benchmark and map each concept
to the code that implements it. They are meant to be read alongside the module
docstrings in `ddg/`.

---

## 1. Simplicial complexes and triangle meshes

A triangle mesh is a **2-dimensional simplicial complex**: a set of vertices
(0-simplices), edges (1-simplices), and triangular faces (2-simplices) glued
along shared sub-simplices. We store it with three arrays
(`ddg/mesh.py`, `TriMesh`):

- `vertices` — `(|V|, 3)` coordinates,
- `faces` — `(|F|, 3)` vertex indices, **consistently oriented** (every face
  wound counter-clockwise as seen from outside, i.e. outward normals),
- `edges` — `(|E|, 2)` unique undirected edges with the convention `i < j`.

For any closed surface the **Euler characteristic** is a topological invariant

$$\chi = |V| - |E| + |F|,$$

with $\chi = 2$ for the sphere. A closed triangulation also satisfies
$3|F| = 2|E|$ (every triangle has three edges, every edge borders two
triangles). Both are checked in `tests/test_ddg.py`.

### Icosphere

We generate meshes by recursively subdividing a regular **icosahedron** and
projecting new vertices onto the sphere (`icosphere`). Each 1→4 subdivision
multiplies the face count by 4, giving the family

$$|V| = 10 \cdot 4^{n} + 2, \qquad |F| = 20 \cdot 4^{n}.$$

Icosphere triangles are nearly equilateral and uniform, which keeps cotangent
weights well-behaved and makes the convergence study clean.

---

## 2. Incidence matrices and the discrete exterior derivative

Functions sampled on $k$-simplices are **discrete $k$-forms** (cochains). The
**discrete exterior derivative** $d_k$ is the signed boundary/incidence operator
(`ddg/topology.py`):

- $d_0 \in \mathbb{R}^{|E| \times |V|}$ maps vertex functions to edge functions.
  Row for oriented edge $i \to j$ has $-1$ at $i$ and $+1$ at $j$, so
  $(d_0 u)_{ij} = u_j - u_i$ — a discrete gradient along edges.
- $d_1 \in \mathbb{R}^{|F| \times |E|}$ maps edge functions to face functions.
  Each face contributes $\pm 1$ per boundary edge depending on whether the
  face's traversal agrees with the stored edge orientation — a discrete curl.

The cornerstone identity of any cochain complex is $d^2 = 0$, the discrete
analogue of "the boundary of a boundary is empty":

$$d_1 d_0 = 0.$$

`verify_d_squared` checks this is exactly zero (to machine precision) — a strong
correctness test on orientation bookkeeping.

---

## 3. Hodge stars, mass matrix, and metric

The exterior derivative is purely topological. **Geometry** enters through the
**Hodge star** operators, which depend on the metric (edge lengths and areas).
In Discrete Exterior Calculus the Hodge stars are diagonal matrices relating
primal $k$-forms to dual $(n-k)$-forms (`ddg/geometry.py`):

- $\star_0$ — the **mass matrix** $M$, a diagonal matrix of **vertex areas**.
  We use the *lumped barycentric* area: each vertex receives one third of the
  area of every incident triangle, so $\operatorname{tr}(M)$ equals the total
  surface area. $M$ is the discrete integration measure $\int_{S^2} u \, dA
  \approx \mathbf{1}^\top M u$.
- $\star_1$ — the **cotangent weights** on edges,
  $w_{ij} = \tfrac12(\cot\alpha_{ij} + \cot\beta_{ij})$, where $\alpha,\beta$ are
  the two angles opposite edge $ij$.

As the mesh refines, $\operatorname{tr}(M) \to 4\pi$ and individual face areas
$\to 0$.

---

## 4. The cotangent Laplacian

The **Laplace–Beltrami operator** generalizes $\Delta = \nabla\cdot\nabla$ to a
curved surface. Its standard discretization is the **cotangent Laplacian**. We
assemble the symmetric **stiffness matrix** (`cotangent_laplacian`)

$$L_{ij} = -w_{ij}\ (i\neq j), \qquad L_{ii} = \sum_{j} w_{ij},$$

equivalently $L = d_0^\top \star_1 d_0$. Key properties (checked in tests):

- **symmetric** ($L = L^\top$),
- **positive semidefinite** ($x^\top L x = \sum_{ij} w_{ij}(x_i-x_j)^2 \ge 0$;
  this is the discrete Dirichlet energy),
- **constant nullspace** ($L\mathbf{1} = 0$),
- **sparse** (one nonzero per incident edge).

The Laplace–Beltrami *operator* is then $\Delta = -M^{-1} L$ (negative
semidefinite). The sign matters: with this convention the heat equation and the
Poisson recovery in the heat method use $L$ directly, but the operator equation
$\Delta\phi = f$ becomes $L\phi = -Mf$.

### Spectrum

The generalized eigenproblem $L\phi_i = \lambda_i M\phi_i$ approximates the
continuous Laplace–Beltrami spectrum. On the unit sphere this is known exactly:

$$\lambda_\ell = \ell(\ell+1), \quad \ell = 0,1,2,\dots, \quad
\text{multiplicity } 2\ell+1,$$

with eigenfunctions the **spherical harmonics** $Y_\ell^m$
(`ddg/spectral.py`). Watching the discrete eigenvalues approach
$0, 2, 2, 2, 6, 6, \dots$ and the eigenfunctions look like spherical harmonics is
the spectral-geometry validation.

---

## 5. Heat equation

We solve $\partial_t u = \Delta u$ from a point heat source. Discretizing time
with **backward (implicit) Euler** gives the unconditionally stable step
(`ddg/heat.py`)

$$(M + \Delta t\, L)\, u^{n+1} = M u^{n}.$$

The system matrix is symmetric positive-definite and constant in time, so we
prefactor it once and reuse the factorization. On a closed surface the total
heat $\mathbf{1}^\top M u$ is conserved — a useful diagnostic.

---

## 6. Geodesic distance: the Heat Method

The **heat method** (Crane, Weischedel & Wardetzky, *Geodesics in Heat*, 2013)
computes geodesic distance with three linear solves (`ddg/geodesic.py`). It
rests on **Varadhan's formula**: for small time $t$,

$$d(x,y) = \lim_{t\to 0} \sqrt{-4t \log k_t(x,y)},$$

where $k_t$ is the heat kernel. Rather than the unstable logarithm, the method
keeps only the *direction* of $\nabla u$, which already aligns with geodesics:

1. **Heat flow.** Solve $(M + t L) u = \delta_{\text{source}}$ with
   $t = m\,h^2$ ($h$ = mean edge length, $m \approx 1$).
2. **Normalized gradient.** Compute the per-face gradient and reverse/normalize:
   $X = -\nabla u / \lVert \nabla u\rVert$. This unit field points along
   geodesics away from the source. The per-triangle gradient is
   $$\nabla u = \frac{1}{2A}\sum_{m} u_m\,(N \times e_m),$$
   with $e_m$ the edge opposite vertex $m$ (exact for linear fields — tested).
3. **Poisson recovery.** Find the scalar field whose gradient best matches $X$
   by solving $\Delta\phi = \nabla\cdot X$, i.e. $L\phi = -\operatorname{div}X$,
   then shift so $\phi(\text{source}) = 0$. The integrated divergence at a
   vertex is
   $$(\nabla\cdot X)_i = \tfrac12\sum_{f \ni i}
     \cot\theta_1\,(e_1\cdot X_f) + \cot\theta_2\,(e_2\cdot X_f).$$

Because only the right-hand side depends on the source, both matrices are
prefactored once and reused for any source.

---

## 7. Analytical validation and convergence

On the unit sphere the exact geodesic distance between $p$ and $q$ is the
central angle (`ddg/analysis.py`)

$$d(p,q) = \arccos(p\cdot q).$$

We report three error norms over all vertices:

$$\mathrm{MAE} = \tfrac1N\sum_i |d_i^{\text{num}} - d_i^{\text{exact}}|, \quad
L_\infty = \max_i |\cdot|, \quad
\mathrm{RMSE} = \sqrt{\tfrac1N\sum_i (\cdot)^2}.$$

Fitting $\text{error} \sim C h^{p}$ in log–log space (`convergence_rate`) gives
the empirical **convergence rate** $p$ as the mesh is refined. The heat method
is first-order accurate ($p \approx 1$) in the smoothed distance, which the
study reproduces.

---

## 8. Performance and scaling

All operators are sparse with $O(|V|)$ nonzeros. Assembly is $O(|V|)$; the
sparse Cholesky/LU solves and the Lanczos eigensolver scale close to linearly on
these well-shaped meshes. `scripts/run_benchmark.py --phase performance`
measures assembly, solve, and eigensolve time plus peak memory, and fits the
scaling — quantifying the practical hardware limits (success criterion 4).
