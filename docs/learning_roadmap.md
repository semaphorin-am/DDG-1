# Learning Roadmap

**Who this is for:** a curious student (around 8th grade) with little or no
coding experience, who wants to learn three things at once:

1. how to write programs and build real software,
2. what **matrices** and **calculus** actually are — not from a textbook, but
   by computing with them, and
3. the geometry of curved surfaces (the grown-up name for this subject is
   *discrete differential geometry*, which means: calculus on shapes made of
   triangles).

**A secret about this subject:** calculus is usually taught with limits and
infinity, and matrices as walls of symbols. But on a shape made of triangles,
both become completely concrete: a *derivative* is just a difference between
neighbors, an *integral* is just a weighted sum, a *matrix* is just a machine
that does many of those at once, and a *differential equation* is just an
update rule you run in a loop. Everything in this roadmap is real calculus
and real linear algebra — you'll simply meet it in a form you can compute by
hand and check by program. When you meet the official versions in high school
and college, they will feel like old friends.

**What this repository is:** a small, finished program written by
professionals. It builds a ball out of triangles, measures it, spreads heat
across it, and computes distances along its curved surface. You won't
understand it at first — it's the *destination*. By the end you'll read it,
use it, test it, and build your own application on top of it.

**The one rule:** every stage ends with something you made that works — a
number, a picture, or a program someone else can run. Each stage has a
**Done when** line. Don't move on until you hit it.

You'll need a helper (parent, teacher, or older student) for about an hour at
the start, to install Python and get this folder onto your computer. After
that, the roadmap is yours.

---

# Part I — Learning to code, using shapes

## Stage 1 — Hello, Python

**Skill: giving a computer instructions.**

- With your helper, install Python and open a terminal in this folder.
- Start Python by typing `python`. Try:

  ```python
  print("hello")
  2 + 2
  10 * 4 ** 3 + 2
  ```

- Learn the four things almost every program is made of: **variables**
  (naming a value), **lists** (a row of values), **loops** (do this for each
  item), and **if** (do this only when something is true). Any beginner
  Python tutorial covers these; spend a few days playing.

**Done when:** you can write a loop that prints the numbers 1 to 100, and
explain each line to your helper.

---

## Stage 2 — A shape is just lists

**Skill: representing real things as data.
Math: Euler's formula — your first theorem.**

Here is a complete description of a cube, as Python lists:

```python
corners = [
    (0,0,0), (1,0,0), (1,1,0), (0,1,0),   # bottom four corners
    (0,0,1), (1,0,1), (1,1,1), (0,1,1),   # top four corners
]
# each face names its corners by position in the list above (counting from 0)
faces = [
    (0,1,2,3), (4,5,6,7),                  # bottom, top
    (0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7),  # the four sides
]
```

That's the big secret of computer graphics and of this whole repository:
**a shape is a list of corner positions plus a list of which corners make
each face.** Nothing more.

- By hand, count the cube's corners (V), edges (E), and faces (F). Compute
  **V − E + F**. Do the same for a triangular pyramid (tetrahedron) and an
  octahedron. Same answer every time: this is **Euler's formula** — for any
  shape that could be inflated into a ball, V − E + F = 2.
- Now write a program: given the `faces` list, count the edges with a loop
  (careful — each edge is shared by two faces, don't count it twice), and
  check V − E + F = 2 *by program*.

**Done when:** your program prints `V - E + F = 2` for the cube, and you
checked the tetrahedron by hand.

---

## Stage 3 — Using someone else's code

**Skill: libraries — building on other people's work.**

Real programmers rarely start from nothing. This repository is a **library**:
a box of tools you can import.

```python
from ddg import icosphere

mesh = icosphere(2)          # a ball made of triangles, detail level 2
print(mesh.n_vertices)       # how many corners?
print(mesh.n_faces)          # how many triangles?
print(mesh.euler_characteristic())   # V - E + F ... it knows the formula!
```

- Try detail levels 0 through 6 in a loop. Confirm Euler's formula gives 2
  every time, even with 81,920 triangles.
- `mesh.vertices` is the same idea as your cube's `corners` list, just
  longer, and every corner sits at distance exactly 1 from the center.
  Check it with the Pythagorean theorem: for a corner (x, y, z), is
  x² + y² + z² equal to 1?

**Done when:** you've verified, with your own loop, that V − E + F = 2 at
every detail level, and that the corners really lie on the ball.

---

## Stage 4 — Pictures

**Skill: making the computer draw.**

- Plot the mesh's corners as dots in 3D (`matplotlib` does this in about
  five lines — your helper or an online example can show you).
- Color each dot by its height (the z value). Now you can *see* the ball.
- Browse the `figures/` folder — pictures the professional program made. By
  the end of this roadmap you'll know what every one means, and you'll have
  made your own versions.

**Done when:** you have a picture file you made yourself showing the
triangle ball.

---

## Stage 5 — Functions, and your first vector

**Skill: functions — naming a recipe so you can reuse it.
Math: a *vector* is just a list of numbers.**

```python
def euler_number(mesh):
    return mesh.n_vertices - mesh.n_faces  # <- deliberately wrong! fix it.
```

- Write three working functions: `euler_number(mesh)`,
  `longest_edge(mesh)`, and `total_area(mesh)` (the library helps:
  `from ddg import face_areas`, then add the results up).
- Put them in your own file `mytools.py` and import them. You now have a
  library of your own.

And learn one piece of vocabulary that unlocks everything later: a list of
numbers is called a **vector**. A corner's position (x, y, z) is a vector
with 3 entries. If you assign a temperature to *every* corner of a mesh with
2,562 corners, that's a vector with 2,562 entries. The library uses `numpy`,
which does arithmetic on whole vectors at once:

```python
import numpy as np
u = np.array([1.0, 2.0, 3.0])
print(u * 10)        # multiplies every entry: [10. 20. 30.]
print(u.sum())       # adds them all: 6.0
```

**Done when:** `from mytools import total_area` works, and
`total_area(icosphere(5))` prints a number close to **12.566...** (remember
that number — Stage 6 explains it).

---

# Part II — Calculus and matrices, made of triangles

## Stage 6 — Integration: adding up little pieces

**Math: what an integral is.**

The true surface area of a ball of radius 1 is 4π ≈ **12.566**. Your
`total_area` computed it by **adding up the areas of thousands of tiny flat
triangles**. That operation — chop a curved thing into tiny pieces, add the
pieces — is called **integration**, and it is half of calculus. You have
already done it.

Two experiments:

1. **Convergence.** The flat triangles cut corners, so the total is a little
   *less* than 4π. Compute the error at detail levels 2–6 and make a table.
   The error shrinks about **4×** per level: the sum of flat pieces
   *converges* to the curved truth. (In school calculus, "the limit of finer
   and finer sums" — you just *watched* that limit happen.)
2. **Who owns the area?** The library can split the surface among the
   corners — `from ddg import vertex_areas` gives each corner the patch it
   "owns" (a third of every triangle touching it). Check the patches add up
   to the same total. Now you can integrate anything: to add up a
   temperature vector `u` over the surface, weight each corner by its patch:
   `(vertex_areas(mesh) * u).sum()`. **Integral = weighted sum.** Keep this
   in your `mytools.py` as `integrate(mesh, u)`.

**Done when:** you have the error table, and your `integrate` of the
all-ones vector (`np.ones(mesh.n_vertices)`) returns the surface area —
explain to your helper why it must.

---

## Stage 7 — Differentiation: differences between neighbors. And what a matrix is

**Math: what a derivative is, and what a matrix is. The most important stage
in this roadmap.**

The other half of calculus asks: **how fast does something change from place
to place?** On a mesh, the honest answer is a subtraction: if corner 5 has
temperature 40 and its neighbor corner 2 has temperature 10, the change along
that edge is 30. That's it. **A derivative is a difference between
neighbors.** (School calculus shrinks the neighbors infinitely close; the
idea is the same.)

A mesh has thousands of edges, so you want all the differences at once. Here
is the machine that does it, for a single triangle with corners 0, 1, 2 and
edges (0→1), (1→2), (0→2):

```
            corner:  0   1   2
edge 0→1  [ -1   1   0 ]
edge 1→2  [  0  -1   1 ]
edge 0→2  [ -1   0   1 ]
```

This table of −1s and +1s is a **matrix**. To apply it to a temperature
vector, take each row, multiply it entry-by-entry against the vector, and
add. Do it by hand with temperatures u = (10, 20, 40):

- row 1: −1·10 + 1·20 + 0·40 = **10**  (change along edge 0→1)
- row 2: 0·10 − 1·20 + 1·40 = **20**  (change along edge 1→2)
- row 3: −1·10 + 0·20 + 1·40 = **30**  (change along edge 0→2)

One matrix multiplication = every difference on the shape, simultaneously.
**A matrix is a machine that eats a vector and produces a vector**, and
matrix multiplication is just "rows times entries, then add" — which you have
now done by hand.

The library builds this exact machine for any mesh:

```python
from ddg import icosphere, build_d0
import numpy as np

mesh = icosphere(3)
D = build_d0(mesh)         # the difference machine (one row per edge)
u = mesh.vertices[:, 2]    # temperature = height of each corner
changes = D @ u            # @ means "apply the matrix" — all differences at once
```

- Check the machine against your hand: `print(D.toarray())` for
  `icosphere(0)` is too big to read, so build a tiny mesh instead — or just
  verify three entries of `D @ u` with your own subtraction.
- One beautiful check: if `u` is the same number everywhere
  (`np.ones(...)`), every difference must be zero. Confirm `D @ u` is all
  zeros. A constant has zero derivative — you've verified a calculus theorem
  with one line.

**Done when:** you've multiplied the 3×3 matrix by hand and gotten
(10, 20, 30), the library's `D @ u` agrees with subtractions you did
yourself, and `D @ ones` is zero.

---

## Stage 8 — Curvature: the missing angle

**Math: what "curved" means — and a quantity that's exact, not approximate.**

Do this first *without* a computer:

- On flat paper, the angles around any point add to exactly **360°**.
- At one corner of a cube, three squares meet: 3 × 90° = **270°**. There are
  90° *missing* — that's why the corner pokes out. The missing angle is
  called the **angle defect**, and it is curvature, concentrated at a point.
- The cube: 8 corners × 90° = **720°** missing in total. The tetrahedron:
  4 corners × 180° = **720°** again. Coincidence?

It is not. **For every shape that can be inflated into a ball, the missing
angles always total exactly 720°** — the **Gauss–Bonnet theorem**. Verify it
on the triangle ball:

- `from ddg.geometry import interior_angles` gives every triangle's three
  angles (in radians: a full turn is 2π, so 720° = 4π ≈ 12.566 — that number
  yet again!). For each corner, add the angles touching it, subtract from
  2π, and total over all corners. Use a loop, or be fancy with
  `np.add.at` like the professionals do in `geometry.py`.
- Try every detail level: **exactly 4π every time**, not approximately.
  Compare with Stage 6, where the area was only ever *close* to 4π. Some
  things converge; some are perfectly true on every mesh. Knowing which is
  which is real mathematical taste. (In calculus language: the total defect
  is the *integral of curvature* — Stage 6's weighted sums and this stage
  are the same idea.)

**Done when:** paper gives 720° for the cube, your program gives 4π for
every icosphere, and you can explain the difference between exact and
approximate to your helper.

---

## Stage 9 — The heat equation: calculus in motion

**Math: what a differential equation is. Skill: simulation.**

Touch a metal ball with a hot needle and warmth spreads. The physical law is
simple enough to say in one sentence:

> **Each point's temperature rises at a rate proportional to (the average of
> its neighbors − itself).**

Hotter than your neighbors → you cool. Cooler → you warm. A rule about
*rates of change* is called a **differential equation** — this one is the
famous **heat equation**. And on a mesh you can solve it with a loop. Try it
by hand first on a tetrahedron (4 corners, everyone neighbors everyone):
start with temperatures (100, 0, 0, 0), and repeatedly apply

```
new_temp = temp + 0.1 * (average_of_neighbors - temp)
```

Watch the 100 melt toward (25, 25, 25, 25). **You just solved a differential
equation** — that's all "solving" means: step forward in time, over and over.

Now the real thing. "Compare me with my neighbors" is — of course — a matrix
(the library calls it `L`; its weights are chosen carefully so the physics
comes out right on a curved surface). And "how big is each corner's patch"
is a matrix too (`M` — the diagonal matrix of Stage 6's vertex areas):

```python
from ddg import icosphere, cotangent_laplacian, mass_matrix
from ddg import HeatSolver, point_source

mesh = icosphere(4)
L, _ = cotangent_laplacian(mesh)   # the neighbor-comparison machine
M = mass_matrix(mesh)              # the patch-size machine (Stage 6!)
solver = HeatSolver(M, L, dt=0.01)
u = point_source(mesh.n_vertices, 0)     # all heat starts at corner 0
snapshots = solver.run(u, n_steps=30)
```

(The professional solver steps time using a clever, never-explodes method —
trust it for now, test it always.)

- Make a strip of pictures of the ball colored by temperature: watch the hot
  spot bloom and spread.
- **Test the physics.** Heat is never created or destroyed, only spread, so
  its integral must stay constant. Use *your own* `integrate(mesh, u)` from
  Stage 6 on every snapshot: same number to many decimal places, every time.
- Connect to Stage 7: `L @ ones` should be all zeros (if everyone has the
  same temperature, nobody changes). Check it.

**Done when:** you did the tetrahedron by hand, your picture strip exists,
and your own integral stays constant across all snapshots.

---

## Stage 10 — Distance from heat

**Math: the gradient, and a convergence study you run yourself.**

An airplane from New York to Madrid flies a curved arc — the shortest path
*along the surface*. Computing such distances is hard, and this repository
implements a beautiful 2013 discovery: **watch heat to find distance.** Heat
reaches nearby points sooner, so early warmth encodes how far everything is.
One more calculus word makes it work: the **gradient** of a temperature
field is, at each spot, the arrow pointing in the direction of fastest
increase (it's built from Stage 7's differences). Heat flows *along* those
arrows, away from the source — follow them backwards and you can recover
distance.

```python
from ddg import icosphere, HeatMethodSolver, exact_sphere_distance

mesh = icosphere(4)
dist = HeatMethodSolver(mesh).distance(0)   # distance from corner 0 to all
```

- Color the ball by `dist`: you should see rings, like ripples around the
  starting corner — bands of "equally far."
- On a unit ball the exact answer is known, and the library provides it
  (`exact_sphere_distance`). Subtract, average the absolute error.
- Repeat at detail levels 2–6 and table the errors, exactly as in Stage 6.
  You have run a **convergence study** — the core experiment of all
  computational science. Open `docs/error_report.md` and compare your table
  with the professionals'. They should agree.

**Done when:** your rings picture exists and your error table shows the
error shrinking as the mesh gets finer.

---

## Optional bonus — the ball's pure tones

A drumhead vibrates in special patterns, each with its own pitch. So does a
ball, and the matrices `L` and `M` know those patterns (they're called
**eigenvectors** — "a matrix's favorite vectors"). The library can find
them: `from ddg import eigenpairs`, then color the ball by a few columns of
the result and compare with `figures/eigenfunctions.png`. The pitches come
out as 0, 2, 2, 2, 6, 6, 6, 6, 6, ... — a pattern with a formula, ℓ(ℓ+1),
that physicists use to describe atoms. No "Done when" here; this is a
postcard from further up the mountain.

---

# Part III — Becoming a software builder

## Stage 11 — Tests: programs that check programs

**Skill: automated testing — the professional habit.**

You've been checking results by eye. Professionals write the checks *as
programs* so they run every time, forever. You own four rock-solid truths:

- V − E + F = 2 (Stage 2),
- a constant vector has zero differences: `D @ ones == 0` (Stage 7),
- angle defects total exactly 4π (Stage 8),
- total heat never changes (Stage 9).

Write each as a test function in `test_mytools.py`:

```python
from ddg import icosphere
from mytools import euler_number

def test_euler():
    assert euler_number(icosphere(3)) == 2
```

Install `pytest` and run `pytest test_mytools.py` — it finds every `test_`
function and reports pass or fail. Break something on purpose (flip a sign)
and watch a test catch it. *That feeling of being caught is the whole
point.* Then read the professionals' tests in `tests/test_ddg.py` — you'll
recognize most of what they check, including your four truths.

**Done when:** four passing tests, and you've seen one fail when you
sabotage the code.

---

## Stage 12 — The capstone: ship a real application

**Skill: putting everything together into something others can use.**

Build **Sphere Explorer**, a program someone else runs from the terminal:

```
python sphere_explorer.py --detail 4 --start 17
```

It should:

1. build the mesh at the requested detail level,
2. print a fact sheet — corners, edges, faces, Euler number, surface area
   (and its distance from 4π), total angle defect (always 4π!),
3. save three pictures: the mesh, heat spreading from the starting corner,
   and the distance rings,
4. refuse politely (a clear message, not a crash) if asked for detail
   level 50.

Requirements for "real software" — all things you now know:

- the logic lives in functions in `mytools.py`, not one giant script,
- your `pytest` tests pass,
- a `README.md` (this repository's is a good model) says what it is and how
  to run it,
- with your helper, put it on GitHub — and ask a friend to follow only your
  README on their computer. Whatever confuses them, fix.

**Done when:** someone who is not you, on a computer that is not yours, runs
Sphere Explorer successfully using only your README.

---

# The dictionary you've built

When you meet these words in school, you already know them:

| Official word | What you did |
|---|---|
| vector | a list of numbers — temperatures at every corner (Stage 5) |
| **integral** | weighted sum of values × patch areas; finer pieces → truer answer (Stage 6) |
| **derivative** | difference between neighbors (Stage 7) |
| **matrix** | a machine that eats a vector and produces a vector; you multiplied one by hand (Stage 7) |
| curvature | the missing angle at a corner (Stage 8) |
| **differential equation** | a rule about rates of change, solved by stepping in a loop (Stage 9) |
| gradient | the arrow pointing uphill fastest (Stage 10) |
| convergence | the error table that shrinks 4× per level (Stages 6 & 10) |
| eigenvector | a shape's pure tone (bonus stage) |

And the software column: variables, lists, loops, functions, libraries,
`numpy` vectors, plotting, automated tests, and a shipped application with a
README.

## Where to go next

- **The math:** Keenan Crane's free lectures and notes, *Discrete
  Differential Geometry: An Applied Introduction* (Carnegie Mellon). The
  first lectures will feel familiar — you've computed half the examples.
- **The code:** read `docs/theory.md`, then the source in `ddg/`. The one
  thing we treated as magic — exactly how the neighbor-comparison matrix `L`
  chooses its weights on a curved surface (look up "cotangent Laplacian") —
  is explained there, and understanding it is the natural next mountain.
