# Learning Roadmap

**Who this is for:** a curious student (around 8th grade) with little or no
coding experience, who wants to learn two things at once:

1. how to write programs and build real software, and
2. the geometry of curved surfaces, using a computer (the grown-up name for
   this subject is *discrete differential geometry*).

**What this repository is:** a small, finished program written by
professionals. It builds a ball out of triangles, measures it, spreads heat
across it, and computes distances along its curved surface. You will not
understand most of it at first. That's fine — it's the *destination*. By the
end of this roadmap you will read it, use it, test it, and build your own
small application on top of it.

**The one rule:** every stage ends with something you made that *works* —
a number you computed, a picture you drew, or a program someone else can run.
Each stage has a **Done when** line. Don't move on until you've hit it.

You'll need a helper (parent, teacher, or older student) for about an hour at
the very start, to install Python and get this folder onto your computer.
After that, the roadmap is yours.

---

# Part I — Learning to code, using shapes

## Stage 1 — Hello, Python

**Coding skill: giving a computer instructions.**

- With your helper, install Python and open a terminal in this folder.
- Start Python by typing `python`. You're now talking to the computer
  directly. Try:

  ```python
  print("hello")
  2 + 2
  10 * 4 ** 3 + 2
  ```

- Learn the four things almost every program is made of: **variables**
  (giving a value a name), **lists** (a row of values), **loops** (do this
  for each item), and **if** (do this only when something is true). Any
  beginner Python tutorial covers these; spend a few days playing.

**Done when:** you can write a loop that prints the numbers 1 to 100, and
explain each line to your helper.

---

## Stage 2 — A shape is just lists

**Coding skill: representing real things as data.
Math idea: a solid shape can be written down as numbers.**

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

- By hand (paper!), count the cube's corners (V), edges (E), and faces (F).
  Compute **V − E + F**.
- Do the same for a triangular pyramid (tetrahedron) and an octahedron.
- You should get the same answer every time. This is **Euler's formula**:
  for any shape that could be inflated into a ball, V − E + F = 2. It was
  discovered in 1758 and it is the first theorem of this whole subject.
- Now write a small program: given the `faces` list, count the edges with a
  loop (careful — each edge is shared by two faces, so don't count it twice).
  Check V − E + F = 2 *by program*.

**Done when:** your program prints `V - E + F = 2` for the cube, and you
checked the tetrahedron by hand.

---

## Stage 3 — Using someone else's code

**Coding skill: libraries — building on other people's work.**

Real programmers rarely start from nothing. This repository is a **library**:
a box of tools you can import. Time to open the box.

```python
from ddg import icosphere

mesh = icosphere(2)          # a ball made of triangles, detail level 2
print(mesh.n_vertices)       # how many corners?
print(mesh.n_faces)          # how many triangles?
print(mesh.euler_characteristic())   # V - E + F ... it knows the formula!
```

- Try detail levels 0 through 6 in a loop. Watch the counts grow. Confirm
  Euler's formula gives 2 *every single time*, even with 81,920 triangles.
- Look at `mesh.vertices` — it's the same idea as your cube's `corners`
  list, just longer, and every corner sits at distance exactly 1 from the
  center (it's a unit ball). Check that with the Pythagorean theorem:
  for a corner (x, y, z), is x² + y² + z² equal to 1?

**Done when:** you've verified, with your own loop, that V − E + F = 2 at
every detail level, and that the corners really lie on the ball.

---

## Stage 4 — Pictures

**Coding skill: making the computer draw.**

Numbers are good; pictures are better. The library `matplotlib` draws.

- Plot the mesh's corners as dots in 3D (your helper or an online tutorial
  can show you `matplotlib`'s 3D scatter plot — it's about five lines).
- Color each dot by its height (the z value). Now you can *see* the ball.
- Open the `figures/` folder in this repository and look at the pictures the
  professional program made. By the end of this roadmap you'll know what
  every one of them means — and you'll have made your own versions.

**Done when:** you have a picture file you made yourself showing the
triangle ball, and you've emailed it to someone or pinned it on the wall.

---

## Stage 5 — Functions: your own tools

**Coding skill: functions — naming a recipe so you can reuse it.**

You've been writing loose lines of code. A **function** wraps a recipe so
you can use it again:

```python
def euler_number(mesh):
    return mesh.n_vertices - mesh.n_edges + mesh.n_faces
```

- Write three functions of your own: `euler_number(mesh)`,
  `longest_edge(mesh)`, and `total_area(mesh)` (for area, the library helps:
  `from ddg import face_areas`, then add up `face_areas(mesh)` with `sum()`).
- Put them in your own file, `mytools.py`, and import them. You now have a
  library *of your own*. That's what software is: tools made of tools.

**Done when:** `from mytools import total_area` works, and
`total_area(icosphere(5))` prints a number close to **12.566...**
(keep that number in mind — Stage 6 explains it).

---

# Part II — The geometry of curved surfaces

## Stage 6 — Flat triangles can't quite make a ball

**Math idea: approximation and convergence — the heart of the subject.**

The true surface area of a ball of radius 1 is 4π ≈ **12.566** (a famous
formula). Your triangle ball is made of *flat* triangles, which cut corners,
so its area is a little *less*.

- Use your `total_area` on detail levels 2, 3, 4, 5, 6. Subtract each from
  4π. Make a table.
- Notice: every time the triangles get smaller, the error shrinks — by about
  **4×** per level. The flat model *converges* to the curved truth.

This is the single most important idea in computational geometry: **we can't
store a perfectly curved surface in a computer, but we can get as close as we
want with enough flat triangles, and we can measure exactly how close.**

**Done when:** you have the table, and you can say in one sentence why the
triangle area is always a little under 12.566.

---

## Stage 7 — Curvature: the missing angle

**Math idea: what "curved" means, with scissors and paper.**

Do this first *without* a computer:

- On flat paper, the angles around any point add to exactly **360°**.
- Now look at one corner of your cube: three squares meet there, each
  contributing 90°. Total: **270°**. There are 90° *missing*. That missing
  angle is why the corner pokes out — it's called the **angle defect**, and
  it is what mathematicians mean by curvature being concentrated there.
- The cube has 8 corners × 90° missing = **720°** in total.
- Try the tetrahedron: 4 corners, each with three 60° angles → 180° missing
  each → total **720°** again. Coincidence?

It is not. **For every shape that can be inflated into a ball, the missing
angles always total exactly 720°** (two full turns). This is the famous
**Gauss–Bonnet theorem**, and your triangle ball obeys it too:

- Write a program: the library gives every triangle's three angles
  (`from ddg.geometry import interior_angles`). For each corner of the mesh,
  add up all the angles that touch it, subtract from 360°, and total the
  results over all corners.
- You should get 720° (the code works in radians: 4π ≈ 12.566 — that number
  again, a fun coincidence of the sphere). Try every detail level: it's
  **exact every time**, not approximate like the area was. Some things in
  this subject converge; some are perfectly true on every mesh. Knowing
  which is which is real mathematical taste.

**Done when:** paper says 720° for the cube, your program says 720° for
every icosphere, and you can explain the difference between this (exact) and
Stage 6 (approximate).

---

## Stage 8 — Heat spreading on a surface

**Math idea: simulating a physical process.**

Touch a metal ball with a hot needle: warmth spreads outward. This
repository simulates exactly that, and you can drive it:

```python
from ddg import icosphere, cotangent_laplacian, mass_matrix
from ddg import HeatSolver, point_source

mesh = icosphere(4)
L, _ = cotangent_laplacian(mesh)    # the "how heat flows" matrix
M = mass_matrix(mesh)               # the "how big each region is" matrix
solver = HeatSolver(M, L, dt=0.01)
u = point_source(mesh.n_vertices, 0)   # all the heat starts at corner 0
snapshots = solver.run(u, n_steps=30)
```

You don't need to know how `L` and `M` work inside yet — using a tool before
you understand its internals is normal and honest, *as long as you test it*.

- Make a picture of the ball colored by temperature at several snapshots
  (Stage 4 skills). Watch the hot spot bloom and spread.
- Test it: physics says heat is never created or destroyed, only spread.
  Use the library's `total_heat(M, u)` on every snapshot — the total should
  stay the same to many decimal places.

**Done when:** you have a strip of pictures showing heat spreading, and a
printout showing total heat staying constant.

---

## Stage 9 — Distance along a curved surface

**Math idea: the shortest path on a sphere isn't a straight line.**

An airplane from New York to Madrid flies a curved arc — the shortest path
*along the surface*. Computing such distances is genuinely hard, and this
repository implements a beautiful 2013 discovery: **you can find distances by
watching heat.** Heat reaches nearby points sooner — the pattern of warmth
encodes distance (think of how you locate a campfire with your eyes closed).

```python
from ddg import icosphere, HeatMethodSolver

mesh = icosphere(4)
solver = HeatMethodSolver(mesh)
dist = solver.distance(0)    # distance from corner 0 to every other corner
```

- Picture time: color the ball by `dist`. You should see rings, like ripples
  around the starting point — those are "equally far" bands.
- Check it against truth. On a unit ball the exact answer is known (the
  library provides it: `exact_sphere_distance`). Subtract, and find the
  average error.
- Repeat at detail levels 2–6 and table the errors, exactly like Stage 6.
  Smaller triangles → smaller error. You have just done a **convergence
  study**, which is precisely what the professional benchmark in this
  repository does — open `docs/error_report.md` and compare your table with
  theirs.

**Done when:** your rings picture exists, and your error table shows the
error shrinking as the mesh gets finer.

---

# Part III — Becoming a software builder

## Stage 10 — Tests: programs that check programs

**Coding skill: automated testing — the professional habit.**

You've been checking results by eye. Professionals write the checks *as
programs* so they run every time, forever. You already know three rock-solid
truths:

- V − E + F = 2 (Stage 2),
- angle defects total 4π exactly (Stage 7),
- total heat never changes (Stage 8).

Write each as a test function in a file `test_mytools.py`:

```python
from ddg import icosphere
from mytools import euler_number

def test_euler():
    assert euler_number(icosphere(3)) == 2
```

Install `pytest` and run `pytest test_mytools.py` — it finds and runs every
`test_` function and reports pass or fail. Now break something on purpose (a
minus sign somewhere) and watch the test catch it. *That feeling of being
caught is the whole point.* Read the professionals' tests in
`tests/test_ddg.py` — you will now recognize most of what they check.

**Done when:** you have at least three passing tests, and you've seen one
fail when you sabotage the code.

---

## Stage 11 — The capstone: ship a real application

**Coding skill: putting everything together into something others can use.**

Build **Sphere Explorer**: a program someone else can run from the terminal:

```
python sphere_explorer.py --detail 4 --start 17
```

It should:

1. build the mesh at the requested detail level,
2. print a fact sheet — corners, edges, faces, Euler number, surface area
   (and how close to 4π it is), total angle defect (always 4π!),
3. save three pictures into a folder: the mesh, heat spreading from the
   starting corner, and the distance rings,
4. refuse politely (a clear message, not a crash) if the user asks for
   detail level 50.

Requirements for "real software," all things you now know how to do:

- the logic lives in functions in your `mytools.py`, not one giant script,
- your `pytest` tests pass,
- a `README.md` (look at this repository's for inspiration) explains in a
  few sentences what it is and how to run it,
- with your helper, put it on GitHub — and ask a friend to follow only your
  README on their computer. Whatever confuses them, fix.

**Done when:** someone who is not you, on a computer that is not yours, runs
Sphere Explorer successfully using only your README.

---

# What you will have learned

| | Coding | Geometry |
|---|---|---|
| 1 | talking to Python: variables, lists, loops | — |
| 2 | data: representing things as lists | Euler's formula V−E+F=2 |
| 3 | using libraries | meshes; points on a sphere |
| 4 | drawing with matplotlib | seeing 3D data |
| 5 | writing functions; your own module | edge lengths, areas |
| 6 | experiments with loops and tables | approximation & convergence; 4π |
| 7 | translating paper math into code | curvature; Gauss–Bonnet (720°) |
| 8 | using solvers; testing physics | the heat equation |
| 9 | a real numerical experiment | geodesic distance; the heat method |
| 10 | automated tests with pytest | exact vs approximate truths |
| 11 | building & shipping an application | all of the above, in one tool |

## Where to go next

- **The math:** Keenan Crane's free video lectures and notes, *Discrete
  Differential Geometry: An Applied Introduction* (Carnegie Mellon). The
  first lectures will feel familiar — you've already computed half the
  examples.
- **The code:** read this repository's `docs/theory.md` and then the source
  files in `ddg/`. The mysterious matrices `L` and `M` from Stage 8 are
  explained there; understanding *how they're built* is the natural next
  mountain, and it's the door into the rest of this subject.
