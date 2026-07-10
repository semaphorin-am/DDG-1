"""Stage 3: Using Someone Else's Code - icosphere()

Challenge: use the library's icosphere() function to generate meshes at different
detail levels, then verify Euler's formula and the sphere equation on each.
"""

from ddg import icosphere

# PART 1: Generate an icosphere and understand its structure
# =========================================================

# Create a simple icosphere (subdivision level 0)
mesh = icosphere(subdivisions=0)

# Let's look at what we get
print("=== ICOSPHERE (subdivisions=0) ===")
print(f"Number of vertices: {len(mesh.vertices)}")
print(f"Number of faces: {len(mesh.faces)}")
print(f"Number of edges (from mesh object): {len(mesh.edges)}")

# PART 2: Count edges yourself (like Stage 2, but with a different mesh)
# =====================================================================
# You already know how to do this! Use a set to avoid double-counting.
# Write your edge-counting code here:

# Create an empty set for edges
edges_set = set()

# Loop through each face
for face in mesh.faces:
    # Each face is a triangle, so it has 3 vertices
    v1, v2, v3 = face

    # Add the three edges (remember: normalize with min/max)
    edges_set.add((min(v1, v2), max(v1, v2)))
    edges_set.add((min(v2, v3), max(v2, v3)))
    edges_set.add((min(v3, v1), max(v3, v1)))

n_edges_counted = len(edges_set)

# PART 3: Verify Euler's formula
# ==============================
v = len(mesh.vertices)
e = n_edges_counted
f = len(mesh.faces)

euler = v - e + f

print(f"\nV = {v}")
print(f"E = {e}")
print(f"F = {f}")
print(f"V - E + F = {euler}")
print(f"✓ Euler's formula holds!" if euler == 2 else "✗ Something's wrong!")


# PART 4: Verify the sphere equation x² + y² + z² = 1
# ===================================================
# Take a few vertices and check if they're on the unit sphere

print("\n=== CHECKING SPHERE EQUATION ===")
print("Taking the first 5 vertices and checking x² + y² + z² = 1:")

for i in range(min(5, len(mesh.vertices))):
    x, y, z = mesh.vertices[i]
    sum_squares = x**2 + y**2 + z**2
    print(f"Vertex {i}: ({x:.4f}, {y:.4f}, {z:.4f}) → x²+y²+z² = {sum_squares:.6f}")


# CHALLENGE: Now try different subdivision levels
# ================================================
print("\n=== TESTING MULTIPLE SUBDIVISION LEVELS ===")
print("Does Euler's formula hold at all detail levels?\n")

for subdiv in range(7):
    mesh = icosphere(subdivisions=subdiv)

    # Count edges
    edges_set = set()
    for face in mesh.faces:
        v1, v2, v3 = face
        edges_set.add((min(v1, v2), max(v1, v2)))
        edges_set.add((min(v2, v3), max(v2, v3)))
        edges_set.add((min(v3, v1), max(v3, v1)))

    v = len(mesh.vertices)
    e = len(edges_set)
    f = len(mesh.faces)
    euler = v - e + f

    print(f"Subdiv {subdiv}: V={v:6d}  E={e:6d}  F={f:6d}  →  V-E+F = {euler}")


# Thought questions (answer in your head or on paper):
# =====================================================
# 1. As the mesh gets more detailed (higher subdivision), do V, E, and F all increase?
# 2. Does Euler's formula (V - E + F = 2) *always* hold, no matter the detail level?
# 3. Why is the icosphere "better" than a cube for representing a sphere?
