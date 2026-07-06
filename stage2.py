"""Stage 2: Count edges from faces, check Euler's formula.

Challenge: given only the faces list, count the edges by yourself,
accounting for the fact that each edge is shared by two faces.
"""

from ddg import icosphere

# Create a small mesh to test with
mesh = icosphere(subdivisions=0)

# Extract the data
vertices = mesh.vertices
faces = mesh.faces
n_vertices = len(vertices)

print(f"Mesh has {n_vertices} vertices and {len(faces)} faces.")
print(f"Vertices: {vertices}")
print(f"Faces:\n{faces}")

# PUZZLE: Count the edges without using mesh.edges or mesh.n_edges.
# You have only the faces list.
#
# Hint: think about storing edges in a set to avoid counting them twice.
# An edge between vertex i and j can be represented as a tuple (i, j)
# where you always put the smaller index first, so (3, 7) and (7, 3)
# are the same edge.
#
# Write your edge-counting logic here:

# Your code goes below this line
# ...


# Once you count the edges, compute the Euler characteristic
# and verify it equals 2.
print(f"\nV = {n_vertices}")
print(f"E = {n_edges}")  # Replace with your answer
print(f"F = {len(faces)}")
print(f"V - E + F = {n_vertices - n_edges + len(faces)}")  # Should be 2
