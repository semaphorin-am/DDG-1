from ddg import icosphere 

for subdivisions in range(0,7):

    mesh = icosphere(subdivisions)
    print(mesh.n_vertices)
    print(mesh.euler_characteristic())

    edges_set = set()

    for face in mesh.faces:
        a = face[0]
        b = face[1]
        c = face[2]
        edges_set.add((min(a,b), max(a,b)))
        edges_set.add((min(b,c), max(b,c)))
        edges_set.add((min(c,a), max(c,a)))
    print(len(edges_set))

    n_vertices = (mesh.n_vertices)
    n_edges = len(edges_set)
    n_faces = len(mesh.faces)
   
    print(f"V = {n_vertices}")
    print(f"E = {n_edges}")
    print(f"F = {n_faces}")
    print(f"V - E + F = {n_vertices - n_edges + n_faces}")

count_on_sphere = 0 

for vertex in mesh.vertices:
    (x,y,z) = vertex
    sum_squares = x**2+y**2+z**2
    if sum_squares > 0.999 and sum_squares < 1.0001:
        count_on_sphere = count_on_sphere + 1

print(f"Vertices on sphere: {count_on_sphere} out of {len(mesh.vertices)}")
