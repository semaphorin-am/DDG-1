from ddg import icosphere

for subdivision in range(0,7):
    mesh = icosphere(subdivision)

    print (mesh.n_vertices)
    print(mesh.n_edges)
    print(mesh.n_faces)

    edges_set = set()

    for face in mesh.faces:
        a = face[0]
        b = face[1]
        c = face[2]
        edges_set.add((min(a,b),max(a,b)))
        edges_set.add((min(b,c),max(b,c)))
        edges_set.add((min(c,a),max(c,a)))

    if mesh.n_edges == len(edges_set):
        print(f"Subdivision level: {subdivision}")
        print(f"Mesh count: {mesh.n_edges}")
        print(f"My count: {len(edges_set)}")
        print("Computer and my count match")
