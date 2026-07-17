from ddg import icosphere 

def euler_number(mesh):
    print(f"{mesh.n_vertices}")
    print(f"{mesh.n_edges}")
    print(f"{mesh.n_faces}")

    V = (mesh.n_vertices)
    E = (mesh.n_edges)
    F = (mesh.n_faces)
    result = V - E + F
    return result

result = euler_number(icosphere(3))
print(result)


def longest_edge(mesh):
    import math

    max_distance = 0 

    for edge in mesh.edges:
        print(edge)
        distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2 - z1)**2)**0.5

        if distance > max_distance: 
            max_distance = distance 
    return max_distance

result_edge = max_distance
print(result_edge)


