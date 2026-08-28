from ddg import icosphere 

def euler_number(mesh):
    V = (mesh.n_vertices)
    E = (mesh.n_edges)
    F = (mesh.n_faces)
   
    result = V - E + F
    return result

result = euler_number(icosphere(3))
print(f"Euler Number: {result}")

def longest_edge(mesh):
    import math

    max_distance = 0 

    for edge in mesh.edges:
        v1_idx, v2_idx = edge
        v1 = mesh.vertices[v1_idx]  # ← Get vertex 1's coordinates
        v2 = mesh.vertices[v2_idx]  # ← Get vertex 2's coordinates
        x1, y1, z1 = v1            # ← Unpack vertex 1
        x2, y2, z2 = v2            # ← Unpack vertex 2
        distance = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

        if distance > max_distance: 
            max_distance = distance 
    return max_distance

result_edge = longest_edge(icosphere(3))
print(f"longest edge = {result_edge}")


def total_area(mesh):
    from ddg import face_areas
    total_area = 0   
    
    for area in face_areas(mesh):
        total_area = total_area + area

    return total_area

result_total_area = total_area(icosphere(3))
print(f"Total area of mesh = {result_total_area}")
