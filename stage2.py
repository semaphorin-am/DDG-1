faces = [(0,1,2,3), (4,5,6,7),
(0,1,5,4), (1,2,6,5), (2,3,7,6), (3,0,4,7),]

edges_set = set()


for face in faces:
    a = face[0]
    b = face[1]
    c = face[2]
    d = face[3]
    edges_set.add((min(a,b), max(a,b)))
    edges_set.add((min(b,c), max(b,c)))
    edges_set.add((min(c,d), max(c,d)))
    edges_set.add((min(d,a), max(d,a)))
print(len(edges_set))

n_vertices = 8
n_faces = 6
n_edges = len(edges_set)

print(f"V = {n_vertices}")
print(f"E = {n_edges}")
print(f"F = {n_faces}")
print(f"V - E + F = {n_vertices - n_edges + n_faces}")