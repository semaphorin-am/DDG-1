from ddg import icosphere

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

mesh = icosphere(subdivisions=3)

print(mesh.vertices)

x_values = []
y_values = []
z_values = []

for vertex in mesh.vertices:
    x = vertex[0]
    y = vertex[1]
    z = vertex[2]
    x_values.append(x)
    y_values.append(y)
    z_values.append(z)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111,projection='3d')
scatter = ax.scatter(x_values, y_values, z_values, c=z_values, cmap='viridis')
ax.set_zlabel('Z')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.colorbar(scatter, ax=ax, label='Yellow represents a higher z valueas purple represents lower')
ax.set_title('Icosphere Mesh-Colored by Height')
plt.savefig('icosphere.png')
plt.close()