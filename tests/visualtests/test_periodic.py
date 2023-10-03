import numpy as np
import matplotlib.pyplot as plt
from evenlyspacedstreamlines import evenly_spaced_streamlines, streamlines_to_tubes

radius = 0.04
n = 20
r1, r2 = 1, 1.3

theta = np.linspace(0, 2*np.pi, n+1)[:-1]
x1, y1 = r1*np.cos(theta), r1*np.sin(theta)
theta += theta[1]/2
x2, y2 = r2*np.cos(theta), r2*np.sin(theta)

ver = np.vstack((
    np.column_stack((x1, y1)),
    np.column_stack((x2, y2))
))
I = np.arange(n)
tri = np.vstack(( 
    np.column_stack((I, (I+1)%n, I+n)),
    np.column_stack(((I+1)%n, (I+1)%n + n, I+n))
)).astype(np.int32)

plt.triplot(ver[:, 0], ver[:, 1], tri, c=(0.8,)*3)

G = (ver[tri[:, 0]] + ver[tri[:, 1]] + ver[tri[:, 2]]) / 3
orient = np.column_stack((G[:, 1]+1.1e-4, -G[:, 0]+1e-4))

I = 2 *G[:, 0] + G[:, 1] < 0
orient[I] = -orient[I]


plt.quiver(G[:, 0], G[:, 1], orient[:, 0], orient[:, 1])

lines, faces, infos = evenly_spaced_streamlines(
    ver, tri, orient, radius, random_seed=1656514220, oriented_streamlines=True)

theta = np.linspace(0, 2*np.pi, 100)
x, y = radius*np.cos(theta), radius*np.sin(theta)
for line in lines:
    plt.plot(line[:, 0], line[:, 1], c='C0')
    plt.plot(line[0, 0]+x, line[0, 1]+y, c='C1')
    plt.plot(line[-1, 0]+x, line[-1, 1]+y, c='C3')


plt.axis('equal')
plt.show()