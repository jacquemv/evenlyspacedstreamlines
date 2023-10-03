import numpy as np
import matplotlib.pyplot as plt
from evenlyspacedstreamlines import evenly_spaced_streamlines
from readmesh import read_off

geom = 'moebius'
radius = 0.05

geom = 'klein'
radius = 0.5

ver, tri = read_off('../data/test_'+geom+'.off')
orient = np.fromfile('../data/test_'+geom+'.orient', sep=' ').reshape((-1, 3))

lines, idtri, infos = evenly_spaced_streamlines(ver, tri, orient, radius, seed_points=1024)

# with open('klein_lines.txt', 'tw') as file:
#     for line in lines:
#         line.ravel().tofile(file, sep=' ')
#         file.write('\n')

plt.figure().add_subplot(projection='3d')
#plt.gca().plot_trisurf(ver[:, 0], ver[:, 1], ver[:, 2], triangles=tri, alpha=0.2, color=(0.7,)*3)
for line in lines:
    plt.plot(line[:, 0], line[:, 1], line[:, 2])

plt.gca().view_init(67, 124)
plt.axis('equal')
plt.axis('off')
plt.show()