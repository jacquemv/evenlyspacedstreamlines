import numpy as np
from readmesh import read_off
from evenlyspacedstreamlines import evenly_spaced_streamlines

ver, tri = read_off('../data/test_finemesh2d.off')
orient = np.fromfile('../data/test_finemesh2d.orient', sep=' ').reshape((-1, 3))
radius = 0.02

i = 0
while True:
    i += 1
    print(f'iteration {i}')
    lines, faces, infos = evenly_spaced_streamlines(
        ver, tri, orient, radius)
