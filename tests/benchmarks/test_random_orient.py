import numpy as np
from readmesh import read_off
from evenlyspacedstreamlines import evenly_spaced_streamlines

ver, tri = read_off('../data/test_finemesh2d.off')
radius = 0.05

i = 0
while True:
    i += 1
    print(f'iteration {i}')
    lines, faces, infos = evenly_spaced_streamlines(
        ver, tri, np.random.normal(size=tri.shape), radius)
    print(f'nl = {len(lines)}')
