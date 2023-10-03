import numpy as np
import matplotlib.pyplot as plt
from evenlyspacedstreamlines import evenly_spaced_streamlines
from readmesh import read_off

radius = 0.2
ver, tri = read_off('../data/test_mesh2d.off')
orient = np.fromfile('../data/test_mesh2d.orient2', sep=' ').reshape((-1, 3))

lines, faces, infos = evenly_spaced_streamlines(
    ver, tri, orient, radius, random_seed=1656514220)

def colinear(x, y, z):
    return np.linalg.norm(np.cross(x-z, y-z)) < 1e-6

def point_in_tri(p_test, p0, p1, p2):
     dX = p_test[0] - p0[0]
     dY = p_test[1] - p0[1]
     dX20 = p2[0] - p0[0]
     dY20 = p2[1] - p0[1]
     dX10 = p1[0] - p0[0]
     dY10 = p1[1] - p0[1]
     s_p = (dY20*dX) - (dX20*dY)
     t_p = (dX10*dY) - (dY10*dX)
     D = (dX10*dY20) - (dY10*dX20)
     if D > 0:
         return ((s_p >= 0) and (t_p >= 0) and (s_p + t_p) <= D)
     else:
         return ((s_p <= 0) and (t_p <= 0) and (s_p + t_p) >= D)
        
# check if each segment is in the expected triangle
for line, face in zip(lines, faces):
    for i in range(face.size):
        x1, x2 = line[i], line[i+1]
        j = face[i]
        a = ver[tri[j, 0]]
        b = ver[tri[j, 1]]
        c = ver[tri[j, 2]]
        b1 = colinear(x1, a, b) or colinear(x1, a, c) or colinear(x1, b, c)
        b2 = colinear(x2, a, b) or colinear(x2, a, c) or colinear(x2, b, c)
        assert b1 or b2
        if not b1:
            assert point_in_tri(x1, a, b, c)
        if not b2:
            assert point_in_tri(x2, a, b, c)

