import numpy as np
from evenlyspacedstreamlines import evenly_spaced_streamlines


radius = 0.2

ver = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0]], dtype=np.float64)
tri = np.array([[0, 1, 2]], dtype=np.int32)
orient = np.array([[1, 1, 1]], dtype=np.float64)

lines, faces, infos = evenly_spaced_streamlines(
    ver, tri, orient, radius, orthogonal=True, allow_tweaking_orientation=True, parallel=False)
assert len(lines) > 0

lines, faces, infos = evenly_spaced_streamlines(
    ver[:, :2], tri, orient[:, :2], radius, orthogonal=True)
assert len(lines) > 0

lines, faces, infos = evenly_spaced_streamlines(
    ver, tri, orient*0, radius)
assert len(lines) == 0

lines, faces, infos = evenly_spaced_streamlines(
    ver, tri, np.full_like(orient, np.nan), radius)
assert len(lines) == 0


# import matplotlib.pyplot as plt
# plt.triplot(ver[:, 0], ver[:, 1], tri, c=(0.8,)*3)
# for line in lines:
#     plt.plot(line[:, 0], line[:, 1], '-')
# plt.show()
