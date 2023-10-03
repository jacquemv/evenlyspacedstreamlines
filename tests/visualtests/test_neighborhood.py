import os
import numpy as np
import matplotlib.pyplot as plt
from readmesh import read_off


ver, tri = read_off('../data/test_mesh2d.off')

with open('test_mesh2d.neightri', 'rt') as file:
    radius = float(file.readline())
    lines = file.readlines()
os.remove('test_mesh2d.neightri')
os.remove('test')

theta = np.linspace(0, 2*np.pi, 200)
cx = radius * np.cos(theta)
cy = radius * np.sin(theta)
i = 0
for line in lines:
    I = np.array([int(s) for s in line.split()])
    plt.triplot(ver[:, 0], ver[:, 1], tri, c=(0.8,)*3)
    plt.triplot(ver[:, 0], ver[:, 1], tri[I], c='r')
    for k in range(3):
        G = ver[tri[i, k]]
        plt.plot(G[0] + cx, G[1] + cy, 'C0')
    plt.axis('equal')
    plt.show()
    i += 1

