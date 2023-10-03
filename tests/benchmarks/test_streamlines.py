import numpy as np
import matplotlib.pyplot as plt
from readmesh import read_off

radius = 0.2
theta = np.linspace(0, 2*np.pi, 32)
cx = radius * np.cos(theta)
cy = radius * np.sin(theta)

ver, tri = read_off('../data/test_finemesh2d.off')
plt.triplot(ver[:, 0], ver[:, 1], tri, c=(0.8,)*3)

i = 0
with open('out', 'rt') as file:
    for line in file.readlines():
        x = np.fromstring(line, sep=' ').reshape((-1, 3))
        i += 1
        plt.plot(x[:, 0], x[:, 1], '-')
        #plt.plot(x[0, 0]+cx, x[0, 1]+cy, c=(0.7,)*3)
        #plt.plot(x[-1, 0]+cx, x[-1, 1]+cy,  c=(0.7,)*3)
plt.axis('equal')
plt.show()
