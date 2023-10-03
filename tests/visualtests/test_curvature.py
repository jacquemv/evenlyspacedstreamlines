import os
import numpy as np
import matplotlib.pyplot as plt

with open('out', 'rt') as file:
    n, r = file.readline().split()
    n, r = int(n), float(r)
    x = np.fromfile(file, count=3*n, sep=' ').reshape((n, 3))
    start, length = file.readline().split()
    start, length = int(start), int(length)


plt.plot(x[:, 0], x[:, 1])
plt.plot(x[start:start+length+1, 0], x[start:start+length+1, 1])
for i in range(x.shape[0]-1):
    p = x[i, :2]
    q = x[i+1, :2]
    T = q - p
    N = np.array([T[1], -T[0]])
    N /= np.linalg.norm(N)
    y = np.row_stack((p-r*N, p+r*N, q+r*N, q-r*N, p-r*N))
    plt.plot(y[:, 0], y[:, 1], c=(0.8,)*3, zorder=-20)
plt.axis('equal')
plt.show()