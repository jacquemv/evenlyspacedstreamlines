import sys
import numpy as np
import matplotlib.pyplot as plt
from readmesh import read_off

ver, tri = read_off(sys.argv[1])
radius = float(sys.argv[2])

theta = np.linspace(0, 2*np.pi, 90)
cx = radius * np.cos(theta)
cy = radius * np.sin(theta)

def read_streamlines(fname):
    curves = []
    with open(fname, "rt") as file:
        for line in file.readlines():
            x = np.fromstring(line, sep=' ').reshape((-1, 3))
            curves.append(x)
    return curves

def plot_curves(C, *args, **kwargs):
    for x in C:
        plt.plot(x[:, 0], x[:, 1], *args, **kwargs)

def plot_mask(i):
    with open('output/mask'+str(i)+'.txt', 'r') as file:
        for line in file.readlines():
            x = np.fromstring(line, sep=' ').reshape((-1, 3))
            plt.fill(x[:, 0], x[:, 1], c=(247/256, 215/256, 215/256))

C = []
for k in range(1, 6):
    C.append(read_streamlines('output/streamlines'+str(k)+'.txt'))

orient = np.fromfile('output/orient.txt', sep=' ').reshape((-1, 6))
Vx = np.column_stack((orient[:, 0], 
                      orient[:, 3], 
                      np.full(orient.shape[0], np.nan)))
Vy = np.column_stack((orient[:, 1], 
                      orient[:, 4], 
                      np.full(orient.shape[0], np.nan)))

blocks = np.fromfile('output/blocks.txt', sep=' ').astype(int).reshape((-1, 2))
unique_blocks, unique_counts = np.unique(np.sort(blocks, -1), 
                                    axis=0, return_counts=True)
assert np.all((unique_counts == 1) | (unique_counts == 2))
Bx, By = [0, 0], [0, 0]
for val in 1, 2:
    I = unique_counts == val
    Bx[val-1] = np.column_stack((ver[unique_blocks[I, 0], 0], 
                        ver[unique_blocks[I, 1], 0], 
                        np.full(I.sum(), np.nan)))
    By[val-1]  = np.column_stack((ver[unique_blocks[I, 0], 1], 
                        ver[unique_blocks[I, 1], 1], 
                        np.full(I.sum(), np.nan)))


for i in range(0, len(C[0])):
    plt.figure(figsize=(14, 10))
    plot_mask(i)
    plt.triplot(ver[:, 0], ver[:, 1], tri, c=(0.8,)*3)
    plt.plot(Vx.ravel(), Vy.ravel(), c=(0.9,)*3)
    plt.plot(Bx[0].ravel(), By[0].ravel(), 'k--')
    plt.plot(Bx[1].ravel(), By[1].ravel(), 'k-')
    plot_curves(C[-1][:i], '.-', c=(0.6,)*3)
    for step in range(5):
        x = C[step][i]
        plt.plot(x[:, 0], x[:, 1], '.-', c='C'+str(step), 
                 label='step '+str(step+1))
    for x0 in C[-1][i]:
        plt.plot(x0[0]+cx, x0[1]+cy, 'k:')
    plt.legend()
    plt.title('Streamline '+str(i))
    plt.show()