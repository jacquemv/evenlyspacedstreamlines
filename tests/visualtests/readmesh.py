import numpy as np

def read_off(file):
    with open(file, 'rt') as file:
        file.readline()
        nv, nt, _ = [int(x) for x in file.readline().split()]
        ver = np.fromfile(file, sep=' ', dtype=np.float32, count=nv*3)
        tri = np.fromfile(file, sep=' ', dtype=np.int32, count=nt*4)
    return ver.reshape((-1, 3)), tri.reshape((-1, 4))[:, 1:]