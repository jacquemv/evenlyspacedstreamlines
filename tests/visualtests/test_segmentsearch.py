import subprocess
import numpy as np
import matplotlib.pyplot as plt

b = 0.96
step = 0.1
radius = 0.06
N = 10000

proc = subprocess.Popen(['./test', str(b), str(step), str(radius), str(N)], stdout=subprocess.PIPE)
out, err = proc.communicate()
out = out.decode().replace('\n', ' ')
x = np.array([float(x) for x in out.split(' ')[:-1]]).reshape((-1, 2))

cut = np.where(x[:, 0] > 1e9)[0][0]
x[cut] = np.nan
x = np.vstack((x[cut-1::-1], x[cut:]))

plt.plot(x[:, 0], x[:, 1], '-')
t = np.linspace(0, 2*np.pi, 50)
plt.plot(x[-1, 0] + radius * np.cos(t), x[-1, 1] + radius * np.sin(t))
plt.plot(x[0, 0] + radius * np.cos(t), x[0, 1] + radius * np.sin(t))
plt.axis('equal')
plt.show()