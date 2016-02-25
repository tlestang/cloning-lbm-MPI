import os
import numpy as np
import matplotlib.pyplot as plt

L = 32
nu = (1.0/3)*(0.502-0.5)

filename = '../output/clone_2/evolution_0_clone_2'
f1 = open(filename, 'rb')
u = np.fromfile(f1, dtype=np.double)

# filename = '../output/clone_1/evolution_0_clone_1'
# f1 = open(filename, 'rb')
# u1 = np.fromfile(f1, dtype=np.double)

# filename = '../output/clone_2/evolution_0_clone_2'
# f1 = open(filename, 'rb')
# u2 = np.fromfile(f1, dtype=np.double)

# filename = '../output/clone_3/evolution_0_clone_3'
# f1 = open(filename, 'rb')
# u3 = np.fromfile(f1, dtype=np.double)

# filename = 'CLONE_2'
# f1 = open(filename, 'rb')
# u2 = np.fromfile(f1, dtype=np.double)
print u

N = np.size(u);
t0 = 5857.0;
dt = 1.0/t0;
tMax = N*dt;
tDomain = np.linspace(0,tMax, N);

plt.plot(tDomain, u)
# plt.plot(tDomain, u1)
# plt.plot(tDomain, u2)
# plt.plot(tDomain, u3)
plt.show()
# plt.savefig('corr_32_zoom.png')
