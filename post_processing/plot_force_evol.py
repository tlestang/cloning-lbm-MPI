import os
import numpy as np
import matplotlib.pyplot as plt

L = 32
nu = (1.0/3)*(0.502-0.5)

for i in range(0,128):
    print i
    filename = '/home/tlestang/tailleur_lecomte_MPI/output_11_03/clone_'+str(i)+'/evolution_1_clone_'+str(i)
    f1 = open(filename, 'rb')
    u = np.fromfile(f1, dtype=np.double)
    N = np.size(u);
    t0 = 17950.0;
    dt = 1.0/t0;
    tMax = N*dt;
    tDomain = np.linspace(0,tMax, N);
    plot(tDomain, u)

plt.show()
# plt.savefig('corr_32_zoom.png')
