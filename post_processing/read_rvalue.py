import os
import numpy as np
import matplotlib.pyplot as plt

filename = '/home/thibault/tailleur_lecomte/lbm/test_eps_0.003_T_200/rvalues.dat'
f1 = open(filename, 'rb')
r = np.fromfile(f1, dtype=np.double)


T = 0
N = np.size(r)
lam = np.zeros(N)

for i in range(1,N):
    s = 0
    for j in range(1,i):
        s = s + np.log(r[j])
    T = T + 0.3
    oneOvT = 1./T
    lam[i] = oneOvT*s
    print lam[i]

tDomain = np.linspace(0,T,N)

#plt.plot(tDomain, lam)          
plt.plot(tDomain, r)
plt.show()
