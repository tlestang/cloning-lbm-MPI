import os
import numpy as np
import matplotlib.pyplot as plt

filename = ''
f1 = open(filename, 'rb')
r = np.fromfile(f1, dtype=np.double)
q

T = 10

N = np.size(u)

for i in range(1,N):
    s = 0
    for j in range(1,i):
        s = s + np.log(r[j])
    lam[i] = oneOvT*s
    

tDomain = np.linspace(1,N,N)

plt.plot(tDomain, lam)
plt.show()
