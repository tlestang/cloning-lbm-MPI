import os
import numpy as np
import matplotlib.pyplot as plt

filename = '../output_test/copies_evolution_1'
f1 = open(filename, 'rb')
u = np.fromfile(f1, dtype=np.int32)

print u[8:16]
