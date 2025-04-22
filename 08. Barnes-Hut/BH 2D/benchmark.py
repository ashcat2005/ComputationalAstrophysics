import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('09. Barnes-Hut/BH 2D/Benchmark.txt',delimiter=',')

plt.figure()
plt.plot(data[:,0], data[:,1])
plt.xlabel(f'Número de cuerpos')
plt.ylabel(f'Tiempo de computo [s]')
plt.title(f'Tiempo de computo de 100 pasos de integración')
plt.show()