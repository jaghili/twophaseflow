import matplotlib.pyplot as plt
import numpy as np

datiF = np.loadtxt("physics_constant_tau.dat", skiprows=1);
datC = np.loadtxt("physics_constant_sat.dat", skiprows=1);

plt.subplot(121)
plt.title("Capillary pressures at Interfaces")
plt.grid(linewidth=1)
plt.plot(datiF[:,3], datiF[:,1],'r-',label="drain/fracture PoV")
plt.plot(datiF[:,5], datiF[:,1],'b-',label="barriere/matrice PoV")
plt.legend()


plt.subplot(122)
plt.title("Capillary pressures in cells")
plt.grid(linewidth=1)
plt.plot(datC[:,0], datC[:,1],'r-',label="drain/frac")
plt.plot(datC[:,0], datC[:,3],'b-',label="barriere/matrice")
plt.legend()

plt.show()
