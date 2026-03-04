#calculate dipole dipole interaction in python
import sys
import numpy as np


dip_1=np.array([float(x) for x in sys.argv[1].split(',')])
dip_2=np.array([float(x) for x in sys.argv[2].split(',')])
r_vec=dip_2 - dip_1
r_mag=np.linalg.norm(r_vec)
r_hat=r_vec/r_mag
eps_0=8.854187817e-12 #C^2/(N m^2) 
k_e=1/(4*np.pi*eps_0) #N m^2/C^2

eps_0_au=55.263494E-3  #e^2/(eV Å) 
k_e_au=1/(4*np.pi*eps_0_au) # eV Å/e^2

U=k_e_au*(np.dot(dip_1,dip_2)/r_mag**3 - 3*(np.dot(dip_1,r_hat))*(np.dot(dip_2,r_hat))/r_mag**3)
print(f"Dipole-dipole interaction energy (eV): {U}")
print(f"Dipole-dipole interaction energy (kJ/mol): {U*96.485}")