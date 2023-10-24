import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import argrelmin
from scipy.signal import find_peaks
     

def coord_Number(x,y,dens,dist_min,dist_max,mols):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*dens*np.trapz(prod,val_x)
    CN_one=CN_all/mols
    return(CN_one)



with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()

with open(sys.argv[2],'r') as ifile2:
    lines2=ifile2.readlines()

x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
y=[float(x.split()[1]) for x in lines if "#" not in x if "@" not in x]

x2=[float(x.split()[0]) for x in lines2 if "#" not in x if "@" not in x]
y2=[float(x.split()[1]) for x in lines2 if "#" not in x if "@" not in x]

grad=np.gradient(y2)

minima=find_peaks([-x for x in grad],width=2)

print(minima)
CN=coord_Number(x,y,0.003495,0,4.5,26)
CN_from_int=y2[minima[0][0]]
print(CN,CN_from_int)

plt.plot(y2,grad)
plt.scatter([y2[x] for x in minima[0]],[grad[x] for x in minima[0]])
plt.show()