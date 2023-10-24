import os
import sys
import numpy as np


def get_dipole_scan(file):
    dip=[]
    dip_conv=[]
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            if "Stationary" in line:
                dip_conv.append(dip[-1])
            if "Dipole  " in line:
                dip.append(line)
    dip_mag_ha=[np.sqrt(sum([float(j.split('=')[-1].rstrip()[i:i+15].replace('D','E')) for i in range(0, len(j.split('=')[-1].rstrip()), 15)])**2) for j in dip_conv ]
    dip_mag_d=[x/0.393456 for x in dip_mag_ha]
    return(dip_conv,dip_mag_ha,dip_mag_d)


dip_conv,dip_mag_ha,dip_mag_d=get_dipole_scan(sys.argv[1])
print(dip_mag_ha)
print(dip_mag_d)

#for j in dip_conv:
 #   y=j.split('=')[-1].rstrip()
  #  print([y[i:i+15].replace('D','E') for i in range(0, len(y), 15) ])