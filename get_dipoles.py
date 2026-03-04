import sys
import numpy as np
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

start=False
results = []
dipoles = []
for line in lines:
    if "@" in line:
       start=False
    if "1\\1\\" in line:
       start=True
    if start:
       results.append(line.strip())
res_string="".join([x for x in results])
res_split=res_string.split("\\")
for i,j in enumerate(res_split):
   if "Dipole=" in j:
       try:
          vals=j.split("=")[1].split(",")
          dipoles.append([float(vals[0]),float(vals[1]),float(vals[2])])
       except:
          print(vals)
          print(j)
dipole_magnitudes=[np.linalg.norm(np.array(x)) for x in dipoles]
print(dipoles)
print(f"Dipole Magnitudes: {dipole_magnitudes}")
avg_dipole=np.mean(np.array(dipoles),axis=0)
print(f"Average Dipole Moment: {avg_dipole}")     
