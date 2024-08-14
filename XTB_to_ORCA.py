import numpy as np
import sys

with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()
with open(sys.argv[2],'r') as template:
    route=template.readlines()

natoms=int(lines[0].strip())
crds=[]
i=0
while True:
    try:
        tmp=[lines[x] for x in range(i+2,i+natoms+2)]
        crds.append(tmp)
        i+=natoms+2
    except:
        break

for i,struct in enumerate(crds):
    with open(f'orca_{i:06d}.inp','w') as ofile:
        for line in route:
            ofile.write(line)
        for line in struct:
            ofile.write(line)
        ofile.write('*\n')




