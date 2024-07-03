import glob
import sys
import os
import dpdata as dp
"""
HW
OW
c
cc
cd
h4
ha
hn
na
o
"""

trajectories=glob.glob(sys.argv[1])
topo = sys.argv[2]
new_type=[0,1,2,2,2,0,0,0,3,1]
new_elem=["H","O","C","C","C","H","H","H","N","O"]


for file in trajectories:
    name=os.path.splitext(file)[0]
    data=dp.System(name,'amber/md',parm7_file=topo)
    data.to('deepmd/npy',f'dpdata_{name}')

    
    


for folder in glob.glob("dpdata*"):
    print(folder)
    newlines=[]
    with open(f"{folder}/type.raw",'r') as ifile:
        lines=ifile.readlines()
    for line in lines:
        tip=line.split()[0]
        newlines.append(new_type[int(tip)])
    with open(f"{folder}/type.raw",'w') as ofile:
        for line in newlines:
            ofile.write(f'{line}\n')
    with open(f"{folder}/type_map.raw",'r') as ifile:
        lines=ifile.readlines()
    with open(f"{folder}/type_map.raw",'w') as ofile:
        for i,line in enumerate(lines):
            ofile.write(f'{new_elem[i]}\n')
    
        
        
        

        
        