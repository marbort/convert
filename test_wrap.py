import dpdata as dp
import numpy as np

box=[15 , 15, 15]
types=["C","Cl","H","O","Mg"]
data=dp.System('lmp.lammpstrj','lammps/dump')
#center
transl=data['coords'][0][0]-np.array(box)/2
for i in range(len(data['coords'])):
    data['coords'][i]-=transl
#wrap
for frame in range(len(data['coords'])):
    for atom in range(len(data['coords'][frame])):
        for dim in range(len(data['coords'][frame][atom])):
            if data['coords'][frame][atom][dim] > box[dim]:
                data['coords'][frame][atom][dim]-=box[dim]
            if data['coords'][frame][atom][dim] < 0:
                data['coords'][frame][atom][dim]=box[dim]+data['coords'][frame][atom][dim]

for i in range(len(data['atom_names'])):
    data['atom_names'][i]=types[i]


data.to('deepmd/npy','dpdata')
                
        