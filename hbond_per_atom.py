import numpy as np
import glob
import sys
import os


atom_numbs={"CHL":400,"CL":400,"GCL":800,"ACP":26,"Clm":400}

def calc_avg_per_atom(atoms_numb,hbonds,atom):
    hb=np.loadtxt(hbonds,unpack=True)
    avg=np.average([x/atoms_numb[atom] for x in hb[1]])
    std=np.std([x/atoms_numb[atom] for x in hb[1]])
    return(avg,std)


files=glob.glob("*hbond*")

for file in files:
    atom=os.path.basename(file).split('_')[3]
    avg,std=calc_avg_per_atom(atom_numbs,file,atom)
    print(file,f"{avg:.3f}",f"{std:.3f}")
    