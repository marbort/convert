import dpdata as dp
import glob
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import regex as re
from natsort import natsorted, ns
import json
import argparse


def change_types(input,type_elm):
    with open(type_elm,'r') as dictfile:
        type_elm_dict=json.load(dictfile)
    data=dp.System(input,'dump')
    for i,item in enumerate(data["atom_names"]):
        data["atom_names"][i]=type_elm_dict[item]
    return(data)


def lmp_to_xyzbox(input,output,ty2elm):
    trj="trj.xyz"
    data=change_types(input,ty2elm)
    box=[[str(x[i][i]) for i in range(3)] for x in data['cells']]
    elms=np.array([data['atom_names'][x] for x in data['atom_types']])
    """
    with open(output,'w') as ofile:
        for i,crd in enumerate(data['coords']):
            ofile.write(f"{len(crd)}\n")
            ofile.write(" ".join(box[i])+"\n")
            #for line in crd:
                #ofile.write(" ".join(line))
            crd_with_elm=np.array(np.c_[elms,crd],dtype='|S2,f8,f8,f8')
            #crd_with_elm=np.c_[elms,crd]
            #print(crd_with_elm[0])
            np.savetxt(ofile,crd_with_elm,fmt='%3s %.6f %.6f %.6f')
        
    
    """
    data.to('xyz',trj)

    with open(trj,'r') as ifile:
        lines=ifile.readlines()
        k=0
        for i,line in enumerate(lines):
            if line =='\n':
                lines[i]=" ".join(box[k])+"\n"
                k=k+1
    with open(output,'w') as ofile:
            for line in lines:
                ofile.write(line)
    subprocess.call(['rm',trj])
    
def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(
        description='Convert DPData files to a single extended xyz formatted file')
    parser.add_argument('--input',  type=str, help='Input lammps trj')
    parser.add_argument('--output', type=str,
                        help='Output xyz file')
    parser.add_argument('--types_elm', dest='ty2elm', default=1,
                    type=str, help='types to element dictionary')
    args = parser.parse_args()
    
    
    lmp_to_xyzbox(args.input,args.output,args.ty2elm)


if __name__=="__main__":
    main()
    