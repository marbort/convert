import json
import os
from re import T
import numpy as np
import argparse
import dpdata
boxl=["A","B","C"]
name="THF_iPrMgCl_CP2K"
data=dpdata.LabeledSystem('./dpdata','deepmd/npy')



with open('cp2k_template.tmpl','r') as ifile:
    lines=ifile.read()
    
    for i,number in enumerate(data['coords']):
        print(i)
        crds=[]
        for j,item in enumerate(number):
                X=item[0]
                Y=item[1]
                Z=item[2]
                crds.append(data['atom_names'][data['atom_types'][j]]+\
                    "\t{:8.5f}\t{:8.5f}\t{:8.5f}".format(X,Y,Z))
                join='\n'.join(crds)
        struct=lines.replace("##coord##",join)
        box=[]
        for k,entry in enumerate(data['cells'][i]):
            X=entry[0]
            Y=entry[1]
            Z=entry[2]
            box.append(boxl[k]+"\t{:8.5f}\t{:8.5f}\t{:8.5f}".format(X,Y,Z))
            bjoin='\n'.join(box)
        struct=struct.replace("##cell##",bjoin)
        struct=struct.replace("##project##","{}_{}".format(name,i))

        if os.path.isdir("inputs"):
            print("Creating input for frame {}".format(i))
        else:
            os.mkdir("inputs")
            print("Creating input for frame {}".format(i))
        with open('inputs/cp2k_{}_{}.cinp'.format(name,i),'w') as ofile:
            ofile.write(struct)
        





