import json
import os
from re import T
import numpy as np
import argparse
import dpdata
boxl=["A","B","C"]
name="THF_iPrMgCl_CP2K"
data=dpdata.LabeledSystem('./dpdata','deepmd/npy')


parser = argparse.ArgumentParser(description='Plot data')
parser.add_argument('--offset', dest='offset', default=1, type=int, help='save only every nth frame')
parser.add_argument('--name', dest='name', default=1, type=str, help='Project name')
parser.add_argument('--template', dest='template', default='cp2k_template.tmpl', type=str, help='File path to template file.')

args = parser.parse_args()


with open(args.template,'r') as ifile:
    lines=ifile.read()
    
    for i in range(0,len(data['coords']),args.offset):
        print(i)
        crds=[]
        for j,item in enumerate(data['coords'][i]):
                X=item[0]
                Y=item[1]
                Z=item[2]
                crds.append("       "+data['atom_names'][data['atom_types'][j]]+\
                    "   {:13.9f}   {:13.9f}   {:13.9f}".format(X,Y,Z))
                join='\n'.join(crds)
        struct=lines.replace("##coord##",join)
        box=[]
        for k,entry in enumerate(data['cells'][i]):
            X=entry[0]
            Y=entry[1]
            Z=entry[2]
            box.append("       " + boxl[k]+"   {:13.9f}   {:13.9f}   {:13.9f}".format(X,Y,Z))
            bjoin='\n'.join(box)
        struct=struct.replace("##cell##",bjoin)
        struct=struct.replace("##project##","{}_{:06d}".format(args.name,i))

        if os.path.isdir("inputs"):
            print("Creating input for frame {}".format(i))
        else:
            os.mkdir("inputs")
            print("Creating input for frame {}".format(i))
        with open('inputs/cp2k_{}_{:06d}.cinp'.format(args.name,i),'w') as ofile:
            ofile.write(struct)
        





