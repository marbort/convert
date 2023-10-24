import json
import os
from re import T
import numpy as np
import argparse
import shutil
import json
import dpdata as dp

parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--input' , dest='input',help='lammps trajectory')
parser.add_argument('--types_elm', dest='ty2elm', default=1,
                    type=str, help='types to element dictionary')
parser.add_argument('--save', dest='save', default=None,
                    type=str, help='Save trajectory to xyz or pdb')
parser.add_argument('--last', dest='last', default=False,action='store_true',
                    help='Save only last frame')



args = parser.parse_args()



def change_types(input,type_elm):
    with open(type_elm,'r') as dictfile:
        type_elm_dict=json.load(dictfile)
    data=dp.System(input,'dump')
    for i,item in enumerate(data["atom_names"]):
        data["atom_names"][i]=type_elm_dict[item]
    data.to('deepmd/npy','dpdata_elm')
    
def save_geo(dpdata,format,last):
    data=dp.System(dpdata,'deepmd/npy')
    if last:
        if format=='xyz':
            data.sub_system(-1).to(format,'last.xyz')
        elif format=='pdb':
            data.sub_system(-1).data.to(format,'last.pdb')
    else:
        if format=='xyz':
            data.to(format,'traj.xyz')
        elif format=='pdb':
            data.to(format,'traj.pdb')
        
        


change_types(args.input,args.ty2elm)
if args.save:
    save_geo('dpdata_elm',args.save,args.last)
        
        


        
    
        
    

        
    

