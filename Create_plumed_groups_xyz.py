import glob
import os
import regex as re
import json
import argparse
import dpdata as dp
import numpy as np 

def create_plumed_input(name):
    with open('input_plmd.json') as ifile:
        input_dict=json.load(ifile)
    with open('plumed.dat','w') as pf:
        #print(tset,args.meta,previous_iteration)
        data=dp.System(name,'xyz')
        print(input_dict['plumed']['print'])
        input_dict['plumed']['print']=input_dict['plumed']['print'].replace('FILE=colvar',f'FILE=colvar_{name}')
        groups=[np.where(data['atom_types']==data['atom_names'].index(input_dict['plumed']['groups'][x][0])) for x in input_dict['plumed']['groups']]
        for i,item in enumerate(input_dict['plumed']['groups']):
            if input_dict['plumed']['groups'][item][1]=="all":
                sel_group=groups[i][0]
                sel_group=[x+1 for x in sel_group]
                sel_group_str=np.char.mod('%0d', sel_group)

                pf.write(item+": GROUP ATOMS="+','.join(sel_group_str)+'\n')
            else:
                sel_group=[groups[i][0][input_dict['plumed']['groups'][item][1][x]] for x in range(len(input_dict['plumed']['groups'][item][1]))]
                sel_group=[x+1 for x in sel_group]
                sel_group_str=np.char.mod('%0d', sel_group)
                pf.write(item+": GROUP ATOMS="+','.join(sel_group_str)+'\n')
            pf.write("\n\n")
        for i in input_dict['plumed']['cvs']:
                pf.write(i+'\n')
        pf.write("\n\n")
        
        pf.write(input_dict['plumed']['print']+'\n\n')

def main():
         parser = argparse.ArgumentParser()
         parser.add_argument('--input', type=str, dest='input', 
                    help='LAMMPS lmp file to create plumed file from')
         args = parser.parse_args()
        
         create_plumed_input(args.input)



if __name__== "__main__":
     main()
