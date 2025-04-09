import json
import argparse
import os
import dpdata as dp
import numpy as np

def add_new_atom_type_lammps(position_dict,tset,subst,name,output):
    pairs_dict = json.load(open(position_dict))
    data=dp.System(tset,'lammps/lmp')
    crd_len=str(len(data['coords'][0]))
    #subst_idx=data['atom_names'].index(subst)
    allgood=False
    data['atom_names'].append(name)
    new_idx=data['atom_names'].index(name)

    if crd_len in pairs_dict.keys():
        changed_idxs=[]
        for i in range(len(pairs_dict[crd_len])):
            allgood=[]
            for j in pairs_dict[crd_len][i]:
                if int(data['atom_types'][j]) == int(subst.split('_')[-1]):
                    allgood.append(True)
                else:
                    allgood.append(False)
            if np.all(allgood):
                for j in pairs_dict[crd_len][i]:
                    data['atom_types'][j]=new_idx
                    changed_idxs.append(j)
                assert changed_idxs == pairs_dict[crd_len][i] , f"Error in {tset} with {crd_len} atoms. Changed atoms {changed_idxs} but should be {pairs_dict[crd_len][i]}"
        print(f"Changed {changed_idxs} atoms")
        data.to('lammps/lmp',os.path.join(output,"temp.data"))
        return(os.path.join(output,"temp.data"),len(data['atom_names']))
    else:
        print("Not found length {crd_len} in the dictionary")
    


def add_masses_lammps(old_datafile,new_datafile,mass,new_index,temp_datafile,total_types):
    with open(old_datafile,'r') as ifile:
        lines=ifile.readlines()
    masses=[]
    lines_with_masses=[]
    start=False
    for line in lines:
        if "Atoms"  in line:
            start=False
        if "Masses" in line:
            start=True
        if start:
            masses.append(line)
    masses.pop(-1)
    masses.append(f"{int(new_index.split('_')[-1])+1} {mass}\n\n")
   
    with open(temp_datafile,'r') as  tempfile:
        newlines=tempfile.readlines()
        for line in newlines:
            if "atom types" in line:
                lines_with_masses.append(line.replace(f"{total_types-1}",f"{total_types}"))
            elif "Atoms" in line:
                masses.append(line)
                lines_with_masses.append(line.replace(line,"".join(masses)))
            else:
                lines_with_masses.append(line)
    with open(os.path.join(new_datafile,f"conf_{new_index}.data"),'w') as ofile:
        for line in lines_with_masses:    
            ofile.write(line)
        
                
        
            
    


def main():
    parser = argparse.ArgumentParser(description='Add new atom type to the tset')
    
    parser.add_argument('--dict', type=str, help='Dictionary of positions')
    parser.add_argument('--datafile', type=str, help='LAMMPS datafile')
    parser.add_argument('--out', type=str, help='output path for the new datafile',default='.')    
    parser.add_argument('--newindex', type=str, help='New atom type index e.g. Type_0, Type_1... etc')
    parser.add_argument('--subst', type=str, help='Old atom type index to substitute e.g. Type_0, Type_1... etc')
    parser.add_argument('--mass' , type=float, help='Mass of the new atom type')
    parser.add_argument('--fromjson', type=str, help='Take the tsets from a json file',default=None)
    
    args=parser.parse_args()
    
    if args.fromjson:
        dct=json.load(open(args.fromjson))
        paths=dct['training']['training_data']['systems']
        for path in paths:
            temp_datafile,total_types=add_new_atom_type_lammps(args.dict,path,args.subst,args.name,args.out)
            add_masses_lammps(args.datafile,args.out,args.mass,args.newindex,temp_datafile)
    
    else:
        temp_datafile,total_types=add_new_atom_type_lammps(args.dict,args.datafile,args.subst,args.newindex,args.out)
        add_masses_lammps(args.datafile,args.out,args.mass,args.newindex,temp_datafile,total_types)
    


if __name__ == "__main__":
    main()    
    
   