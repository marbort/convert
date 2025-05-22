import dpdata as dp
import json
import argparse
import os
import numpy as np


def add_new_atom_type(position_dict,tset,subst,name,output):
    pairs_dict = json.load(open(position_dict))
    data=dp.LabeledSystem(tset,'deepmd/npy')
    crd_len=str(len(data['coords'][0]))
    subst_idx=data['atom_names'].index(subst)
    allgood=False
    data['atom_names'].append(name)
    new_idx=data['atom_names'].index(name)
    print(crd_len)
    print(pairs_dict.keys())
    if crd_len in pairs_dict.keys():
        changed_idxs=[]
        for i in range(len(pairs_dict[crd_len])):
            allgood=[]
            for j in pairs_dict[crd_len][i]:
                if data['atom_types'][j] == subst_idx:
                    allgood.append(True)
                else:
                    allgood.append(False)
            if np.all(allgood):
                for j in pairs_dict[crd_len][i]:
                    data['atom_types'][j]=new_idx
                    changed_idxs.append(j)
                assert changed_idxs == pairs_dict[crd_len][i] , f"Error in {tset} with {crd_len} atoms. Changed atoms {changed_idxs} but should be {pairs_dict[crd_len][i]}"
        print(f"Changed {changed_idxs} atoms")
        data.to('deepmd/npy',os.path.join(output,f"{os.path.basename(tset)}_{name}"))
    else:
        print(f"Not found length {crd_len} in the dictionary")
    
def fix_type_map(tset,new_type_map):
    # Read type_map.raw
    type_map_dict = {}
    with open(os.path.join(tset, 'type_map.raw'), 'r') as f:
        for i, line in enumerate(f):
            type_map_dict[str(i)] = line.strip()

    # Read type.raw
    types = np.loadtxt(os.path.join(tset, 'type.raw'), dtype=int)


    # Read new type_map.raw
    type_map_dict_new = {}
    with open(new_type_map, 'r') as f:
        for i, line in enumerate(f):
            type_map_dict_new[str(i)] = line.strip()

    inv_type_map_dict_new = {v: k for k, v in type_map_dict_new.items()}
    # Create new type.raw
    types_new = np.zeros_like(types)
    for i, t in enumerate(types):
        types_new[i] = int(inv_type_map_dict_new[type_map_dict[str(t)]])

    # Write new type.raw and type_map.raw
    np.savetxt(os.path.join(tset, 'type.raw'), types_new, fmt='%d')
    np.savetxt(os.path.join(tset, 'type_map.raw'),
            np.array(list(type_map_dict_new.values())), fmt='%s')



def main():
    parser = argparse.ArgumentParser(description='Add new atom type to the tset')
    
    parser.add_argument('--dict', type=str, help='Dictionary of positions')
    parser.add_argument('--tset', type=str, help='Training set')
    parser.add_argument('--out', type=str, help='output path for the new tset',default='.')    
    parser.add_argument('--name', type=str, help='Name of the atom type')
    parser.add_argument('--subst', type=str, help='Old atom tye to substitute')
    parser.add_argument('--fromjson', type=str, help='Take the tsets from a json file',default=None)
    parser.add_argument('--type_map', type=str, help='Type map to fix eventual inconsistencies',default=None)
    
    args=parser.parse_args()
    
    if args.fromjson:
        dct=json.load(open(args.fromjson))
        paths=dct['training']['training_data']['systems']
        for path in paths:
            add_new_atom_type(args.dict,path,args.subst,args.name,args.out)
            fix_type_map(path,args.type_map)
    
    else:
        add_new_atom_type(args.dict,args.tset,args.subst,args.name,args.out)
        fix_type_map(args.tset,args.type_map)


if __name__ == "__main__":
    main()    
