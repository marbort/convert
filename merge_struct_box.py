import os
import dpdata
import numpy as np
import argparse
import re
import subprocess

mols=["PhMgCl_dimer","Ph-Cl-bridge","Ph-bridge","PhMgCl","Ph2Mg"]
box={"PhMgCl_dimer":[30.080,30.080,30.080],"Ph-Cl-bridge":[30.080,30.080,30.080],"Ph-bridge":[30.080,30.080,30.080],
     "PhMgCl":[17.58,17.58,17.58],"Ph2Mg":[17.58,17.58,17.58]}
root="/run/media/marco/SHARED/RATIO/WP1/ML/fox/PheMgCl/"
for mol in mols:
    with open(root+"BOX_"+mol+".xyz",'r') as ifile:
        lines=ifile.readlines()
        size_box=int(lines[0])
        crds_box=[x for x in lines[2:]]
    with open(root+mol+"_only.xyz",'r') as ifile:
        lines=ifile.readlines()
        size_mol=int(lines[0])
        crds_mol=[x for x in lines[2:]]
    with open(root+mol+"_BOX.xyz",'w') as ofile:
        ofile.write("{}\nCreated by Me\n".format(size_box+size_mol))
        for line in crds_mol:
            ofile.write(line)
        for line in crds_box:
            ofile.write(line)
     #Read type_map.raw 
    
    data=dpdata.System(root+mol+"_BOX.xyz",'xyz')
    for i,x in enumerate(data['cells'][0]):
        data['cells'][0][i][i]=box[mol][i]
    data.to('deepmd/npy',root+mol+"_BOX_tset")

    type_map_dict = {}
    with open(os.path.join(root+mol+'_BOX_tset','type_map.raw'), 'r') as f:
        for i, line in enumerate(f):
            type_map_dict[str(i)] = line.strip()

    # Read type.raw
    types = np.loadtxt(os.path.join(root+mol+'_BOX_tset','type.raw'), dtype=int)
    

    # Read new type_map.raw
    type_map_dict_new = {}
    with open(root+'type_map.raw', 'r') as f:
        for i, line in enumerate(f):
            type_map_dict_new[line.strip()] = i

    # Create new type.raw
    types_new = np.zeros_like(types)
    for i, t in enumerate(types):
        types_new[i] = type_map_dict_new[type_map_dict[str(t)]]

    # Write new type.raw and type_map.raw
    np.savetxt(os.path.join(root+mol+'_BOX_tset', 'type.raw'), types_new, fmt='%d')
    np.savetxt(os.path.join(root+mol+'_BOX_tset', 'type_map.raw'),
            np.array(list(type_map_dict_new.keys())), fmt='%s')
    
    subprocess.call(["rm","-f",root+mol+'_BOX_tset/nopbc'])

    
    data=dpdata.System(root+mol+"_BOX_tset",'deepmd/npy')
    # Convert to lammps data file
    nframes = len(data['coords'])
    frames = [0]
    name = "conf.data"
    data.to_lammps_lmp(os.path.join(root+mol+'_BOX_tset',name))

    masses = [dpdata.periodic_table.Element(
        name).mass for name in data['atom_names']]
    text = " Masses\n\n"
    for i, m in enumerate(masses):
        text += "{} {}\n".format(i+1, m)
    text += '\n'
    text += ' Atoms'
    with open(os.path.join(root+mol+'_BOX_tset',name), 'r') as fp:
        content = fp.read()
        content = re.sub("Atoms", text, content)
    with open(os.path.join(root+mol+'_BOX_tset',name), 'w') as fp:
        fp.write(content)
        
    