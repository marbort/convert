import numpy as np
import json

with open('sys_spec.json','r') as ifile:
    sys_spec=json.load(ifile) 

input_crd=np.loadtxt('It1_000000-HOMO_centers_s1-1_0.xyz',skiprows=2,usecols=(1,2,3))
input_atoms=np.loadtxt('It1_000000-HOMO_centers_s1-1_0.xyz',skiprows=2,usecols=(0),dtype=str)

elm_atoms=[x for x in input_atoms if x in sys_spec["elements"]]
atoms=[input_crd[i] for i,x in enumerate(input_atoms) if x in sys_spec["elements"]]
wc=[input_crd[i] for i,x in enumerate(input_atoms) if x[0] == "X"]


mols=[[] for x in atoms]
#print(mols)
for center in wc:
    dists=[]
    for i,atom in enumerate(atoms):
        dist=np.linalg.norm(center-atom)
        #print(dist)
        dists.append((dist,i))
    closest=min(dists)
    #print(closest)
    mols[closest[1]].append(center)
#print(mols[1])
for k,mol in enumerate(atoms):
    centroid_x=sum([x[0] for x in mols[k]])/len(mols[k])
    centroid_y=sum([x[1] for x in mols[k]])/len(mols[k])
    centroid_z=sum([x[2] for x in mols[k]])/len(mols[k])
    
    
    with open(f"{elm_atoms[k]}{k}_WC.xyz",'w') as ofile:
        ofile.write(f"{len(mols[k])}\n")
        ofile.write(f"Wannier Centers for {elm_atoms[k]}{k}\n")
        for  ctr  in mols[k]:
            ofile.write("X   " + "   ".join([str(x) for x in ctr])  + "\n")
    with open(f"{elm_atoms[k]}{k}_COM_WC.xyz",'w') as ofile:
        ofile.write(f"1\n")
        ofile.write(f"Centroid of Wannier Centers for {elm_atoms[k]}{k}\n")
        ofile.write(f"X  {centroid_x}  {centroid_y}  {centroid_z}\n")
        

#print(mols[1])
            
        
    
        



#atoms=[x for x in input if ]

