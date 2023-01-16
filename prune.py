import json
import os
from re import T
import numpy as np
import argparse
import shutil
import dpdata as dpd

def prune_coords(prop: list,
                thresh: list,
                path: str):
    """
    Remove coordinates and relative data based
    on selected property

    Parameters
    ----------

    prop: list
    selected properties. Multiple selection
    will afford frames that satisfy ALL the 
    requested criteria:
    'energies'
    #'coords'
    'frc_tot'
    'frc_per_atom' 
    

    tresh: treshold for each selected property.
    Frames with a higher threshold will be deleted

    path: str
    path to the dpdata folder
    """
    data = dpd.LabeledSystem(path,'deepmd/npy',)
    print(data['energies'])
    
    #data_pruned=data.copy()
    test_frc=[]
    test_nrg=[]
    for i,item in enumerate(prop):
        print(thresh[i])
        if prop=='frc_tot':
            test_frc=[np.linalg.norm(x) < thresh[i] for x in data['forces'] ]
            print(test)
            forces_pruned=[data['forces'][k] for k in range(len(test)) if test[k]]
        if prop=='frc_per_atom':
            #test=[np.linalg.norm(sum(x)) < thresh[0] for i in range(len(data['forces']))  
             #for x in data['forces'][i][j]]
            test=[[np.linalg.norm(m) < thresh[0] for m in l] for l in data['forces'] ] 
            test_frc=[x for x in test if all(x)]
            #print(test_atom)
            forces_pruned=[data['forces'][k] for k in range(len(test)) if test[k]]
        else:
            test_nrg=[x < thresh[i] for x in data['energies']]
    
    print(test_nrg)
    if not test_frc:
        print("CIAO")
        print(test_nrg)
        test_final=test_nrg
    if not test_nrg:
        test_final=test_frc
    elif test_frc and test_nrg:
        print("CIAO2")
        test_final=test_nrg and test_frc
    print(test_final)
    print("Initial Number of Frames {}".format(len(data['energies'])))
    for i in prop:
        print("Frames pruned for {}  {}".format(i,len(data['energies'])-sum(test_final)))
    print("Final Number of Frames {}".format(sum(test_final)))
    print(test_final)
    energies_pruned=[data['energies'][k] for k in range(len(test_final)) if test_final[k]]
    forces_pruned=[data['forces'][k] for k in range(len(test_final)) if test_final[k]]
    box_pruned=[data['cells'][k] for k in range(len(test_final)) if test_final[k]]
    coord_pruned=[data['coords'][k] for k in range(len(test_final)) if test_final[k]]
    
        
    #print(forces_pruned[0])
    check=np.load(path+"/set.000/force.npy")
    #print(check[0])
    #print(data['forces'][0])
    isExist = os.path.exists(path+"_pruned/set.000")
    if not isExist:
        os.makedirs(path+"_pruned/set.000")
    
    np.save(path+"_pruned/set.000/force.npy",forces_pruned)
    np.save(path+"_pruned/set.000/energy.npy",energies_pruned)
    np.save(path+"_pruned/set.000/box.npy",box_pruned)
    np.save(path+"_pruned/set.000/coord.npy",coord_pruned)
    np.savetxt(path+"_pruned/type.raw",data['atom_types'],fmt='%s')
    np.savetxt(path+"_pruned/type_map.raw",data['atom_names'],fmt='%s')

        #data_pruned['forces']=[data['forces'][i] for i in range(len(test)) if test[i]]

        





    
    
