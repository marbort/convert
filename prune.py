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
    test_frc=[]
    test_nrg=[]
    for i,item in enumerate(prop):
        print(thresh[i])
        if prop=='frc_tot':
            test_frc=[np.linalg.norm(x) < thresh[i] for x in data['forces'] ]
            print(test)
            forces_pruned=[data['forces'][k] for k in range(len(test)) if test[k]]
        if prop=='frc_per_atom':
            test=[[np.linalg.norm(m) < thresh[0] for m in l] for l in data['forces'] ] 
            test_frc=[x for x in test if all(x)]
            forces_pruned=[data['forces'][k] for k in range(len(test)) if test[k]]
        else:
            test_nrg=[x < thresh[i] for x in data['energies']]
    
    print(test_nrg)
    if not test_frc:
        test_final=test_nrg
    if not test_nrg:
        test_final=test_frc
    elif test_frc and test_nrg:
        print("CIAO2")
        test_final=test_nrg and test_frc
    print("Initial Number of Frames {}".format(len(data['energies'])))
    for i in prop:
        print("Frames pruned for {}  {}".format(i,len(data['energies'])-sum(test_final)))
    print("Final Number of Frames {}".format(sum(test_final)))
    filtered_indices=[i for i,item in enumerate(test_final) if item]

    tset_folder=path+"_pruned_auto/"
    isExist = os.path.exists(tset_folder)
    if not isExist:
        os.makedirs(tset_folder)
    
    pruned_dpdata = data.sub_system(filtered_indices)
    pruned_dpdata.to('deepmd/npy', tset_folder)

    

        





    
    
