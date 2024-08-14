import glob
import dpdata
import os
import numpy as np

input_dict={"model_devi":{"max_force":100}}

for tset_folder in glob.glob('tset_*'):

            try:
                unfiltered_dpdata = dpdata.LabeledSystem(
                    tset_folder, fmt='deepmd/npy')
            except:
                print("Empty training set: {}.".format(tset_folder))
                print("removing {tset_folder}".format(tset_folder=tset_folder))
                os.system(
                    'rm -r {tset_folder}'.format(tset_folder=tset_folder))
                continue

            # Loop through all configurations, and keep track of the indices of configurations with forces less than input_dict["model_devi"]["max_force"]
            # Use the indices to filter out configurations with forces greater than input_dict["model_devi"]["max_force"]
            # Write the filtered dpdata to a new folder
            filtered_indices = []
            for i in range(len(unfiltered_dpdata['forces'])):
                if np.linalg.norm(unfiltered_dpdata['forces'][i], axis=1).max() < input_dict["model_devi"]["max_force"]:
                    filtered_indices.append(i)
            if len(filtered_indices) == 0:
                print("No configurations with forces less than {max_force} found in {tset_folder}".format(
                    max_force=input_dict["model_devi"]["max_force"], tset_folder=tset_folder))
                print("removing {tset_folder}".format(tset_folder=tset_folder))
                os.system(
                    'rm -r {tset_folder}'.format(tset_folder=tset_folder))
                continue

            filtered_dpdata = unfiltered_dpdata.sub_system(filtered_indices)
            filtered_dpdata.to('deepmd/npy', tset_folder)

