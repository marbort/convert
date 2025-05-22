import json
import dpdata
import numpy as np
from copy import deepcopy
import os

np.random.seed(1)
with open('Fit.json') as f:
    deepmd_dict = json.load(f)
training_set_dirs = deepcopy(deepmd_dict['training']['training_data']['systems'])
#training_set_dirs+= deepcopy(deepmd_dict['training']['validation_data']['systems'])
def train_val_split(data_dir,val_frac=0.05):
    data = dpdata.LabeledSystem(data_dir, fmt='deepmd/npy')
    numframes=len(data)

    indices = list(range(numframes))
    split = int(np.round(numframes * (1-val_frac)))
    train_indices = sorted(np.random.choice(indices, size=split, replace=False))
    val_indices = sorted(list(set(indices) - set(train_indices)))

    if len(train_indices)>0:
        train_data = data[train_indices]
    else:
        train_data = None
    if len(val_indices)>0:
        val_data = data[val_indices]
    else:
        val_data = None
    return train_data, val_data
common_name='/cluster/projects/nn4654k/sigbjobo/projects/ML-Cascella/'
curr_dir=os.getcwd()
deepmd_output_dict = deepcopy(deepmd_dict)
deepmd_output_dict['training']['training_data']['systems'] = []
deepmd_output_dict['training']['validation_data']={'systems':[]}
n=0
for f in deepcopy(training_set_dirs):
    folder=deepcopy(f)
    #name=folder.replace(common_name,'')
    it=os.path.basename(os.path.dirname(folder))
    name=os.path.basename(folder)
    print(n)
    train_name=curr_dir+'/training_set/'+str(n)+'_'+it+'_'+name
    val_name=curr_dir+'/validation_set/'+str(n)+'_'+it+'_'+name
    train_data, val_data = train_val_split(folder)
    if train_data is not None:
        deepmd_output_dict['training']['training_data']['systems'].append(train_name)
        train_data.to_deepmd_npy('./training_set/'+str(n)+'_'+it+'_'+name)
    if val_data is not None:
        deepmd_output_dict['training']['validation_data']['systems'].append(val_name)
        val_data.to_deepmd_npy('./validation_set/'+str(n)+'_'+it+'_'+name)
    del train_data, val_data
    n+=1
with open('Fit_split.json', 'w') as f:
    json.dump(deepmd_output_dict, f, indent=4)

