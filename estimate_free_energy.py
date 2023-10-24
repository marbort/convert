import sys

free_files=sys.argv[1:-1]
n_profiles=len(free_files)

import numpy as np

data_dicts={}
for i, free_file in enumerate(free_files):
    data_dicts[i]=np.loadtxt(free_file)

profiles=[]
for key in data_dicts.keys():
    if key != 'time':
        profiles.append(data_dicts[key][:,1])

profiles=np.array(profiles)
# Compute mean_profile
mean_profile=np.mean(profiles,axis=0)
# Compute standard deviation
error_profile=np.std(profiles,axis=0)*1.96/np.sqrt(n_profiles)
print(error_profile)
x=data_dicts[0][:,0]
print(x)

import matplotlib.pyplot as plt
# set style to seaborn colorblind
plt.style.use('seaborn-colorblind')

# set up figure
fig, ax = plt.subplots(1, 1, figsize=(9, 9),dpi=150)
plt.plot(x, mean_profile, color='black', linewidth=2)
plt.fill_between(x, mean_profile-error_profile, mean_profile+error_profile, color='black', alpha=0.2)
plt.xlabel('Time (ns)')
plt.ylabel('Free Energy (kJ/mol)')
plt.tight_layout()
plt.savefig('profile.png', dpi=300)
