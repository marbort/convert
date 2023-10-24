#%%

import numpy as np
import matplotlib.pyplot as plt
import os
import glob

with open('/home/marco/dottorato/QMmodels/cluster/param/SeMET/dihed/MSE_fit_out_energy_genetic.dat','r') as ifile:
    lines=ifile.readlines()
    num=[float(x.split()[0]) for x in lines if "#" not in x]
    AMBERK=[float(x.split()[1]) for x in lines if "#" not in x]
    QM=[float(x.split()[2]) for x in lines if "#" not in x]
    
fig=plt.figure(figsize=(5,5),dpi=150)
plt.plot(num,AMBERK,num,QM)
plt.legend(["Amber+K","QM"])
plt.xlabel("Structure Number")
plt.ylabel("$Energy\ /\ kcal\ mol^{-1}$")

# %%
