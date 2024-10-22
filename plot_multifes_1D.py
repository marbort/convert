import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np
from scipy.optimize import minimize
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def get_fes_limited(paths,max):
    fes_dict={}
    biases=[]
    for path in paths:
        fes=np.loadtxt(path,unpack=True)
        fes[1][fes[1]>max]=max
        bias=np.exp(fes[1]/2.34)
        fes_dict[path]=fes
        biases.append(bias)
    return(fes_dict,biases)

paths=glob.glob("*[0-9].dat")
fes_dict,biases=get_fes_limited(paths,25)
font = {'family' : 'Formular',
    'weight' : 'normal',
    'size'   : 46}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3
fig=plt.figure(figsize=(20,20),dpi=150)
for i,item in enumerate(fes_dict):
    #plt.plot(fes_dict[item][0],biases[i],label=item)
    plt.plot(fes_dict[item][0],fes_dict[item][1],label=item)
plt.legend()
plt.savefig('test.png',format='png',dpi=150)


    