import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import argparse
import sys
import numpy as np
from scipy.optimize import minimize
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def extract_data(input):
    data=np.loadtxt(input,unpack=True)
    return(data)

def plot_lcurve(data):
    font = {'family' : 'Formular',
    'weight' : 'normal',
    'size'   : 46}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    fig=plt.figure(figsize=(30,10),dpi=150)
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    titles=["Total","Energy","Forces"]
    
    for i in range(3):
        plt.subplot(1,3,i+1)
        plt.plot(data[0],data[2*i+2],label="Training RMSE")
        plt.plot(data[0],data[2*i+1],label="Validation RMSE")
        plt.xlabel("Step Number")
        plt.legend()
        plt.title(titles[i])
    plt.tight_layout()
    
    plt.savefig("lcurve.png",format='png')

def main():
    input=sys.argv[1]
    data=extract_data(input)
    plot_lcurve(data)
    
if __name__=="__main__":
    main()