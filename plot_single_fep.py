import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np


def extract_data(input):
    data=[]
    
    
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile)
    cv=np.array([x[0] for x in file])
    fes=np.array([x[1] for x in file])
    err=np.array([x[2] for x in file])
    return(cv,fes,err)


def plot_single(cv,fes,err,file,labx):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    plt.plot(cv,fes,label="CV")
    plt.fill_between(cv,fes+err,fes-err,alpha=0.2)
    plt.xlabel(labx)
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,4.0])
    plt.ylim([-1.,70])
    plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig('{}_single.png'.format(file),format='png')

def main():
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    
    cmap_active=mpl.colormaps['rainbow']
    colors=[cmap_active(0),cmap_active(0.5),cmap_active(0.75),cmap_active(1.0)]

    cv,fes,err=extract_data(sys.argv[1])
    #cv1_state,cv2_state,free_grid_state,min_pt_state=extract_data(sys.argv[2])
    #print(cv1[min_pt[1][0]],cv2[min_pt[0][0]],free_grid[min_pt[0][0],min_pt[1][0]])
    file=os.path.splitext(sys.argv[1])[0]
    labx=sys.argv[2]
   
    plot_single(cv,fes,err,file,labx)
   

if __name__ == "__main__":
    main()