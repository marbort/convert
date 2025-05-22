import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse

import numpy as np


def extract_data(input,zmax):
    data=[]
    
    
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    #print(file)
    file[1][file[1]>float(zmax)]=float(zmax)
    return(file[0],file[1],file[2])


def plot_single(cv,fes,err,file,labx,fig,zmax=100,labels=None,idx=0):
    
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    fig=fig
    plt.plot(cv,fes,label=labels[idx],color=pres_colors[idx],linewidth=3)
    plt.fill_between(cv,fes+err,fes-err,alpha=0.2,color=pres_colors[idx])
    plt.xlabel(labx)
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,4.0])
    plt.ylim([-1.,float(zmax)])
    if labels is not None:
        plt.legend(loc='upper right')
    #plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    

def main():
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    
    cmap_active=mpl.colormaps['rainbow']
    colors=[cmap_active(0),cmap_active(0.5),cmap_active(0.75),cmap_active(1.0)]
    
    parser = argparse.ArgumentParser(description='Plot free energy surface')
    parser.add_argument('--input', type=str, help='Input file',default=None)
    parser.add_argument('--paths', type=str, help='File with paths to inputs',default='paths.txt')
    parser.add_argument('--labx', type=str, help='Label for x-axis',default='CV1')
    parser.add_argument('--zmax', type=float, help='Maximum z value',default=100)
    parser.add_argument('--labels', type=str, help='Labels for the curves',default="FES", nargs='+')
    args = parser.parse_args()
    
    
    fig=plt.figure(figsize=(16,10),dpi=150)
    try:
        files=glob.glob(args.input)
    except:
        try:
            print("No input file found, reading from paths.txt")
            with open(args.paths,'r') as ifile:
                files=ifile.readlines()
                files=[x.strip() for x in files]
            print(files)
        except:
            print("No input file found, exiting")
            exit()
        
    for i,file in enumerate(files):
        cv,fes,err=extract_data(file,args.zmax)
        plot_single(cv,fes,err,file,args.labx,fig,args.zmax,args.labels,i)
    plt.savefig('fes_single.png',format='png')
        
        
        
        
   

if __name__ == "__main__":
    main()