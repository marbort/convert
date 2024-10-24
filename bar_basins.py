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
import math


def bar_percentages(inputs):
    fig=plt.figure(figsize=(48,16),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 60}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    shades=["#1c1c1c","#ffd7d7","#ffd7d7","#99ccff"]
    bars={}
    x_vals=[]
    structs=["R=Me","R=Et","R=$i$Pr","R=$t$Bu"]
    titles=["(1,1)","(1,2)","(2,1)","(2,2)"]
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    for input in inputs:
        path=os.path.join(input.rstrip(),"basin_analysis.txt")
        structs.append(path.split('/')[1])
        with open(path.rstrip(),'r') as ifile:
            lines=ifile.readlines()    
            for i,line in enumerate(lines):
               
                if "Composition" in line:
                    basin=line.split()[4]+line.split()[5]
                    print(basin)
                    bars_tmp=[]
                    for j in range(5):
                        x_vals.append(lines[i+j+1].split(':')[0].split()[-1])
                        bars_tmp.append(float(lines[i+j+1].split(':')[-1].rstrip()))
                        bars_tmp=[0 if math.isnan(x) else x for x in bars_tmp]
                    if basin in list(bars.keys()):
                        bars[basin].append(bars_tmp)
                    else:
                        bars[basin]=[bars_tmp]
    x_vals=np.unique(x_vals)
   
    width=0.2
    print(len(bars))
    for k,bar in enumerate(bars):
        """
        if k==2:
            continue
        elif k==3:
            plt.subplot(1,len(bars),k)
            multiplier=0
            for j,mol in enumerate(bars[bar]):
                offset=width*multiplier
                plt.bar([float(x)+offset for x in x_vals[1:]],mol[1:],width=width,label=structs[j])
                multiplier+=1
                xtix=[float(x)+width*len(bars[list(bars.keys())[0]][0])/2 for x in x_vals[1:]]
                plt.xticks(xtix,labels=[0,1,2])
                plt.xlabel("# of Bridging Cl")
                plt.ylim([0,1.1])
            plt.title(titles[k-1])
        else:
        """
        plt.subplot(1,4,k+1)
        ax=fig.gca()
        #ax.set_facecolor(shades[k])
        ax.patch.set_alpha(0.5)
        multiplier=0
        for j,mol in enumerate(bars[bar]):
            offset=width*multiplier
            plt.bar([float(x)+offset for x in x_vals],mol,width=width,label=structs[j],color=pres_colors[j])
            multiplier+=1
        xtix=[float(x)+(width*(len(bars)-1))/2 for x in x_vals]
        print(xtix)
        plt.xticks(xtix,labels=[x-2 for x in range(len(bars)+1)])
        plt.xlabel("# of Bridging Cl")
        plt.xlim([1.5,5])
        plt.ylim([0,1.1])
        print(k)
        plt.title(titles[k])
    plt.subplot(1,4,1)
    plt.ylabel("Molecule fraction")
    #plt.subplot(1,len(bars),len(bars))
    plt.legend()
    plt.tight_layout()
    plt.savefig('bars.png')
            
                        
                    
                    
        
    
    

def main():
    with open('basin_analysis_paths.txt','r') as ifile:
        paths=ifile.readlines()
    bar_percentages(paths)
    
if __name__=="__main__":
    main()
        