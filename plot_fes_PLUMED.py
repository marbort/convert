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



def extract_data(input):
    data=[]
    
    
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    #cv1_temp=np.array([x[0] for x in file])
    #cv2_temp=np.array([x[1] for x in file])
    cv1=np.unique(file[0])
    cv2=np.unique(file[1])
    val=np.array(file[2])
    val_shift=val-min(val)
    free_grid=val_shift.reshape(len(cv1),len(cv2))
    
    return(cv1,cv2,free_grid)

def plot2d(x,y,maxz,value,file,labx,laby,cmap,minima):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 30}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX=int(maxz)
    
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #kjmol/plot
    lev=range(0,MAX+5,5)
    
    CLines=plt.contour(x, y,value,levels=range(0,MAX,20),vmin=0,vmax=MAX,linewidths=1,colors='black')
    plt.clabel(CLines,levels=range(0,MAX,20), inline=True, fontsize=10,colors='black')
    plt.contourf(x, y,value,lev,vmin=0,vmax=MAX,cmap=cmap)
    ###
    
    #kcal/mol plot
    #lev=[x/2 for x in range(0,21,1)]
    #val_kcal=[x/4.184 for x in value]
    #plt.contourf(x, y,val_kcal,lev,vmin=0,vmax=10,cmap=cmap)
    
    plt.xlabel(labx)
    plt.ylabel(laby)
    plt.xticks(np.arange(min(x),max(x)+0.5,0.5))
    #bounds=[1,2,3,4]
    #cbarkcal
    #cbar=plt.colorbar(label="$\Delta A\ (kcal\ mol^{-1})$",ticks=range(0,11,1))
    #cbar.ax.set_ylim(0,10)
    #cbar kj/mol
    cbar=plt.colorbar(label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(0,MAX+20,20))
    cbar.ax.set_ylim(0,MAX)
    #for i in minpath:
    #    plt.scatter(x[i[0]],y[i[1]],color='black')
    #plt.scatter(x[minima[0]],y[minima[1]])
    #plt.scatter(x[maxima[0]],y[maxima[1]],color='red')
    #plt.xlim([0.75,3.25])
    #plt.ylim([0.75,3.25])
    if len(minima) > 0:
        min_crd=[]
        with open(minima,'r') as ifile:
            lines=ifile.readlines()
        for line in lines:
            min_crd.append((int(line.split()[1]),float(line.split()[5]),float(line.split()[-1])))
        for pt in min_crd:
            plt.scatter(pt[1],pt[2],color='white')
    good_minima=[]
    try:
        for i,item in enumerate(min_pt[0]):
            if value[item,min_pt[1][i]] > 60:
                pass
            else:
                good_minima.append((x[min_pt[1][i]],y[item],value[item][min_pt[1][i]]))
                #plt.scatter(x[min_pt[1][i]],y[item],color='red')
        with open('minima.dat','w') as ofile:
            for i in good_minima:
                ofile.write(" ".join([f"{j:10.4f}" for j in i])+"\n")
    except:
        print("No minima found. Skipping")
        pass
        
    plt.savefig('{}_PLUMED.png'.format(file),format='png')

def main():

    cv1,cv2,free_grid=extract_data(sys.argv[1])
    file=os.path.splitext(sys.argv[1])[0]
    labx=sys.argv[2]
    laby=sys.argv[3]
    MAX=sys.argv[4]
    try:
        minima=sys.argv[5]
    except:
         minima=""
    cmap_active='rainbow'
    plot2d(cv1,cv2,MAX,free_grid,file,labx,laby,cmap_active,minima)
   

if __name__ == "__main__":
    main()