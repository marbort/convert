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
from plot2d_reweight import extract_data

def plot2d(x,y,value,file,labx,laby,cmap,minima,min_pt):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 30}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX=50
    MIN=-50
    
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #kjmol/plot
    lev=range(MIN-2,MAX+2,2)
    plt.contourf(x, y,value,lev,vmin=MIN,vmax=MAX,cmap=cmap)
    ###
    
    #kcal/mol plot
    #lev=[x/2 for x in range(0,21,1)]
    #val_kcal=[x/4.184 for x in value]
    #plt.contourf(x, y,val_kcal,lev,vmin=0,vmax=10,cmap=cmap)
    
    plt.xlabel(labx)
    plt.ylabel(laby)
    #bounds=[1,2,3,4]
    #cbarkcal
    #cbar=plt.colorbar(label="$\Delta A\ (kcal\ mol^{-1})$",ticks=range(0,11,1))
    #cbar.ax.set_ylim(0,10)
    #cbar kj/mol
    cbar=plt.colorbar(label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(MIN-MIN//10,MAX+MAX//10,abs(MIN//10)))
    cbar.ax.set_ylim(MIN,MAX)
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
            if value[item,min_pt[1][i]] > 80:
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
        
    plt.savefig('{}.png'.format(file),format='png')


def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--ref', dest='ref', 
                        type=str, help='ref FES data',nargs="+")
    parser.add_argument('--reflab', dest='reflab', 
                        type=str, help='ref label names',nargs="+")
    parser.add_argument('--target', dest='target', 
                        type=str, help='target FES data',nargs="+")
    parser.add_argument('--targetlab', dest='targetlab', 
                        type=str, help='target label names',nargs="+")
    parser.add_argument('--fac', dest='fac', 
                        type=float, help='cv conversion factors',nargs="+")
    parser.add_argument('--xlabel', dest='xlabel', default='CV', type=str,
                         help='x label')
    parser.add_argument('--ylabel', dest='ylabel', default='Free Energy ($kJ\ mol^{-1}$)', type=str,
                         help='y label')
    parser.add_argument('--labels', dest='labels', default=None, type=str, nargs='+',
                         help='Plot labels label')
    parser.add_argument('--pres', dest='pres', action='store_true',
                         help='use presentation colors')
    parser.add_argument('--color', dest='color', type=str,
                         help='choose single color')
    args = parser.parse_args()
    
    
    cmap='jet'
    minima=""
    min_pt=""
    refs=[extract_data(x) for x in args.ref]
    targets=[extract_data(x) for x in args.target]
    results=[(i,j,np.subtract(y[2],x[2])) for i,y in enumerate(targets) for j,x in enumerate(refs)]
    
    for i,result in enumerate(results):
        print(result[0])
        name=f"{args.targetlab[result[0]]}_{args.reflab[result[1]]}"
        plot2d(targets[result[0]][0],targets[result[0]][1],results[result[0]][2],name,args.xlabel,args.ylabel,cmap,minima,min_pt)
    
    
    
        
    
    


if __name__ == "__main__":
    main()
    