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



def plot3d(x,y,value,input,labx,laby):
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile)
    fig=plt.figure(figsize=(20,20),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    ax = fig.add_subplot(projection='3d')
    x=np.array([x[0] for x in file]).reshape(len(x),len(y),order='F')
    y=np.array([x[1] for x in file]).reshape(len(x),len(y),order='F')
    #x=np.unique(np.array([x[0] for x in file]))
    #y=np.unique(np.array([x[1] for x in file]))
    #ax.plot_surface(x,y,value,cmap='rainbow')
    ax.set(xlim=(-.1, 3.1), ylim=(-.1, 3.1), zlim=(-40, 120))
    ax.plot_surface(x, y, value, edgecolor='royalblue', lw=0.5, 
                alpha=0.7,cmap='rainbow',vmin=0,vmax=MAX)
    
    ax.contourf(x, y, value, zdir='z', offset=-20, cmap='rainbow',vmin=0,vmax=MAX)
    ax.contour(x, y,value,levels=range(0,int(MAX),20),offset=-20,vmin=0,vmax=MAX,linewidths=1.5,colors='black')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_axis_off()
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.tight_layout()
    plt.savefig('{}_3D.png'.format(input),format='png')

def extract_data(input):
    data=[]
    
    
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile)
    cv1_temp=np.array([x[0] for x in file])
    cv2_temp=np.array([x[1] for x in file])
    cv1=np.unique(cv1_temp)
    cv2=np.unique(cv2_temp)
    val=np.array([x[2] for x in file])
    free_grid=val.reshape(len(cv1),len(cv2))
    free_grid_limit=free_grid.copy()
    free_grid_limit[free_grid_limit > 120] = 120
    return(cv1,cv2,free_grid_limit)
    

cv1,cv2,free_grid=extract_data(sys.argv[1])
#cv1_state,cv2_state,free_grid_state,min_pt_state=extract_data(sys.argv[2])
#print(cv1[min_pt[1][0]],cv2[min_pt[0][0]],free_grid[min_pt[0][0],min_pt[1][0]])
file=os.path.splitext(sys.argv[1])[0]
labx=sys.argv[2]
laby=sys.argv[3]
MAX=sys.argv[4]
plot3d(cv1,cv2,free_grid,file+'.dat',labx,laby)