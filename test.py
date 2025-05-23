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
        file=np.loadtxt(ifile)
    cv1_temp=np.array([x[0] for x in file])
    cv2_temp=np.array([x[1] for x in file])
    cv1=np.unique(cv1_temp)
    cv2=np.unique(cv2_temp)
    val=np.array([x[2] for x in file])
    free_grid=val.reshape(len(cv1),len(cv2))
    #min_pt=detect_local_minima(free_grid)1
    return(cv1,cv2,free_grid)


def symmetryze(free_grid,cv1,cv2):
    symm_free_grid=np.empty_like(free_grid)
    vals=(2.5,2.0)
    print(f"Idxs {vals}: {np.where(cv2==vals[1])[0][0]}")
    #print(free_grid[60][50])
    print(free_grid[np.where(cv2==vals[1])[0][0]][np.where(cv1==vals[0])[0][0]])
    print(free_grid[np.where(cv2==vals[0])[0][0]][np.where(cv1==vals[1])[0][0]])
    target=(free_grid[50][60]+free_grid[60][50])/2
    print(f"Target={target}")
    #print(f"Target: {(free_grid[f][k]+free_grid[-1-k][-1-f])/2}")
    for i in range(len(free_grid)):
        for j in range(len(free_grid[i])):
            symm_free_grid[i][j]=(free_grid[i][j]+free_grid[j][i])/2
    print(symm_free_grid[60][50])
    print(symm_free_grid[50][60])
    return(symm_free_grid)
            
        
        
    
 
 
    

cv1,cv2,free_grid=extract_data('fes-rew_square_sparse.dat')
#print(len(cv1),len(cv2),len(free_grid))
symmetryze(free_grid,cv1,cv2)
