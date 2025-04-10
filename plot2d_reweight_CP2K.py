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

"""
Takes a fes.dat input from CP2K graph tool and rescales the energy 
making the minimum a 0 and converting the values from Ha to Kj/mol
"""
ha_to_kjmol=2625.5
    

sys.path.append('/home/marco/SHARED/GitHub/mepfinder/mepfinder')
from  flooder import Flooder

def detect_local_minima(arr):
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = scipy.ndimage.generate_binary_structure(len(arr.shape),2)
    # apply the local minimum filter; all locations of minimum value 
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = (scipy.ndimage.minimum_filter(arr, footprint=neighborhood)==arr)
    # local_min is a mask that contains the peaks we are 
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    # 
    # we create the mask of the background
    background = (arr==0)
    # 
    # a little technicality: we must erode the background in order to 
    # successfully subtract it from local_min, otherwise a line will 
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = scipy.ndimage.binary_erosion(
        background, structure=neighborhood, border_value=1)
    # 
    # we obtain the final mask, containing only peaks, 
    # by removing the background from the local_min mask
    detected_minima = local_min ^ eroded_background
    return np.where(detected_minima)       

def extract_data(input,ncvs:int):
    data=[]
    cvs=[]
    
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    
    for i in range(ncvs):
        cvs.append(np.unique(file[i]))
    min=np.min(file[-1])
    val=np.subtract(file[-1],min)*ha_to_kjmol
    free_grid=val.reshape(len(cvs[0]),len(cvs[1]))
    min_pt=detect_local_minima(free_grid)
        
    return(cvs,free_grid,min_pt)
        

def get_minimum_path(cv1,cv2,free_grid):
    min_free_cv1=[(min(x),np.where(x==min(x))[0][0]) for x in free_grid]
    min_free_cv2=[(min([x[i] for x in free_grid]),np.where([x[i] for x in free_grid]==min([x[i] for x in free_grid]))) 
                   for i in range(len(cv2))]
    min_path_cv1=[(cv1[i],cv2[x[1]],x[0]) for i,x in enumerate(min_free_cv1)]
    min_path_cv2=[(cv1[x[1]],cv2[i],x[0]) for i,x in enumerate(min_free_cv2)]
    return(min_path_cv1,min_path_cv2)


def minpath(input,x_start,y_start,x_end,y_end):
    # Define your 2D free energy surface as a NumPy array
    # Replace this with your actual data
    free_energy_surface = input

    # Define the starting and ending points on the surface
    start_point = (x_start, y_start)  # Replace with your actual starting point
    end_point = (x_end, y_end)        # Replace with your actual ending point

    # Define the number of intermediate points for the string
    num_intermediate_points = 9

    # Initialize the string as a set of linearly spaced points between start and end points
    path = np.linspace(start_point, end_point, num_intermediate_points + 2)

    # Define the energy function to be minimized (distance along the string)
    def energy_function(x):
        # Calculate the sum of distances between consecutive points
        distances = np.linalg.norm(np.diff(x, axis=0), axis=1)
        return np.sum(distances)

    # Optimize the string using the minimize function from scipy
    #result = minimize(energy_function(free_energy_surface), path, method='L-BFGS-B', options={'disp': True})

    # Extract the optimized string points
    #optimized_string = result.x.reshape(-1, 2)

    # Now `optimized_string` contains the optimized path along the free energy surface
    return(optimized_string)



def dim_red(cv1,cv2,free_grid,temp):
    kbt=0.0083144621*temp #in kj/mol
    prob_grid=np.exp(-free_grid/kbt)
    reduced_prob_2=[np.trapz(y, cv2) for y in prob_grid]
    reduced_fes_2=-kbt*np.log(reduced_prob_2)
    offset2=min(reduced_fes_2)
    reduced_fes_2-=offset2
    reduced_prob_1=[np.trapz([y[i] for y in prob_grid] , cv1) for i in range(len(cv1)) ]
    reduced_fes_1=-kbt*np.log(reduced_prob_1)
    offset1=min(reduced_fes_1)
    reduced_fes_1-=offset1
    return(reduced_fes_1,reduced_fes_2)
    
    


def plot2d(x,y,value,file,labx,laby,cmap,minima):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 30}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX=100
    
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #kjmol/plot
    lev=range(0,MAX+5,5)
    plt.contourf(y, x,value,lev,vmin=0,vmax=MAX,cmap=cmap)
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
    cbar=plt.colorbar(label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(0,MAX+MAX//10,MAX//10))
    cbar.ax.set_ylim(0,MAX)
    #for i in minpath:
    #    plt.scatter(x[i[0]],y[i[1]],color='black')
    #plt.scatter(x[minima[0]],y[minima[1]])
    #plt.scatter(x[maxima[0]],y[maxima[1]],color='red')
    #plt.xlim([0.75,3.25])
    #plt.ylim([0.75,3.25])
    if minima:
        min_crd=[]
        with open(minima,'r') as ifile:
            lines=ifile.readlines()
        for line in lines:
            min_crd.append((int(line.split()[1]),float(line.split()[5]),float(line.split()[-1])))
        print(min_crd)
        for pt in min_crd:
            plt.scatter(pt[1],pt[2])
    plt.savefig('{}.png'.format(file),format='png')

def plot3d(x,y,value,input,labx,laby):
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile)
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    ax = fig.add_subplot(projection='3d')
    x=np.array([x[0] for x in file]).reshape(len(x),len(y),order='F')
    y=np.array([x[1] for x in file]).reshape(len(x),len(y),order='F')
    #x=np.unique(np.array([x[0] for x in file]))
    #y=np.unique(np.array([x[1] for x in file]))
    ax.plot_surface(x,y,value)
    ax.axes.set_zlim3d(bottom=0, top=50) 
    plt.savefig('{}_3D.png'.format(input),format='png')
    



def plotminpath(min_path_cv1,min_path_cv2,file):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    try:
        cv1=np.loadtxt("fes-rew_square_sparse_cv1.dat")
        #print(cv1)
    except:
        print("No cv1 reweight found")
    try:
        cv2=np.loadtxt("fes-rew_square_sparse_cv2.dat")
    except:
        print("No cv2 reweight found")
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    #lev=range(0,1000,5)
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #plt.plot([x[0] for x in min_path_cv1],[x[2] for x in min_path_cv1],label="CV1")
    plt.plot([x[0] for x in cv1],[x[1] for x in cv1],label="CV1 minpath")
    #plt.plot([x[1] for x in min_path_cv2],[x[2] for x in min_path_cv2],label="CV2")
    plt.plot([x[0] for x in cv2],[x[1] for x in cv2],label="CV2 minpath")
    plt.xlabel("CV")
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,2.0])
    #plt.ylim([-1.,50.0])
    plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig('{}_minpath.png'.format(file),format='png')
    
def plot_reduced(fes1,fes2,file):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    with open(fes1,'r') as ifile1:
        file1=np.loadtxt(ifile1,unpack=True)
    with open(fes2,'r') as ifile2:
        file2=np.loadtxt(ifile2,unpack=True)
    scaled1=np.subtract(file1[1],np.min(file1[1]))*ha_to_kjmol
    scaled2=np.subtract(file2[1],np.min(file2[1]))*ha_to_kjmol
    plt.plot(file1[0],scaled1,label="CV1")
    plt.plot(file2[0],scaled2,label="CV2")
    plt.xlabel("CV")
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,4.0])
    #plt.ylim([-1.,50.0])
    plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig('{}_reduced.png'.format(file),format='png')
    


def main():

    cvs,free_grid,min_pt=extract_data(sys.argv[1],2)
    file=os.path.splitext(sys.argv[1])[0]
    labx=sys.argv[2]
    laby=sys.argv[3]
    minima=sys.argv[4]
    print(minima)
    cmap_active='jet'
    #cmap_active=ListedColormap(np.linspace([0.16862745098, 0.219607843137,1,1],[1, 1,1,1],12)) #blue
    #cmap_active=ListedColormap(np.linspace([1.0, 0.40784313725490196,0,1],[1, 1,1,1],12)) #orange
    #cmap_active=ListedColormap(np.linspace([0.1960784313725490, 0.7686274509803922,0.4980392156862745,1],[1, 1,1,1],12)) #green
    0.1960784313725490
    
    #0.16862745098, 0.219607843137,1
    #cmap_active=LinearSegmentedColormap.from_list("mycmap",["#2b38ff","#FFFFFF"])#"#24BC99","#D6BB61","#E18F26","#FF6800"])
                                                  #oldcmap"#000000","#2b38ff","#17d9ff","#f7059b"])
                                                  
    #cmap_active='Blues_r'
    #path=minpath(free_grid,1,1,2,2)
    #min_path_cv1,min_path_cv2=get_minimum_path(cv1,cv2,free_grid)
    #reduced_fes_1,reduced_fes_2=dim_red(cv1,cv2,free_grid,300)
    #minima=detect_local_minima(free_grid)
    #maxima=detect_local_minima(-free_grid)
    plot2d(cvs[0],cvs[1],free_grid,file,labx,laby,cmap_active,minima)
    #plot3d(cv1,cv2,free_grid,file+'.dat',labx,laby)
    #plotminpath(min_path_cv1,min_path_cv2,file)
    plot_reduced('fes1.dat','fes2.dat',file)

if __name__ == "__main__":
    main()


