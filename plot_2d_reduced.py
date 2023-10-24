import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np


def plot_reduced(cv1,cv2,min_path_cv1,min_path_cv2,file):
    
    plt.plot(cv1,min_path_cv1,label="CV1")
    plt.plot(cv2,min_path_cv2,label="CV2")
  
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
    
def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data',nargs="+")
    parser.add_argument('--labels', dest='labels', default=None, type=str, nargs='+',
                         help='Plot labels label')
    parser.add_argument('--min', dest='min', 
                        type=float, help='CV value for 0')
    
    args = parser.parse_args()
    pres_colors=["#2b38ff","#f7059b","#17d9ff","#000000","#4cb944"]
    data={}
    err_cv1=[]
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 22}
    mpl.rc('font', **font)
    for i,input in enumerate(args.input):
        try:
            cv1,min_path_cv1=np.loadtxt(input,unpack=True)
        except:
            cv1,min_path_cv1,err_cv1=np.loadtxt(input,unpack=True)
        if args.min:
             cv_min=find_nearest(cv1,args.min)
             plt.plot(cv1,[x-min_path_cv1[np.where(cv1==cv_min)] for x in min_path_cv1],label=args.labels[i],linewidth=3,color=pres_colors[i])
        else:
            plt.plot(cv1,min_path_cv1,label=args.labels[i],color=pres_colors[i])
        
        if err_cv1:
            err=[x*2 for x in err_cv1]
            plt.fill_between(cv1,min_path_cv1+err,min_path_cv1-err,linewidth=0,alpha=0.1,color=pres_colors[i])
            
    plt.xlabel("CV1")
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,3.0])
    #plt.ylim([-5.,40.0])
    plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig('2D_reduced_compare.png',format='png')
    #print(len(data['free']),len(data['err']))
    
if __name__=="__main__":
    main()