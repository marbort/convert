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
    parser.add_argument('--output', dest='output', 
                        type=str, help='output file',default="2D_reduced_compare.png")
    parser.add_argument('--title', dest='title', 
                        type=str, help='plot title',default="")
    parser.add_argument('--labels', dest='labels', default=None, type=str, nargs='+',
                         help='Plot labels label')
    parser.add_argument('--xlab', dest='xlab', type=str, 
                         help='X axis label',default='CV1')
    parser.add_argument('--min', dest='min', 
                        type=float, help='CV value for 0')
    parser.add_argument('--limy', dest='limy', 
                        type=float, help='y MAX value')
    parser.add_argument('--errors', dest='errors', 
                        help='plot errors',action='store_true',default=False)
    parser.add_argument('--exterrors', dest='exterrors', 
                        help='take errors from external file',default=None,type=str)
    parser.add_argument('--range', dest='range', 
                        help='set X range',nargs='+',type=float)
    parser.add_argument('--Emax', dest='Emax', 
                        help='set max of E',type=float)
    
    args = parser.parse_args()
    #pres_colors=["#2b38ff","#f7059b","#17d9ff","#000000","#4cb944"]
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    cmap = mpl.colormaps['rainbow']
    cmap_colors=[cmap(0),cmap(0.5),cmap(0.75),cmap(1.0)]
    data={}
    err_cv1=[]
    lw=4
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    for i,input in enumerate(args.input):
        cv_min=0
        try:
            print(input)
            cv1,min_path_cv1=np.loadtxt(input,unpack=True)
        except:
            print(input)
            cv1,min_path_cv1,err_cv1=np.loadtxt(input,unpack=True)
        
        min_path_cv1[min_path_cv1>args.Emax]=args.Emax
        
        if args.min:
             cv_min=find_nearest(cv1,args.min)
             plt.plot(cv1,[x-min_path_cv1[np.where(cv1==cv_min)] for x in min_path_cv1],label=args.labels[i],linewidth=lw,color=pres_colors[i])
        else:
            plt.plot(cv1,min_path_cv1,label=args.labels[i],color=pres_colors[i],linewidth=lw)
        
        if args.errors:
            if args.exterrors:
                try:
                    exterrs=np.loadtxt(args.exterrors,unpack=True)
                    err=[x for x in  exterrs[1]]
                except:
                    print("Error reading external errors file.")
                    sys.exit()
            else:
                err=[x for x in err_cv1]
            if args.min:
                plt.fill_between(cv1,[float(x-min_path_cv1[np.where(cv1==cv_min)]+err[i]) for i,x in enumerate(min_path_cv1)],
                                 [float(x-min_path_cv1[np.where(cv1==cv_min)]-err[i]) for i,x in enumerate(min_path_cv1)],linewidth=0,alpha=0.1,color=pres_colors[i])
            else:
                plt.fill_between(cv1,min_path_cv1+err,min_path_cv1-err,linewidth=0,alpha=0.1,color=pres_colors[i])
            
                
    plt.xlim(args.range[0],args.range[1])  
    plt.xlabel(args.xlab)
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.8,2.2])
    plt.ylim([-5,args.limy])
    plt.title(args.title,y=1.05)
    plt.legend(loc='center left',bbox_to_anchor=(1, 0.5),frameon=False)
    plt.tight_layout()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig(args.output,format='png')
    #print(len(data['free']),len(data['err']))
    
if __name__=="__main__":
    main()