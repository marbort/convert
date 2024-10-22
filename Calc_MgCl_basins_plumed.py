import dpdata as dp
import glob
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import regex as re
from natsort import natsorted, ns
import json
import argparse


def get_data(inputs,symm):
    percent={}
    histo={}
    paths=sorted(glob.glob(inputs))
    symmetrized={}
    symm=[(x.split(',')[0],x.split(',')[1]) for x in symm]+[(x.split(',')[1],x.split(',')[0]) for x in symm]
    for input in paths:
        point=(input.split('_')[-2],input.split('_')[-1])
        with open(input, 'r') as ifile:
            data=np.loadtxt(ifile,unpack=True)
        histo[point]=[x for x in data[1]]
        hist_perc=[x/sum(data[1]) for x in data[1]]
        percent[point]=[data[0],hist_perc]
    for point in percent:
        if point in symm:
            symmetrized[point]=[(x+histo[point[1],point[0]][i])/(sum(histo[point]+histo[point[1],point[0]])) for i,x in enumerate(histo[point]) ]
            #(np.average(percent[(point[0],point[1])][1]+percent[(point[1],point[0])][1])/2
    for point in symmetrized:
        percent[point][1]=symmetrized[point]
    
    return(percent)

def plot_histo(histo,maxcv3,width,symm_minima,histo_symm_sum):
    fig=plt.figure(figsize=(16,16),dpi=150)
    
    for i,hist in enumerate(histo):
        plt.subplot(3,3,i+1)
        plt.xlim(0,maxcv3+width)
        if symm_minima:
            if i in symm_minima:
                plt.bar([x+width/2 for x in histo[hist][1][:-1]],histo_symm_sum,edgecolor='black',width=width)
            else:
                plt.bar([x+width/2 for x in histo[hist][1][:-1]],histo[hist][0],edgecolor='black',width=width)
        else:
            plt.bar([x+width/2 for x in histo[hist][1][:-1]],histo[hist][0],edgecolor='black',width=width)
        plt.title(hist)
    plt.savefig("MgCl_hist.png",format='png')


        

def write_percentage(cvs,percent):
    with open('basin_analysis.txt','w') as ofile:
        ofile.write(f"CVS:{' '.join(cvs)}\n")
        for k,point in enumerate(percent):
            """
            if symm_minima:
                if k in symm_minima:
                    ofile.write(f"Composition of minimum at {point} after symmetryzation \n")
                    for i,entry in enumerate(percent_symm_sum):
                        ofile.write(f"CN_MgCl {percent[point][1][i]:.1f}: {entry:.3f} \n")
                else:
                    ofile.write(f"Composition of minimum at {point}\n")
                    for i,entry in enumerate(percent[point]):
                        ofile.write(f"CN_MgCl {percent[point][1][i]:.1f}: {entry:.3f} \n")
            else:
            """
            ofile.write(f"Composition of minimum at {point}\n")
            for i,entry in enumerate(percent[point][0]):
                ofile.write(f"CN_MgCl {entry:.1f}: {percent[point][1][i]:.3f} \n")
    
                    


def main():
    
    parser = argparse.ArgumentParser(description='Calculate percentage of CV3 in CV1 CV2 basins')

    parser.add_argument('--inputs',dest='inputs',type=str,default='COLVAR',help='plumed histogram files')
    parser.add_argument('--cvs',dest='cvs',type=str,default='cv1',help='name of cvs to analyze',nargs='+')
    parser.add_argument('--cv1min',dest='mincv1',type=float,help='minima of cv1',nargs='+')
    parser.add_argument('--cv2min',dest='mincv2',type=float,help='minima of cv2',nargs='+')
    parser.add_argument('--cv3max',dest='maxcv3',type=float,help='max of cv3')
    parser.add_argument('--symm',dest='symm',type=str,help='Symmetric minima',nargs="+")
    parser.add_argument('--temp',dest='temp',type=float,help='Temperature for reweighting')
    parser.add_argument('--width',dest='width',type=float,help='width of the hist bins')
    
    
    args = parser.parse_args()
    
    
    
    histo={}
    percent=get_data(args.inputs,args.symm)
    print(percent)
    #plot_histo(histo,args.maxcv3,args.width,args.symm,histo_symm_sum)
    #print(points,percent)
    write_percentage(args.cvs,percent)


if __name__ == "__main__":
    main()
        
        
        