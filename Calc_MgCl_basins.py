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



def calc_percent_basins(colvar,cvs,mincv1,mincv2,tol,maxcv3,width,symm_minima,temp):
    kB=8.314462618E-3 
    basins=[tuple([x,mincv2[i]]) for i,x in enumerate(mincv1)]
    points={x:[] for x in basins}
    histo={}
    percent={}
    cv_split={x:[] for x in basins}
    biases={x:[] for x in basins}
    histo_symm_sum=None
    percent_symm_sum=None

    with open(colvar,'r') as ifile:
        lines=ifile.readlines()
    header=lines[0].split()[2:]
    bias_idx=[i for i,item in enumerate(header) if ".bias" in item]
    cv1=header.index(cvs[0])
    cv2=header.index(cvs[1])
    cv3=header.index(cvs[2])
    for j,line in enumerate(lines[1000:]):
        split=[float(x) for x in line.split()]
        for k,basin in enumerate(basins):
            if basin[0]-tol[k] < split[cv1] < basin[0]+tol[k] and basin[1]-tol[k] < split[cv2] < \
                basin[1]+tol[k]:
                bias=[split[x] for x in bias_idx]
                bias_total=sum(bias)
                points[basin].append((split[cv3],split[cv1],split[cv2],j))
                if temp == 0:
                    biases[basin].append(1)
                else:       
                    #weights[basin].append(-kB*temp*np.log(np.exp(bias_total/(kB*temp))))       
                    biases[basin].append(bias_total)
                cv_split[basin].append(line)
    for item in cv_split:
        item_name="_".join([str(x) for x in item])
        with open(f'colvar_{item_name}.dat','w') as ofile:
            ofile.write(lines[0])
            for line in cv_split[item]:
                ofile.write(line)
                    

    
    
    for  point in points:
        #print(points)
        structs=[x[0] for x in points[point]]
        #print(structs)
        bins=int(maxcv3/width)
        biases_shift=[x-max(biases[point]) for x in biases[point]]
        weights=[np.logaddexp.reduce(x) for x in biases_shift]
        #histo[point]=np.histogram(structs,bins=bins,range=(0,maxcv3),density=True)
        histo[point]=np.histogram(structs,bins=bins,range=(0,maxcv3),weights=weights,density=True)
        percent[point]=[x/sum(histo[point][0]) for x in histo[point][0]]
    
    if symm_minima:
        histo_symm=[histo[list(histo.keys())[i]][0] for i in symm_minima]
        histo_symm_sum=[sum(x) for x in zip(*histo_symm)]
        percent_symm_sum=[x/sum(histo_symm_sum) for x in histo_symm_sum]
        #print(histo,histo_symm_sum)
        return(points,histo,percent,histo_symm_sum,percent_symm_sum)
    
    else:
        #print(histo,histo_symm_sum)
        return(points,histo,percent,histo_symm_sum,percent_symm_sum)

def plot_histo(histo,maxcv3,width,symm_minima,histo_symm_sum):
    fig=plt.figure(figsize=(16,16),dpi=150)
    print(histo)
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


        

def write_percentage(cvs,histo,percent,symm_minima,percent_symm_sum):
    with open('basin_analysis.txt','w') as ofile:
        ofile.write(f"CVS:{' '.join(cvs)}\n")
        for k,point in enumerate(percent):
            if symm_minima:
                if k in symm_minima:
                    ofile.write(f"Composition of minimum at {point} after symmetryzation \n")
                    for i,entry in enumerate(percent_symm_sum):
                        ofile.write(f"CN_MgCl {histo[point][1][i]:.1f}: {entry:.3f} \n")
                else:
                    ofile.write(f"Composition of minimum at {point}\n")
                    for i,entry in enumerate(percent[point]):
                        ofile.write(f"CN_MgCl {histo[point][1][i]:.1f}: {entry:.3f} \n")
            else:
                ofile.write(f"Composition of minimum at {point}\n")
                for i,entry in enumerate(percent[point]):
                    ofile.write(f"CN_MgCl {histo[point][1][i]:.1f}: {entry:.3f} \n")
        
                    


def main():
    
    parser = argparse.ArgumentParser(description='Calculate percentage of CV3 in CV1 CV2 basins')

    parser.add_argument('--input',dest='input',type=str,default='COLVAR',help='full COLVAR file')
    parser.add_argument('--cvs',dest='cvs',type=str,default='cv1',help='name of cvs to analyze',nargs='+')
    parser.add_argument('--cv1min',dest='mincv1',type=float,help='minima of cv1',nargs='+')
    parser.add_argument('--cv2min',dest='mincv2',type=float,help='minima of cv2',nargs='+')
    parser.add_argument('--cv3max',dest='maxcv3',type=float,help='max of cv3')
    parser.add_argument('--tol',dest='tol',type=float,help='tolerance for basins',nargs='+')
    parser.add_argument('--width',dest='width',type=float,help='width of the hist bins')
    parser.add_argument('--symm',dest='symm',type=int,help='Symmetric minima',nargs="+")
    parser.add_argument('--temp',dest='temp',type=float,help='Temperature for reweighting')
    
    
    
    args = parser.parse_args()
    
    
    points,histo,percent,histo_symm_sum,percent_symm_sum=calc_percent_basins(args.input,args.cvs,args.mincv1,\
                                                                            args.mincv2,args.tol,args.maxcv3,args.width,args.symm,args.temp)
    #print(histo_symm_sum)
    #plot_histo(histo,args.maxcv3,args.width,args.symm,histo_symm_sum)
    #print(points,percent)
    #write_percentage(args.cvs,histo,percent,args.symm,percent_symm_sum)


if __name__ == "__main__":
    main()
        
        
        