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


#Only implemented for 2 CVs. CVs must be specified in input in the order x y and not otherwise
def find_values(input,CV,valfile,tol):
    cv_ok_full=[]
    cv_both_full=[]
    with open(input,'r') as ifile:
        first_line=ifile.readline()
    header=first_line.split()
    idx = [i-2 for x in CV for i in range(len(header)) if header[i] == x ]
    results=np.loadtxt(input,unpack=False)
    print(f"Number of points: {len(results)}")
    with open(valfile,'r') as val_input:
        val_lines=val_input.readlines()
    for line in val_lines:
        val=[float(x) for x in line.split()]
        cv_ok=[[j for j,x in enumerate(results) if val[i]-tol[i] < x[k] < val[i]+tol[i]] for i,k in enumerate(idx) ]
        cv_ok_full.append((val[0],val[1],cv_ok))
        cv_both=set(cv_ok[0])&set(cv_ok[1])
        diff_sum=[abs(results[i][idx[0]]-val[0])+abs(results[i][idx[1]]-val[1]) for i in cv_both]
        cv_both=list(zip(cv_both,diff_sum))
        cv_both.sort(key=lambda a: a[1])
        cv_both_full.append((val[0],val[1],cv_both))
    return(results,idx,cv_ok_full,cv_both_full)



parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--input' , dest='input',help='lammps trajectory')
parser.add_argument('--CVs', dest='CVs', default=['CV1','CV2'],
                     help='CVs to check',nargs='+')
parser.add_argument('--val', dest='val', 
                    type=str,help='File with values of the CVs')
parser.add_argument('--tol', dest='tol', default=[0.1,0.1],
                    type=float,help='Tolerance on the values of the CVs',nargs='+')




args = parser.parse_args()




results,idx,cv_ok_full,cv_both_full=find_values(args.input,args.CVs,args.val,args.tol)
with open('structure_minima.dat','w') as ofile:
    with open('structures.dat','w') as ofile2:
        for i,item in enumerate(cv_ok_full):
            print(f"#### CV1={item[0]}  CV2={item[1]}")
            print(f"Number of points that satisfy {args.CVs[0]}: {len(item[2][0])}")
            print(f"Number of points that satisfy {args.CVs[1]}: {len(item[2][1])}")
            print(f"Structures that satisfy both criteria: {len(cv_both_full[i][2])}")
            if len(cv_both_full[i][2]) > 0:
                for j,k in enumerate(cv_both_full[i][2]):
                    if j==0:
                        ofile.write(f"Structure {k[0]} with {args.CVs[0]} value {results[k[0]][idx[0]]:.4f} "
                            f"and {args.CVs[1]} value {results[k[0]][idx[1]]:.4f}\n")
                        ofile2.write(f"{k[0]}\n")
                    print(f"Structure {k[0]} with {args.CVs[0]} value {results[k[0]][idx[0]]:.4f} "
                        f"and {args.CVs[1]} value {results[k[0]][idx[1]]:.4f}")

    

