import os
import glob
import sys
import numpy as np
import argparse

def get_closest(data,max):
    diffs=[abs(x-float(max)) for x in data]
    limit=diffs.index(min(diffs))
    return(limit)


parser = argparse.ArgumentParser(description='Plot data')
parser.add_argument('--input', dest='input', 
                    type=str, help='input data')
parser.add_argument('--output', dest='output', 
                    type=str, help='output file',default="fes_slice")
parser.add_argument('--cv1', dest='cv1', 
                    type=str, help='slices values for CV1',default=[],nargs='+')
parser.add_argument('--cv2', dest='cv2', 
                    type=str, help='slices values for CV2',default=[],nargs='+')

    
args = parser.parse_args()


data=np.loadtxt(args.input,unpack=True)

cv1=np.unique(data[0])
cv2=np.unique(data[1])
slices={"cv1":[],"cv2":[]}
for i in args.cv1:
    slices['cv1'].append(get_closest(cv1,i))
for i in args.cv2:
    slices['cv2'].append(get_closest(cv2,i))

graphs={"cv1":{},"cv2":{}}
for slice in slices['cv1']:
    slice_tmp=[[],[]]
    for i,item in enumerate(data[0]):
           if cv1[slice] == item:
            slice_tmp[0].append(data[1][i])
            slice_tmp[1].append(data[2][i])
    slice_tmp[1]=slice_tmp[1]-min(slice_tmp[1])
    graphs['cv1'][slice]=slice_tmp

for slice in slices['cv2']:
    slice_tmp=[[],[]]
    for i,item in enumerate(data[1]):
           if cv2[slice] == item:
            slice_tmp[0].append(data[0][i])
            slice_tmp[1].append(data[2][i])
    slice_tmp[1]=slice_tmp[1]-min(slice_tmp[1])
    graphs['cv2'][cv2[slice]]=slice_tmp
for cv in graphs:
    for slice in graphs[cv]:
        with open(f'{args.output}_{cv}_{slice}.dat','w') as ofile:
            for i,line in enumerate(graphs[cv][slice][0]):
                ofile.write(f"{line} {graphs[cv][slice][1][i]}\n")
            
        

    
    
    



