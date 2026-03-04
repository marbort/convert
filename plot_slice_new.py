import os
import glob
import sys
import numpy as np
import argparse
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import matplotlib as mpl



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
parser.add_argument('--minima', type=float, nargs='+',
                    help='Plot other cv slices for minima of specified cv value',default=False)
parser.add_argument('--title', type=str, default="Free energy slices",
                    help='Title for the plot')
parser.add_argument('--xlabel', type=str, default="CV1",
                    help='Label for the x-axis')
parser.add_argument('--ylabel', type=str, default="Free energy (kJ/mol)",
                    help='Label for the y-axis')
parser.add_argument('--cv2label', type=str, default="CV2",
                    help='Label for the second CV')
parser.add_argument('--kcal', action='store_true',
                    help='Convert free energy to kcal/mol')

    
args = parser.parse_args()


data=np.loadtxt(args.input,unpack=True)

cv1=np.unique(data[0])
cv2=np.unique(data[1])
slices={"cv1":[],"cv2":[],"minima":[]}
zero=cv1[get_closest(cv1,0.0)]
zero_idx=get_closest(cv1,0.0)
print(f"Zero index: {zero_idx} for cv1 value {zero}")   
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
    graphs['cv1'][cv1[slice]]=slice_tmp

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


if args.minima:
    zero=cv1[get_closest(cv1,0.0)]
    slices['minima']=argrelextrema(graphs['cv1'][zero][1], np.less)[0]
    selected=[]
    for i in args.minima:
        selected.append(get_closest(cv2,i))
    print(f"Selected minima: {[cv2[x] for x in selected]}")
        
    
    graphs={"cv1":{},"cv2":{},'minima':{}}
    for slice in selected:
    #for slice in slices['minima'][:3]:
        print(selected)
        slice_tmp=[[],[]]
        for i,item in enumerate(data[1]):
            if cv2[slice] == item:
                slice_tmp[0].append(data[0][i])
                slice_tmp[1].append(data[2][i])
        slice_tmp[1]=slice_tmp[1]-min(slice_tmp[1])
        graphs['minima'][cv2[slice]]=slice_tmp
    for cv in graphs:
        for slice in graphs[cv]:
            with open(f'{args.output}_{cv}_{slice}.dat','w') as ofile:
                for i,line in enumerate(graphs[cv][slice][0]):
                    ofile.write(f"{line} {graphs[cv][slice][1][i]}\n")
    try:
        minpath=np.loadtxt(f"path_clean.dat", unpack=True)
        has_minpath=True
    except:
        print("No path.dat found, skipping minimum path plot.")
        has_minpath=False
    
    fig = plt.figure(figsize=(16, 12), dpi=150)
    font = {"family": "Formular", "weight": "normal", "size": 46}
    mpl.rc("font", **font)
    mpl.rcParams["axes.linewidth"] = 3
    mpl.rcParams["lines.linewidth"] = 3
    
    for slice in selected:
        y= graphs['minima'][cv2[slice]][1]- graphs['minima'][cv2[slice]][1][zero_idx]
        if args.kcal:
            y = y / 4.184
        
        plt.plot(cv1, y, label=f"{args.cv2label} {round(cv2[slice],0)}")
    if has_minpath:
        zero_idx = get_closest(minpath[0], zero)
        y= minpath[2] - minpath[2][zero_idx]
        if args.kcal:
            y = y / 4.184
        plt.plot(minpath[0], y, 'r--', label="Minimum path")
    plt.xlabel(args.xlabel)
    if args.kcal:
        plt.ylabel("Free energy (kcal/mol)")
    else:
        plt.ylabel("Free energy (kJ/mol)")
    plt.title(args.title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{args.output}_minima.png", dpi=150)
    plt.close(fig)
    
    
    
    
        
    
            
        

    
    
    



