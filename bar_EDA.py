import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def extract_paths(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    paths=[x.rstrip() for x in lines]
    return(paths)


def bar_EDA(paths,labx,laby):
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    structs=[]
    labels=["Me","Et","iPr"]
    tix=["Pauli","Elstat","OI","Total1"]
    fig=plt.figure(figsize=(32,32),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 60}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    for path in paths:
        sheet=pd.read_excel(path,sheet_name="EDA")
        structs.append(sheet)
    
    width=0.2
    multiplier=0
    idx=[7,9,13,26]
    total=[]
    for j,mol in enumerate(structs):
            offset=width*multiplier
            plt.bar([x+offset for x in range(len(idx))],[mol["kcal/mol"][k] for k in idx],width=width,color=pres_colors[j],label=labels[j])
            multiplier+=1.1
            total.append(mol["kcal/mol"][idx[-1]])
    plt.xticks(range(len(idx)),labels=tix)
    plt.legend()
    plt.tight_layout()
    plt.savefig("EDA.png")
    total=[x-min(total) for x in total]
    print(total)

paths=extract_paths('/run/media/marco/T7 Shield/SHARED/RATIO/WP1/ML/fox/EDA_paths.txt')
bar_EDA(paths,"A","B")
