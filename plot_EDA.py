import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def extract_data(path):
    data=pd.read_excel(path,sheet_name="EDA", skiprows=47,usecols=[1,2,3])
    return(data)

def plot_EDA(data,geoms):
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
    symbols=["o","s","v"]
    found=0
    for row in data.iterrows():
        #print(row[1][0])
        if row[1][0] == "PAULI":
            styl="dotted"
        if row[1][0] == "ELSTAT":
            styl="dashed"
        if row[1][0] == "OI":
            styl="dashdot"
        if row[1][0] == "TOTAL":
            styl="solid"
        if row[1][0] in geoms:
            symbol=symbols[found%len(geoms)]
            found+=1
            plt.plot(range(len(row[1])-1),[x for x in row[1][1:]],marker='s',linestyle=styl,label=geoms[geoms.index(row[1][0])])
    plt.legend()
    plt.savefig('test_eda.png',format='png')
        
        
geoms=["Me_geo2","Et_geo2","iPr_geo2"]
data=extract_data('/home/marco/OneDrive_UiO/2024/ML_Schlenk/EDA_summary.ods')
plot_EDA(data,geoms)

