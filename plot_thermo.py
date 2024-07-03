import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import argparse
import sys
import plumed
import pandas as pd

def plot_data(input):
    fields=["temp", "etotal", "pe", "ke"]
    data=pd.read_table(input,sep=" ",index_col=False)
    data.columns=list(data.columns[1:])+ [" "]
    
    fig=plt.figure(figsize=(16,10),dpi=150)
    for i,field in enumerate(fields):
        plt.subplot(2,2,i+1)
        plt.plot(data['time'],data[field])
        plt.xlabel("Time (ps)")
        plt.ylabel(field)
    plt.suptitle(os.path.abspath(input))
    plt.savefig('thermo.png',format='png')

    


def main():
    plot_data('thermo.out')



if __name__=="__main__":
    main()
    
    
    
    