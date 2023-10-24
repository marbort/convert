import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import os


parser = argparse.ArgumentParser(description='Plot data')
parser.add_argument('--input' , dest='input',help='string to match input files')
args = parser.parse_args()


with open(args.input,'r') as ifile:
    data=np.loadtxt(ifile)



fig=plt.figure(figsize=(15,10),dpi=150)
plt.errorbar([x[0] for x in data],[x[1] for x in data],yerr=[x[2]*1.96 for x in data],
                        ecolor="#D3D3D3", marker='o',markersize='3') #ecolor="#D3D3D3"
plt.savefig('umbrella_fes+error.png',format='png')

