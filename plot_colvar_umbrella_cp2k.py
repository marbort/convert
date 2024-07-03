import pandas
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

def read_colvar(input):
    inputs=glob.glob(f'*{input}*COLVAR*[0-9][0-9]')
    print(inputs)
    colvars=[]
    for file in inputs:
        with open(file,'r') as ifile:
            header=ifile.readline().split()[2:]
        #print(header)
        colvars.append(np.loadtxt(file,unpack=True))
    return(colvars)

def plot_colvar(data,input):
    fig=plt.figure()
    for i,colvar in enumerate(data):
        print(len(colvar[0]))
        plt.plot(range(len(colvar[0])),colvar[1],label=i)
    plt.title(input)
    plt.savefig(f'colvars_{input}.png',format='png',dpi=150)
    





input=sys.argv[1]
data=read_colvar(input)

#print(colvars)
#print(data['cvMg1IPR'])
plot_colvar(data,input)