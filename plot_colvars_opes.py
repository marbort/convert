import pandas
import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_colvar(input):
    with open(input,'r') as ifile:
        header=ifile.readline().split()[2:]
    #print(header)
    colvar=pandas.read_csv(input,sep=" ",header=None,names=header,skiprows=1)
    return(colvar)

def plot_colvar(colvars,data):
    fig=plt.figure()
    for colvar in colvars:
        plt.plot(range(len(data[colvar])),data[colvar])
    plt.legend(colvars)
    plt.savefig(f'colvars_{"_".join(colvars)}.png',format='png',dpi=150)
    
def histo_colvar(colvars,data):
   
    fig=plt.figure(figsize=(16,9))
    for i,colvar in enumerate(colvars):
        plt.subplot(1,len(colvars),i+1)
        plt.hist(data[colvar],bins=len(colvars)*10,label=colvars[i])
        print(colvars[i])
        plt.legend()
    plt.savefig(f'colvars_HIST_{"_".join(colvars)}.png',format='png',dpi=150)


    

def get_colvars_meta(plumed):
    parser = argparse.ArgumentParser(description='Extract colvars from plumed.dat')
    parser.add_argument('--colvars', type=str, help='Colvars to extract', nargs='+', default=None)
    
    args=parser.parse_args()
    
    colvars=[]
    if args.colvars:
        colvars=args.colvars
    else:
        with open(plumed, 'r') as ifile:
            lines=ifile.readlines()
        for line in lines:
            if "opes:" in line or "opes1:" in line  or "opes2:" in line:
                colvars+=([x.split('=')[-1].split(',') for x in line.split() if "ARG" in x][0])
    return(colvars)

data=read_colvar('colvar')
colvars=get_colvars_meta('plumed.dat')
#print(colvars)
#print(data['cvMg1IPR'])
plot_colvar(colvars,data)
histo_colvar(colvars,data)