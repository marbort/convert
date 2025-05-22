import pandas
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob

def read_colvar(input):
    with open(input,'r') as ifile:
        header=ifile.readline().split()[2:]
    #print(header)
    colvar=pandas.read_csv(input,sep=" ",header=None,names=header,skiprows=1,comment="#")
    return(colvar)

def plot_colvar(colvars,data,input):
    fig=plt.figure()
    for colvar in colvars:
        plt.plot(range(len(data[colvar])),data[colvar])
    plt.legend(colvars,loc='upper right')
    plt.savefig(f'{input}_{"_".join(colvars)}.png',format='png',dpi=150)

def plot_colvar_multi(colvars,data):
    fig=plt.figure(figsize=(16,9))
    for i,replica in enumerate(data):
        plt.subplot(2,len(data)//2+len(data)%2,i+1)
        for colvar in colvars:
            plt.plot(range(len(replica[colvar])),replica[colvar])
        
        plt.legend(loc='upper right')
    
    plt.savefig(f'colvars_all_{"_".join(colvars)}.png',format='png',dpi=150)
    
def histo_colvar(colvars,data,input):
   
    fig=plt.figure(figsize=(16,9))
    for i,colvar in enumerate(colvars):
        plt.subplot(1,len(colvars),i+1)
        plt.hist(data[colvar],bins=len(colvars)*10,label=colvars[i])
        print(colvars[i])
        plt.legend(loc='upper right')
    plt.savefig(f'{input}_HIST_{"_".join(colvars)}.png',format='png',dpi=150)


    

def get_colvars_meta(plumed,colvars):
    
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

parser = argparse.ArgumentParser(description='Extract colvars from plumed.dat')
parser.add_argument('--colvars', type=str, help='Colvars to extract', nargs='+', default=None)
parser.add_argument('--input', type=str, help='Input file', default='colvar')
parser.add_argument('--multi', action='store_true', help='Plot multiple colvar files in same plot', default=False)
    
args=parser.parse_args()
    



colvars=get_colvars_meta('plumed.dat',args.colvars)
#print(colvars)
#print(data['cvMg1IPR'])
if args.multi:
    data=[]
    inputs=glob.glob(f'{args.input}.[0-9]')
    for input in inputs:
        data.append(read_colvar(input))
    plot_colvar_multi(colvars,data)
    
    
else:
    data=read_colvar(args.input)   
    plot_colvar(colvars,data,args.input)
    histo_colvar(colvars,data,args.input)