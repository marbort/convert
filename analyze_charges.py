import numpy as np
import glob
import argparse
import matplotlib
import matplotlib.pyplot as plt
import os
from operator import add
import json



def load_data(input,Vor):
    with open(input,'r') as ifile:
        data=json.load(ifile)
    with open('order.json','r') as ifile:
        order=json.load(ifile)
    if Vor:
        with open(Vor,'r') as ifile:
            VDD=json.load(ifile)
        return(data,VDD,order)
    else:
        return(data,order)

def mean_charge(data,frag):
    
    for i in data:
        frg=[]
        for j in data[i][frag]:
            tmp=[float(x.split()[-1]) for x in j[1:] ]
            frg.append(tmp)
        data[i][frag].append(np.mean(frg,axis=0))
        data[i][frag].append(np.std(frg,axis=0))
    #print(data['2.7']['combo'][-2])
    return(data)

def mean_VDD(data,frag):
    frg=[]
    for i in data:
        for j in data[i][frag]:
            tmp=[float(x.split()[-1]) for x in j ]
            frg.append(tmp)
        data[i][frag].append(np.mean(frg,axis=0))
        data[i][frag].append(np.std(frg,axis=0))
        
    return(data)

def mean_charge_O(data,frag,order):
    O1_mean=[]
    O2_mean=[]
    wind=[]
    Ochrg={}
    for i in data:
        wind.append(i)
        O1=[]
        O2=[]
        for j in data[i][frag]:
            O1_tmp=float(j[order[i]["O1"]].split()[-1])
            O2_tmp=float(j[order[i]["O2"]].split()[-1])
            O1.append(O1_tmp)
            O2.append(O2_tmp)
        O1_mean.append(np.mean(O1_tmp))
        O2_mean.append(np.mean(O2_tmp))
    for i,item in enumerate(wind):
        Ochrg[item]={"O1":O1_mean[i],"O2":O2_mean[i]}
        
    return(Ochrg)
    
def mean_VDD_O(data,frag,order):
    O1_mean=[]
    O2_mean=[]
    wind=[]
    OVDD={}
    for i in data:
        wind.append(i)
        O1=[]
        O2=[]
        for j in data[i][frag]:
            O1_tmp=float(j[order[i]["O1"]-1].split()[-1])
            O2_tmp=float(j[order[i]["O2"]-1].split()[-1])
            O1.append(O1_tmp)
            O2.append(O2_tmp)
        #print("order",order[i]["O2"])
        #print("J",j[order[i]["O2"]-1].split()[-1])
        if i == "3.3":
            for j in O2:
                print(j)
        O1_mean.append(np.mean(O1_tmp))
        O2_mean.append(np.mean(O2_tmp))
    for i,item in enumerate(wind):
        OVDD[item]={"O1":O1_mean[i],"O2":O2_mean[i]}
        
    return(OVDD)
            
def write_data(data,frag,nat_frag,Ochrg):
        windows=[x for x in list(data.keys())]
        windows.sort()
        lines=["{:8s}{}".format("Atoms","  ".join(windows))]
        lines_O=["{:8s}{}".format("Atoms","  ".join(windows))]
        #lines.append(data[list(data.keys())[0]][frag][0][2].split()[1])
        for i in range(nat_frag):
            line=[]
            line.append("{:8s}".format(data[list(data.keys())[0]][frag][0][i+1].split()[1]))   
            for j in windows:
                line.append("{:8.4f}({:8.4f})".format(data[j][frag][-2][i],data[j][frag][-1][i]))
            lines.append(line)
        for i in range(2):
            line=[]
            line.append("O{}".format(i+1))
            for j in windows:
                line.append("{:8.4f}".format(Ochrg[j]["O{}".format(i+1)]))
            lines_O.append(line)
        with open('mean_charges_{}.dat'.format(frag),'w') as ofile:
            for item in lines:
                ofile.write("{}\n".format("".join(item)))
        with open('mean_charges_O_{}.dat'.format(frag),'w') as ofile:
            for item in lines_O:
                ofile.write("{}\n".format("".join(item)))

def write_VDD(data,frag,nat_frag,OVDD):
        windows=[x for x in list(data.keys())]
        windows.sort()
        lines=["{:8s}{}".format("Atoms","  ".join(windows))]
        lines_O=["{:8s}{}".format("Atoms","  ".join(windows))]
        #lines.append(data[list(data.keys())[0]][frag][0][2].split()[1])
        for i in range(1,nat_frag):
            line=[]
            line.append("{:8s}".format(data[list(data.keys())[0]][frag][0][i].split()[1]))   
            for j in windows:
                line.append("{:8.4f}({:8.4f})".format(data[j][frag][-2][i],data[j][frag][-1][i]))
            lines.append(line)
        for i in range(2):
            line=[]
            line.append("O{}".format(i+1))
            for j in windows:
                line.append("{:8.4f}".format(OVDD[j]["O{}".format(i+1)]))
            lines_O.append(line)
        with open('mean_VDD_{}.dat'.format(frag),'w') as ofile:
            for item in lines:
                ofile.write("{}\n".format("".join(item)))
        with open('mean_VDD_O_{}.dat'.format(frag),'w') as ofile:
            for item in lines_O:
                ofile.write("{}\n".format("".join(item)))
        
        
            


                



def main():
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--input' , dest='input',help='string to match input files')
    parser.add_argument('--VDD' , dest='VDD',help='VDD json file')
    parser.add_argument('--plot' , dest='plot', action='store_true',help='plot ')
    parser.add_argument('--all' , dest='all', action='store_true',help='plot all energies in one graph')
    parser.add_argument('--frag' , dest='frag', type=str, help='Fragment charges to extract')
    parser.add_argument('--natfrag' , dest='natfrag', type=int, help='Number of atoms in the fragment')
    
    args = parser.parse_args()
    
    if args.VDD:
        input,VDD,order=load_data(args.input,args.VDD)
        Ochrg=mean_charge_O(input,args.frag,order)
        OVDD=mean_VDD_O(VDD,args.frag,order)
        data=mean_charge(input,args.frag)
        VDD=mean_VDD(VDD,args.frag)
        write_data(data,args.frag,args.natfrag,Ochrg)
        write_VDD(VDD,args.frag,args.natfrag,OVDD)
    else:
        input,order=load_data(args.input,args.VDD)
        Ochrg=mean_charge_O(input,args.frag,order)
        data=mean_charge(input,args.frag)
        write_data(data,args.frag,args.natfrag,Ochrg)
    
   
   
                
       

if __name__=="__main__":
    main()
    