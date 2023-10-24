import numpy as np
import glob
import argparse
import matplotlib
import matplotlib.pyplot as plt
import os
from operator import add
import json



def extract_data(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    data={}
    
    for line in lines[1:]:
        if float(line.split()[0].split('/')[-3].split('_')[-1]) in data:
            try:
                data[float(line.split()[0].split('/')[-3].split('_')[-1])]['E'].append(float(line.split()[-1]))   #Energy in kj/mol
            except:
                try:
                    print("{} not a float".format(line.split()[-1]))
                    data[float(line.split()[0].split('/')[-3].split('_')[-1])]['E'].append(float(line.split()[-2]))   #Energy in kj/mol
                except:
                    print("Energy not found {}".format(line.split()[-1]))
        else:
            data[float(line.split()[0].split('/')[-3].split('_')[-1])]={'E':[]}
            try:
                data[float(line.split()[0].split('/')[-3].split('_')[-1])]['E'].append(float(line.split()[-1]))   #Energy in kj/mol
            except:
                try:
                    print("{} not a float".format(line.split()[-1]))
                    data[float(line.split()[0].split('/')[-3].split('_')[-1])]['E'].append(float(line.split()[-2]))   #Energy in kj/mol
                except:
                    print("Energy not found {}".format(line.split()[-1]))
                    
            
            
    for i in data:
        data[i]['mean']=np.mean(data[i]['E'])
        data[i]['std']=np.std(data[i]['E'])
        data[i]['std_mean']=data[i]['std']/np.sqrt(len(data[i]['E']))
        try:
            with open('colvar_window_{}.xyz'.format(i),'r') as cvfile:
                lines=cvfile.readlines()                 
            data[i]['cv']=[float(x.split()[-1]) for x in lines[1:]]
        except:
            print('no colvar file colvar_{}.xyz'.format(i))
    #print(data[2.3]['cv'])
    return(data)

def calculate_strain(data,Efg1,Efg2):
    kcal_to_kj=4.184
    data['Strain.txt']={}
    data['DeltaE.txt']={}
    for window in data['Interaction.txt']:
        data['Strain.txt'][window]={}
        data['DeltaE.txt'][window]={}
        strain_fg1=[(x-Efg1)*kcal_to_kj for x in data['Fg1.txt'][window]['E']]
        strain_fg2=[(x-Efg2)*kcal_to_kj for x in data['Fg2.txt'][window]['E']]
        strain=list(map(add,strain_fg1,strain_fg2))
        DeltaE=list(map(add,strain,data['Interaction.txt'][window]['E']))
        strain_mean=np.mean(strain)
        DeltaE_mean=np.mean(DeltaE)
        strain_std=np.std(strain)
        DeltaE_std=np.std(DeltaE)
        data['Strain.txt'][window]['E']=strain
        data['Strain.txt'][window]['mean']=strain_mean
        data['Strain.txt'][window]['mean_fg1']=np.mean(strain_fg1)
        data['Strain.txt'][window]['mean_fg2']=np.mean(strain_fg2)
        data['Strain.txt'][window]['std']=strain_std
        data['DeltaE.txt'][window]['E']=DeltaE
        data['DeltaE.txt'][window]['mean']=DeltaE_mean
        data['DeltaE.txt'][window]['std']=DeltaE_std
    return(data)
    

def plot_data(data,input,wd=15,hg=10):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    plt.errorbar([float(list(data.keys())[x]) for x in  range(len(data)) ], [data[x]['mean'] for x in data], yerr=[data[x]['std'] for x in data])
    plt.savefig('{}_plot.png'.format(input),format='png')

def plot_data_all(data,xlab,wd=15,hg=10):
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    matplotlib.rc('font', **font)
    
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    for i in data:
        if i=="Fg1.txt" or i=="Fg2.txt":
            continue
        else:
            #plt.subplot(len(data),1,j+1)
            plt.errorbar([float(list(data[i].keys())[x]) for x in  range(len(data[i])) ], [data[i][x]['mean'] for x in data[i]], 
                            yerr=[data[i][x]['std'] for x in data[i]],label=i.replace('.txt',''))
    plt.legend(ncol=5)
    plt.ylabel("Energy ($Kj\ mol^{-1}$)")
    plt.xlabel(xlab)
    plt.savefig('all_plot.png',format='png')

def plot_data_ASA(data,fragments,xlab,inverted,wd=15,hg=10):
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    matplotlib.rc('font', **font)
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    clrs=["black","blue","red"]
    for j,i in enumerate(["Strain.txt","Interaction.txt","DeltaE.txt"]):
            #plt.subplot(len(data),1,j+1)
            cvsort=[float(list(data[i].keys())[x]) for x in  range(len(data[i])) ]
            cvsort.sort()
            plt.errorbar(cvsort, [data[i][x]['mean'] for x in cvsort], 
                            yerr=[data[i][x]['std'] for x in cvsort],label=i.replace('.txt',''),color=clrs[j])
            if fragments:
                if i=="Strain.txt":
                    plt.plot(cvsort, [data[i][x]['mean_fg1'] for x in cvsort],'--',
                                label='Fg1',color=clrs[j])
                    plt.plot(cvsort, [data[i][x]['mean_fg2'] for x in cvsort],'--',
                                label='Fg2',color=clrs[j])
    if inverted:
        plt.gca().invert_xaxis()
    plt.legend(ncol=5)
    plt.ylabel("Energy ($Kj\ mol^{-1}$)")
    plt.xlabel(xlab)
    plt.savefig('ASA_plot.png',format='png')

def plot_data_EDA(data,xlab,inverted,wd=15,hg=10):
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    matplotlib.rc('font', **font)
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    clrs=["blue","purple","green","brown"]
    for j,i in enumerate(["Interaction.txt","Pauli.txt","Velst.txt","OI.txt"]):
            #plt.subplot(len(data),1,j+1)
            cvsort=[float(list(data[i].keys())[x]) for x in  range(len(data[i])) ]
            cvsort.sort()
            plt.errorbar(cvsort, [data[i][x]['mean'] for x in cvsort], 
                            yerr=[data[i][x]['std'] for x in cvsort],label=i.replace('.txt',''),color=clrs[j])
            #plt.errorbar([float(list(data[i].keys())[x]) for x in  range(len(data[i])) ], [data[i][x]['mean'] for x in data[i]], 
            #                yerr=[data[i][x]['std'] for x in data[i]],label=i.replace('.txt',''),color=clrs[j])
    if inverted:
        plt.gca().invert_xaxis()
    plt.legend(ncol=5)
    plt.ylabel("Energy ($Kj\ mol^{-1}$)")
    plt.xlabel(xlab)
    plt.savefig('EDA_plot.png',format='png')

def plot_data_window(data,wd=15,hg=10):
    for i in data:
        fig=plt.figure(figsize=(wd,hg),dpi=150)
        for j in data[i]:
            print(j)
            plt.plot(data[i][j]['cv'],data[i][j]['E'],'o',markersize=3,)
        plt.savefig('{}_all.png'.format(i),format='png')
    
     #   print('no data')
                



def main():
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--input' , dest='input',help='string to match input files')
    parser.add_argument('--plot' , dest='plot', action='store_true',help='plot ')
    parser.add_argument('--all' , dest='all', action='store_true',help='plot all energies in one graph')
    parser.add_argument('--frags' , dest='frags', action='store_true',help='plot strain separately for each fragment')
    parser.add_argument('--ASA' , dest='ASA', action='store_true',help='plot ASA')
    parser.add_argument('--EDA' , dest='EDA', action='store_true',help='plot EDA')
    parser.add_argument('--invert' , dest='invert', action='store_true',help='invert x axis')
    parser.add_argument('--wd' , dest='wd', help='plot width',default=None)
    parser.add_argument('--hg' , dest='hg', help='plot height',default=None)
    parser.add_argument('--fac' , dest='fac', help='CV conversion factor',default=1, type=float)
    parser.add_argument('--fg1' , dest='fg1', help='Fragment 1 optimized energy',default=1, type=float)
    parser.add_argument('--fg2' , dest='fg2', help='Fragment 2 optimized energy',default=1, type=float)
    parser.add_argument('--xlab' , dest='xlab', help='X axis label',default="CV", type=str)
    args = parser.parse_args()
    
    data={}
    for file in glob.glob(args.input):
            data[file]=extract_data(file)

    print(max(data['Interaction.txt'][3.3]['E']))
    #data=calculate_strain(data,args.fg1,args.fg2)
    with open('data.json','w') as ofile:
        json.dump(data,ofile,indent=6)
    
        
        
        
        
    #print(data['Pauli.txt'])
    
    if args.plot:
        if args.all:
            plot_data_all(data,args.xlab)
            #plot_data_window(data)
        if args.ASA:
            data=calculate_strain(data,args.fg1,args.fg2)
            plot_data_ASA(data,args.frags,args.xlab,args.invert)
            with open('data.json','w') as ofile:
                json.dump(data,ofile,indent=6)
        if args.EDA:
            plot_data_EDA(data,args.xlab,args.invert)
        else:
            for i in data: 
                plot_data(data[i],i)
                
       

if __name__=="__main__":
    main()
    