import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import scipy as sp
import papersize as psz
from decimal import *

#clrs=["#f21d2f","#2104d9","#f2ae2e"]
#clrs=['#5764D9','#D97C57','#86D957']
clrs=['#0000FF','#FF0000','#000000']+list(mpl.colors.TABLEAU_COLORS.keys())
mpl.rcParams['axes.linewidth'] = 3
lnwd=5
font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth']=3 

def set_size(paper,fraction):
    width=psz.convert_length(psz.parse_papersize(paper)[1],'pt','in')/fraction
    height=width/Decimal(1.5)
    return(float(width),float(height))
    

def extract_data(file):
    data={}
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        data['file']=ifile.name
    for line in lines:
        if line[0] == '#':
            continue
        if 'title' in line:
            data['title']=" ".join(line.split()[2:len(line.split())])
        if 'xaxis  label' in line:
            data['xlabel']=" ".join(line.split()[3:len(line.split())])
        if 'yaxis  label' in line:
            data['ylabel']=" ".join(line.split()[3:len(line.split())])
        if 's0 legend' in line:
            data['label']=" ".join(line.split()[3:len(line.split())])
    print("Found data with legend: {}".format(data['label']))
    data['x']=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
    data['y']=[float(x.split()[1]) for x in lines if "#" not in x if "@" not in x]
    return(data)

def split_data(data,min,max):
        data_split_x=[x for x in data['x'] if min < float(x) < max]
        data_split_y=[data['y'][y] for y,x in enumerate(data['x']) if min < float(x) < max]
        data_split=[data_split_x,data_split_y]
        rest=[[],[]]
        for i,val in enumerate(data['x']):
            if float(val) < min:
                rest[0].append(min+(min-val))
                rest[1].append(data['y'][i])
            elif float(val) > max:
                rest[0].append(max-(val-max))
                rest[1].append(data['y'][i])
        rest[0],rest[1]=zip(*sorted(zip(rest[0], rest[1])))
        mean=[[(x+rest[0][i])/2 for i,x in enumerate(data_split[0])],[(y+rest[1][j]) for j,y in enumerate(data_split[1])]]
        return(data_split,rest,mean)
    
def stack_data(data,magnitude):
    data_stack=[]
    for i,entry in enumerate(data):
        entry['y']=[y+magnitude*(len(data)-i) for y in entry['y']]
        data_stack.append(entry)
    return(data_stack)

def plot_data(data,files,legend,limx,limy,shade_min,shade_max,transparent,paper,fraction,clrs=clrs):
    #wd,hg=set_size(paper,fraction)
    wd=15
    hg=15
    lnwd=5
    name="_".join(files).replace('_density','')
    if transparent:
        fig=plt.figure(figsize=(wd,hg),dpi=150,frameon=False)
    else:
        fig=plt.figure(figsize=(wd,hg),dpi=150,facecolor='red')
        fig.patch.set_alpha(0.0)
        ax=plt.subplot()
        ax.patch.set_facecolor('white')
    
    
    for x,i in enumerate(data):
        plt.plot(i['x'],i['y'],'-',label=i['label'][1:-1],linewidth=lnwd,color=clrs[x])
    if shade_min:
        for k,val in enumerate(shade_min):
            plt.axvspan(val, shade_max[k], color='#989898', alpha=0.5, lw=0)
    #plt.text(0.25,0.05,"DES",horizontalalignment='center', transform = ax.transAxes)
    #plt.text(0.75,0.05,"THF",horizontalalignment='center', transform = ax.transAxes)
    #plt.xticks(np.arange(int(min(data[0]['x'])),int(max(data[0]['x'])),2))
    if legend:
        if legend[0] == "file":
            plt.legend()
        else:
            plt.legend(args.legend,loc='center right')
            #plt.legend(args.legend,loc='best')
            
    if limx:
        plt.xlim(limx)
    if limy:
        plt.ylim(limy)
    plt.xlabel("{}".format(data[0]['xlabel'][1:-1].replace('\S','$^{').replace('\\N','}$').replace('\\','}$')))
    plt.ylabel("{}".format(data[0]['ylabel'][1:-1].replace('\S','$^{').replace('\\N','}$')))
    plt.tight_layout()
    plt.savefig('{}_plot.png'.format(name),format='png',facecolor=fig.get_facecolor(),transparent=transparent)

def plot_data_mean(data,mean,files,legend,limx,limy,shade_min,shade_max,transparent,wd=15,hg=10,clrs=clrs):
    lnwd=5
    name="_".join(files).replace('_density_mean','')
    if transparent:
        fig=plt.figure(figsize=(wd,hg),dpi=300,frameon=False)
    else:
        fig=plt.figure(figsize=(wd,hg),dpi=300,frameon=False)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth']=3
    ax=plt.subplot()
    for i,x in enumerate(data):
        plt.plot(mean[i][2][0],mean[i][2][1],'-',label=x['label'][1:-1],linewidth=lnwd,color=clrs[i])
    #if shade_min:
     #   for k,val in enumerate(shade_min):
      #      plt.axvspan(val, shade_max[k], color='#989898', alpha=0.5, lw=0)
    #plt.text(0.25,0.05,"DES",horizontalalignment='center', transform = ax.transAxes)
    #plt.text(0.75,0.05,"THF",horizontalalignment='center', transform = ax.transAxes)
    #plt.xticks(np.arange(int(min(data[0]['x'])),int(max(data[0]['x'])),2))
    if legend:
        if legend == "file":
            plt.legend()
        else:
            plt.legend(args.legend,loc='best')
    if limx:
        plt.xlim(limx)
    if limy:
        plt.ylim(limy)
    plt.xlabel("{}".format(data[0]['xlabel'][1:-1].replace('\S','$^{').replace('\\N','}$').replace('\\','}$')))
    plt.ylabel("{}".format(data[0]['ylabel'][1:-1].replace('\S','$^{').replace('\\N','}$')))
    plt.savefig('{}_plot_mean.png'.format(name),format='png',transparent=True,facecolor=None,edgecolor=None)

def plot_integrated_density(data,files,time,wd=15,hg=10):
    name="_".join(files).replace('_density','')
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    DES_int=[x[0] for x in data]
    THF_int=[x[1]+x[2] for x in data]
    plt.plot([i*time for i in range(len(data))],[x-min(DES_int) for x in DES_int],label="DES")
    plt.plot([i*time for i in range(len(data))],[x-min(THF_int) for x in THF_int],label="THF")
    plt.xlabel("Time / ns")
    plt.ylabel("$\Delta N$ AcPh")
    plt.legend()
    plt.savefig('{}_Int_density_plot.png'.format(name),format='png')

def integrate_density(data,min,max):
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    """
    for i,item in enumerate(data['x']):
        if min <= item <= max:
            x1.append(item)
            y1.append(data['y'][i])
        else:
            x2.append(item)
            y2.append(data['y'][i])
    integral_1=np.trapz(y1,[x for x in x1])
    integral_2=np.trapz(y2,[x for x in x2])
    """     
    x1=[(x,i) for i,x in enumerate(data['x']) if min <= x <= max]
    y1=[data['y'][x[1]] for x in x1 ]
    x2=[(x,i) for i,x in enumerate(data['x']) if x < min ]
    y2=[data['y'][x[1]] for x in x2 ]
    x3=[(x,i) for i,x in enumerate(data['x']) if x > max ]
    y3=[data['y'][x[1]] for x in x3 ]
    integral_1=np.trapz(y1,[x[0] for x in x1])
    integral_2=np.trapz(y2,[x[0] for x in x2])
    integral_3=np.trapz(y3,[x[0] for x in x3])
    return(integral_1,integral_2,integral_3)

def plot_cumul_int(data,clrs,lw,wd=15,hg=10):
    name="_".join(files).replace('_cint','')
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)

    cint=[sp.integrate.cumulative_trapezoid(entry['y'],entry['x']) for entry in data]
    for i,entry in enumerate(data):
        print(entry)
        plt.plot(entry['x'][:-1],[cint[i][j]/x for j,x in enumerate(entry['x'][:-1])],linewidth=lw,color=clrs[i],label=entry['file'])
    #plt.xlim(6,19)
    plt.xlabel("Z / nm")
    plt.ylabel("Integrated number density")
    plt.legend()
    plt.savefig('{}_Int_density_cumul_plot.png'.format(name),format='png')
    
def plot_interval_integral(data,area,labels,limx,limy,wd=15,hg=10):
    name="_".join(files).replace('_iint','')
    intervals=[]
    for entry in data:
        print(entry[0])
        intg_interval=[np.trapz([entry[0][1][val],entry[0][1][val+1]],[entry[0][0][val],entry[0][0][val+1]]) for val in range(len(entry[0][0])-1)]
        intervals.append([entry[0][0],intg_interval])
    fig=plt.figure(figsize=(wd,hg),dpi=150,frameon=False)
    for i,entry in enumerate(intervals):
        plt.plot(entry[0][:-1],[entry[1][i]*area for i,x in enumerate(entry[0][:-1])],linewidth=lnwd,color=clrs[i])
    if labels:
        plt.legend(labels)
    plt.xlim(limx)
    plt.ylim(limy)
    plt.xlabel("Z / nm")
    plt.ylabel("Average number of molecules")
    plt.savefig('{}_Int_density_interval_plot.png'.format(name),format='png',transparent=True)


        
            
    
    
    
    
    

def integrate_interface_density(data,min,max):
    x1=[]
    y1=[]
    len_interface=sum([abs(x-max[i]) for i,x in enumerate(min)])
    """
    for i,item in enumerate(data['x']):
        if min <= item <= max:
            x1.append(item)
            y1.append(data['y'][i])
        else:
            x2.append(item)
            y2.append(data['y'][i])
    integral_1=np.trapz(y1,[x for x in x1])
    integral_2=np.trapz(y2,[x for x in x2])
    """
    for i,ex in enumerate(min):     
        x1.append([(x,j) for j,x in enumerate(data['x']) if ex <= x <= max[i]])
    for k in x1:
        y1.append([data['y'][x[1]] for x in k]) 
    integral=[np.trapz(y1[i],[x[0] for x in k]) for i,k in enumerate(x1)]
    sum_integral=sum(integral) 
    return(sum_integral,len_interface)


parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--input', dest='input',nargs='+')
parser.add_argument('--min_intg', dest='min_intg',type=float)
parser.add_argument('--max_intg', dest='max_intg',type=float)
parser.add_argument('--shade_min', dest='shade_min',type=float, nargs='+')
parser.add_argument('--shade_max', dest='shade_max',type=float, nargs='+')
parser.add_argument('--integrate', dest='intg',action='store_true')
parser.add_argument('--trans', dest='trans',action='store_true')
parser.add_argument('--cumint', dest='cumint',action='store_true')
parser.add_argument('--legend', dest='legend',nargs='+')
parser.add_argument('--xlim', dest='xlim',type=float,nargs=2)
parser.add_argument('--ylim', dest='ylim',type=float,nargs=2)
parser.add_argument('--stack', dest='stack',type=float,default=None)
parser.add_argument('--split', dest='split',type=float,nargs=2)
parser.add_argument('--area', dest='area',type=float)

args = parser.parse_args()
lw=3
data=[]
inputs=[]
integrals=[]
split_mean=[]
files=[]
for input in args.input:
    inputs+=glob.glob(input)
for input in inputs:
    files.append(os.path.basename(input).split('.')[0])
    data.append(extract_data(input))
    if args.split:
        split_mean.append(split_data(data[-1],args.split[0],args.split[1]))

if args.stack:
    data=stack_data(data,args.stack)

if args.intg:
    for i in data:
        integrals.append(integrate_density(i,args.min_intg,args.max_intg))
        print(i['file'])
        print("Integral_1: {:.2f} Integral_2: {:.2f} Integral_3: {:.2f} Min: {} Max: {}".format(integrals[-1][0],integrals[-1][1],integrals[-1][2],args.min_intg,args.max_intg))
        print(np.trapz(i['y'],i['x']))
        integral_interface,len_interface=integrate_interface_density(i,args.shade_min,args.shade_max)
        print(f"Integrated density at the interface(s):{integral_interface:.4f} nm-3")
        print(f"Total length of the interface(s):{len_interface:.4f} nm")
        print("###########################")
    plot_integrated_density(integrals,files,50)
if args.cumint:
    plot_cumul_int(data,clrs,lw)
    #print(split_mean)
    plot_interval_integral(split_mean,args.area,args.legend,args.xlim,args.ylim)


plot_data(data,files,args.legend,args.xlim,args.ylim,args.shade_min,args.shade_max,args.trans,'a4',3)
if args.split:
    plot_data_mean(data,split_mean,files,args.legend,[6,19],args.ylim,args.shade_min,args.shade_max,args.trans)


        
            
            