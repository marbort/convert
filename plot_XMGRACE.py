import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


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

def plot_data(data,files,legend,limx,limy,shade_min,shade_max,wd=15,hg=10):
    name="_".join(files).replace('_density','')
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 26}
    mpl.rc('font', **font)
    ax=plt.subplot()
    for i in data:
        plt.plot(i['x'],i['y'],'-',label=i['label'][1:-1])
    if shade_min:
        for k,val in enumerate(shade_min):
            plt.axvspan(val, shade_max[k], color='#989898', alpha=0.5, lw=0)
    #plt.text(0.25,0.05,"DES",horizontalalignment='center', transform = ax.transAxes)
    #plt.text(0.75,0.05,"THF",horizontalalignment='center', transform = ax.transAxes)
    #plt.xticks(np.arange(int(min(data[0]['x'])),int(max(data[0]['x'])),2))
    if legend:
        if legend == "file":
            plt.legend()
        else:
            plt.legend(args.legend)#,loc='center right')
    if limx:
        plt.xlim(limx)
    if limy:
        plt.ylim(limy)
    plt.xlabel("{}".format(data[0]['xlabel'][1:-1].replace('\S','$^{').replace('\\N','}$').replace('\\','}$')))
    plt.ylabel("{}".format(data[0]['ylabel'][1:-1].replace('\S','$^{').replace('\\N','}$')))
    plt.savefig('{}_plot.png'.format(name),format='png')

def plot_integrated_density(data,files,time,wd=15,hg=10):
    name="_".join(files).replace('_density','')
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 26}
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
    


parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--input', dest='input',nargs='+')
parser.add_argument('--min_intg', dest='min_intg',type=float)
parser.add_argument('--max_intg', dest='max_intg',type=float)
parser.add_argument('--shade_min', dest='shade_min',type=float, nargs='+')
parser.add_argument('--shade_max', dest='shade_max',type=float, nargs='+')
parser.add_argument('--integrate', dest='intg',action='store_true')
parser.add_argument('--legend', dest='legend',nargs='+')
parser.add_argument('--xlim', dest='xlim',type=float,nargs=2)
parser.add_argument('--ylim', dest='ylim',type=float,nargs=2)

args = parser.parse_args()
data=[]
inputs=[]
integrals=[]
files=[]
for input in args.input:
    inputs+=glob.glob(input)
for input in inputs:
    files.append(os.path.basename(input).split('.')[0])
    data.append(extract_data(input))
if args.intg:
    for i in data:
        integrals.append(integrate_density(i,args.min_intg,args.max_intg))
        print("Integral_1: {:.2f}, Integral_2: {:.2f},Integral_3: {:.2f}, Min: {}, Max: {}".format(integrals[-1][0],integrals[-1][1],integrals[-1][2],args.min_intg,args.max_intg))
        print(np.trapz(i['y'],i['x']))
    plot_integrated_density(integrals,files,50)
plot_data(data,files,args.legend,args.xlim,args.ylim,args.shade_min,args.shade_max)


        
            
            