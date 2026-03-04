import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse


def plot2d(inputs,labx,laby,labels):
    fig=plt.figure(figsize=(16,12),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 46}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
   
    pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]
    xmin=0
    xmax=0
    ymin=0
    ymax=0
    for i,item in enumerate(inputs):
       data=np.loadtxt(item,unpack=True)
       plt.plot(data[0],data[1],'o-',color=pres_colors[i],label=labels[i])
       if max(data[0])>xmax:
           xmax=max(data[0])
       if min(data[0])<xmin:
           xmin=min(data[0])
       if max(data[1])>ymax:
           ymax=max(data[1])
       if min(data[1])<ymin:
           ymin=min(data[1])
    
    
    
    plt.xlim([0,2])
    plt.ylim([0.2,1.0])
    plt.xlabel(labx)
    plt.ylabel(laby)
    plt.legend()
    plt.tight_layout()
    plt.savefig('Compare_minpaths.png',format='png')

def main():
    parser = argparse.ArgumentParser(description='Plot minimum paths')
    parser.add_argument('--inputs', type=str, help='Input files', nargs='+', required=True)
    parser.add_argument('--labx', type=str, help='Label for x-axis',default='CV1')
    parser.add_argument('--laby', type=str, help='Label for y-axis',default='CV2')
    parser.add_argument('--labels', type=str, help='Labels for the paths', nargs='+', required=True)
    args = parser.parse_args()
    
    plot2d(args.inputs,args.labx,args.laby,args.labels)
    
if __name__ == "__main__":
    main()
