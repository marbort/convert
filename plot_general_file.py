import numpy as np
import matplotlib.pyplot as plt
import os
import argparse




def extract_data(input):
    data={}
    filename=os.path.splitext(input)[0]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    header=lines[0].split()
    units=lines[1].split()
    legend=lines[2].split()
    for i,var in enumerate(legend):
            data[var]=[float(x.split()[i]) for x in lines[3:]]
    print(data,header)
    return(data,filename,units,legend,header)

def plot_data(data,filename,units,legend,header,log,x=None,y=None,wd=9,hg=9):
    if log:
        if x is not None:
            for i,item in enumerate(x):
                plt.rc('font', size=24)
                fig=plt.figure(figsize=(wd,hg),dpi=150)
                plt.plot([np.log10(x) for x in data[list(data.keys())[item]]],[np.log10(x) for x in data[list(data.keys())[y[i]]]])
                plt.xlabel(header[item]+" ({})".format(units[i]))
                plt.ylabel(header[y[i]]+" ({})".format(units[y[i]]))
                plt.xticks(ticks=[np.log10(x) for x in data[list(data.keys())[item]]],labels=data[list(data.keys())[item]])
                plt.yticks(ticks=[np.log10(x) for x in data[list(data.keys())[y[i]]]],labels=data[list(data.keys())[y[i]]])
                plt.tight_layout()
                plt.savefig(filename+"_Plot_"+list(data.keys())[y[i]]+".png",format='png')
                
        else:
            fig=plt.figure(figsize=(wd,hg),dpi=150)
            #plt.gca().set_aspect('equal')
            for i in range(1,len(data)):
                plt.rc('font', size=24)
                plt.plot([np.log10(x) for x in data[list(data.keys())[0]]],[np.log10(x) for x in data[list(data.keys())[i]]])
                plt.xlabel(header[0])
                plt.ylabel(header[1]+" ({})".format(units[1]))
                
            plt.xticks(ticks=[np.log10(x) for x in data[list(data.keys())[0]]],labels=["{:.0f}".format(x) for x in data[list(data.keys())[0]]])
            plt.yticks(ticks=[np.log10(x) for x in data[list(data.keys())[i]]],labels=["{:.1f}".format(x) for x in data[list(data.keys())[i]]])
            plt.legend(legend[1:])
            plt.tight_layout()
            plt.savefig(filename+"_Plot.png",format='png')
    else:
        if x is not None:
            for i,item in enumerate(x):
                plt.rc('font', size=24)
                fig=plt.figure(figsize=(wd,hg),dpi=150)
                plt.plot(data[list(data.keys())[item]],data[list(data.keys())[y[i]]])
                plt.xlabel(header[item]+" ({})".format(units[i]))
                plt.ylabel(header[y[i]]+" ({})".format(units[y[i]]))
                plt.tight_layout()
                plt.savefig(filename+"_Plot_"+list(data.keys())[y[i]]+".png",format='png')
                
        else:
            fig=plt.figure(figsize=(wd,hg),dpi=150)
            #plt.gca().set_aspect('equal')
            for i in range(1,len(data)):
                plt.rc('font', size=24)
                plt.plot(data[list(data.keys())[0]],data[list(data.keys())[i]])
                plt.xlabel(header[0])
                plt.ylabel(header[1]+" ({})".format(units[1]))
            
            plt.legend(legend[1:])
            plt.tight_layout()
            plt.savefig(filename+"_Plot.png",format='png')
    
    
             
            

def main():
    parser = argparse.ArgumentParser(description='Plot data')


    parser.add_argument('--input', dest='input', default=None,
                    type=str, help='input file')

    parser.add_argument('--xcols', dest='xcols', default=None,nargs='+',
                    type=int, help='x columns')
    parser.add_argument('--ycols', dest='ycols', default=None,nargs='+',
                    type=int, help='y columns')
    parser.add_argument('--log', dest='log', default=False,action='store_true',
                     help='set log plot')
    args = parser.parse_args()
    
    data,filename,units,legend,header=extract_data(args.input)
    plot_data(data,filename,units,legend,header,args.log,args.xcols,args.ycols)

if __name__ == "__main__":
   main()
    
