import plumed
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob

    




def extract_cvs(input,cvs):
    data={}
    files=glob.glob('*/{}'.format(input))
    files.sort()
    for i in files:
        window=i.split('/')[-2]
        with open(i,'r') as ifile:
            lines=ifile.readlines()
        cols=[i-2 for i,x in enumerate(lines[0].split()) if x in cvs ]
        vals=[[float(x.split()[j]) for x in lines[1:]] for j in ([0]+cols)]
        data[window]=vals
    return(data)

def plot_cvs(data,cvs,cvmin,cvmax,timeseries=False,wd=15,hg=10):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    name=""
    print("0")
    for i,window in enumerate(data):
        print("0")
        plt.subplot(len(data)//4+1,4,i+1)    
        if timeseries:
            print(0)
            name="time_cvs"
            for i in range(1,len(data[window])):
                plt.plot(data[window][0],data[window][i],label=cvs[i-1])
            plt.xlabel('time')
            plt.ylabel('CV')
            plt.ylim(min(cvmin),max(cvmax))
            plt.title(window)
            plt.legend()
        else:
            print(1)
            name="cvs"
            plt.plot(data[window][1],data[window][2],'o',markersize=1)
            plt.title(window)
            plt.xlim(cvmin[0],cvmax[0])
            plt.ylim(cvmin[1],cvmax[1])
            plt.xlabel(cvs[0])
            plt.ylabel(cvs[1])
    plt.tight_layout()
    plt.savefig('{}_{}_{}.png'.format(name,cvs[0],cvs[1]),format='png')
    #return(data)
        

    




def main():
    parser = argparse.ArgumentParser(description='Format colvar file for use with WHAM')

    parser.add_argument('--input',dest='input',type=str,default='COLVAR',help='full COLVAR file')
    parser.add_argument('--min',dest='min',type=float,default=0,help='min value of CV',nargs='+')
    parser.add_argument('--max',dest='max',type=float,default=0,help='max value of CV', nargs='+')
    parser.add_argument('--timeseries',dest='ts',action='store_true',help='Plot each CV vs time')
    

    parser.add_argument('--cv', dest='cv', default="",nargs='+',
                        type=str, help='define CVs to be plotted')
    args = parser.parse_args()
    
    data=extract_cvs(args.input,args.cv)
    #print(data['window_2.0'][0])
    plot_cvs(data,args.cv,args.min,args.max,args.ts)
 
    
    
if __name__=="__main__":
    main()




    
    
    
