import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import os


def extract_data(input,split):
    data={}
    print(sorted(glob.glob('*/{}'.format(input))))
    for i in sorted(glob.glob('*/{}'.format(input))):
        #print(os.path.dirname(i))
        with open(i,'r') as ifile:
            lines=ifile.readlines()
            data[i]={'hist':[],'edges':[]}
        for line in lines[1:]:
            data[i]['hist'].append(float(line.split()[1]))
            data[i]['edges'].append(float(line.split()[0]))
        with open("{}/plumed.dat".format(os.path.dirname(i)),'r') as ifile:
            lines=ifile.readlines()
            for line in lines:
                if 'restraint:' in line:
                    data[i]['pos']=float(line.split()[-2].split('=')[-1])
                    data[i]['kappa']=float(line.split()[-1].split('=')[-1])
    if split:
        for i in data:
            colvar=i.split('-')[0].split('_')[-1]
            with open('{}_histo.dat'.format(colvar),'w') as ofile:
                ofile.write("""CPUS    Time 
#       s       
CPUS    Actual 
""")
                for j,line in enumerate(data[i]['hist']):
                    ofile.write('{:8.3f}{:8.3f}\n'.format(data[i]['edges'][j],line))
    else:
        with open('histo_all.dat','w') as ofile:
           ofile.write(" ".join(["{:s}1  {:s}2".format(x,x) for x in data])) 
           ofile.write("\n")
           ofile.write(" ".join(["{:s}3  {:s}4".format(x,x) for x in data])) 
           ofile.write("\n")  
           ofile.write(" ".join(["{:s}5  {:s}6".format(x,x) for x in data])) 
           ofile.write("\n")
           for i,line in enumerate(data[list(data.keys())[0]]['hist']):
               ofile.write("".join(["{:8.3f}{:8.3f}".format(data[x]['edges'][i],data[x]['hist'][i]) for x in data]))
               ofile.write("\n")
        
    return(data)

def plot_data(data,fac,xmin,xmax,wd=15,hg=10):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    for j,i in enumerate(data):
        kt=2.479
        mu=data[i]['pos']
        sigma=np.sqrt(kt/(data[i]['kappa']))
        print("{} Sigma= {}".format(i,sigma))
        gauss_range=np.linspace((mu-0.5)*fac,(mu+0.5)*fac,100)
        gauss=1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (gauss_range - mu)**2 / (2 * sigma**2) )
        width=(max(data[i]['edges'])*fac-min(data[i]['edges'])*fac)/len(data[i]['edges'])
        plt.bar([x*fac for x in data[i]['edges']],data[i]['hist'],width=width,alpha=0.3,color="C{}".format(j))
        plt.plot(gauss_range,gauss,color="C{}".format(j))
    plt.legend([x for x in data])
    plt.xlim(xmin,xmax)
    plt.savefig('histo_all_plot.png',format='png')


def main():
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--split' , dest='split',action='store_true',help='output single histo files')
    parser.add_argument('--input' , dest='input',help='string to match input files')
    parser.add_argument('--plot' , dest='plot', action='store_true',help='plot all histograms in one graph')
    parser.add_argument('--wd' , dest='wd', help='plot width',default=None)
    parser.add_argument('--hg' , dest='hg', help='plot height',default=None)
    parser.add_argument('--fac' , dest='fac', help='CV conversion factor',default=1, type=float)
    parser.add_argument('--min' , dest='min', help='x min',default=0, type=float)
    parser.add_argument('--max' , dest='max', help='x max',default=1, type=float)
    
    
   
    args = parser.parse_args()
    data=extract_data(args.input,args.split)
    if args.plot:
        plot_data(data,args.fac,args.min,args.max)
        

if __name__=="__main__":
    main()
    