import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl


def extract_data(input,split,sep,field):
    data={}
    #print(glob.glob('*.metadynLog'))
    for i in glob.glob('{}'.format(input)):
        with open(i,'r') as ifile:
            lines=ifile.readlines()
            data[i]={'cv':[]}
        for line in lines:
            data[i]['cv'].append(float(line.split()[1]))
        pos=i.split('-')[0].split(sep)[field]
        with open("cp2k_{}.inp".format(pos)) as cpinp:
            lines=cpinp.readlines()
            for line in lines:
                if 'K [' in line:
                    kappa=line.split()[-1]
        data[i]['pos']=float(pos.split('_')[0])
        data[i]['kappa']=float(kappa)
    for i in data:
        bins=100
        hist,edges=np.histogram(data[i]['cv'],bins=bins,density=True)
        data[i]['hist']=hist
        data[i]['edges']=edges[:-1]
    
        
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

def plot_data(data,fac,label,plot_gauss,wd=15,hg=10):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    mpl.rc('font', **font)
    for j,i in enumerate(data):
        kt=2.479
        mu=data[i]['pos']
        sigma=np.sqrt(kt/(data[i]['kappa']))
        gauss_range=np.linspace((mu-0.5)*fac,(mu+0.5)*fac,100)
        gauss=1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (gauss_range - mu)**2 / (2 * sigma**2) )
        #print(np.exp((data[i]['edges'] - mu)**2 / (2 * sigma)))
        #gauss=[len(data[i]['hist'])*np.exp(-data[i]['kappa']*(x-data[i]['pos'])**2) for x in data[i]['edges']*fac]
        width=(max(data[i]['edges'])*fac-min(data[i]['edges'])*fac)/len(data[i]['edges'])
        plt.bar([x*fac for x in data[i]['edges']],data[i]['hist'],width=width,alpha=0.3,color="C{}".format(j))
        if plot_gauss:
            plt.plot(gauss_range,gauss,color="C{}".format(j))
    if label==False:
        plt.legend([x for x in data],bbox_to_anchor=(1,1))
    plt.xlabel("CV")
    plt.ylabel("Counts")
    plt.tight_layout()
    plt.savefig('histo_all_plot.png',format='png')


def main():
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--split' , dest='split',action='store_true',help='output single histo files')
    parser.add_argument('--input' , dest='input',help='string to match input files')
    parser.add_argument('--plot' , dest='plot', action='store_true',help='plot all histograms in one graph')
    parser.add_argument('--nolabel' , dest='nolabel', action='store_true',help='do not plot legend')
    parser.add_argument('--gauss' , dest='gauss', action='store_true',help='plot theoretical gaussians')
    parser.add_argument('--wd' , dest='wd', help='plot width',default=None)
    parser.add_argument('--hg' , dest='hg', help='plot height',default=None)
    parser.add_argument('--fac' , dest='fac', help='CV conversion factor',default=1, type=float)
    parser.add_argument('--sep' , dest='sep',type=str,help='file name separator to extract position of restraint. Default -',default='-')
    parser.add_argument('--field' , dest='field',type=int,help='field number where position of restraint is in file name. Default 4.',default=4)
    
    
   
    args = parser.parse_args()
    data=extract_data(args.input,args.split,args.sep,args.field)
    if args.plot:
        plot_data(data,args.fac,args.nolabel,args.gauss)
        

if __name__=="__main__":
    main()
    