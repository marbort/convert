import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import glob
import matplotlib as mpl
from matplotlib import font_manager


# Specify the path to your custom font file
custom_font_path = "/home/marco/.fonts/f/Formular.ttf"

# Register the custom font
prop = font_manager.FontProperties(fname=custom_font_path)

# Use the custom font in your plot
plt.rcParams["font.family"] = prop.get_name()

  
  

def plot_contacts(paths,labels,title,fig_hist=None,fig_tseries=None,index=0,groups=1):
    font = {"weight": "normal", "size": 24}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=4, marker="o", markersize=6)
    mpl.rcParams["axes.linewidth"] = 3
    clr=['g','k','m','c','r','y']
    if fig_hist is None:
        fig_hist=plt.figure(figsize=(16,9))
    if fig_tseries is None:
        fig_tseries=plt.figure(figsize=(16,9))
    rows=int(np.sqrt(len(paths)))
    if np.sqrt(len(paths))-rows > 0:
        cols=rows + 1
    else:
        cols=rows
    for i,path in enumerate(paths):
        dihedrals = np.loadtxt(path,unpack=True,comments="#")#[float(line.split()[1]) for line in f if not line.startswith('#')]
        #print(dihedrals)
        
        plt.figure(fig_hist.number)
       
        plt.subplot(rows,cols,i+1)
        plt.bar([(groups+1)*i+index for i in range(len(dihedrals[0]))], dihedrals[1], width=1-.05, alpha=0.5, label=path,color=clr[index],
                edgecolor='black',linewidth=0.5)
        plt.xticks([(groups+1)*i-groups//2+index for i in range(len(dihedrals[0]))],dihedrals[0])
        plt.legend(labels,loc='best')
        plt.title(os.path.basename(path))
       
    plt.suptitle(title)
    plt.tight_layout()
    return(fig_hist,fig_tseries)
    

def main():
    parser = argparse.ArgumentParser(description='Calculate dihedral statistics')
    
    parser.add_argument('dihedral_files',  nargs='+', help='Folder containing the dihedral files to analyze')
    parser.add_argument('--pattern',  help='Pattern of data file names', default=None, type=str)
    parser.add_argument('--labels', type=str, help='Legend labels', nargs='+', default=None)
    parser.add_argument('--title', type=str, help='Figure title', default=None)
    
    
    args=parser.parse_args()
    
    paths = args.dihedral_files
   
    fig_hist=None
    fig_tseries=None
    for i,path in enumerate(paths):
        
        files=glob.glob(os.path.join(path,args.pattern))
        print(f"Found {len(files)} files in {path}")
        fig_hist,fig_tseries=plot_contacts(files,args.labels,args.title,fig_hist,fig_tseries,index=i,groups=len(paths))
    
    
    
    plt.figure(fig_hist.number)
    plt.savefig('contacts_histogram.png',dpi=300)
    #plt.figure(fig_tseries.number)
    #plt.savefig('dihedral_timeseries.png',dpi=300)
   
    
        
if __name__ == '__main__':
    main()
        