import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import matplotlib.gridspec as gridspec
from matplotlib import font_manager



# Specify the path to your custom font file
custom_font_path = "/home/marco/.fonts/f/Formular.ttf"

# Register the custom font
prop = font_manager.FontProperties(fname=custom_font_path)

# Use the custom font in your plot
plt.rcParams["font.family"] = prop.get_name()

def calc_dihedral_angles(file,topology,idxs):
    u = mda.Universe(file,topology)
    for item in idxs:
        atom1 = u.select_atoms(f"index {idxs[0]}")[0]
        atom2 = u.select_atoms(f"index {idxs[1]}")[0]
        atom3 = u.select_atoms(f"index {idxs[2]}")[0]
        atom4 = u.select_atoms(f"index {idxs[3]}")[0]

        # Initialize a list to store dihedral angles
        dihedral_angles = []

        # Iterate over the trajectory and calculate dihedral angles
        for ts in u.trajectory:
            # Get the positions of the selected atoms
            pos1 = atom1.position
            pos2 = atom2.position
            pos3 = atom3.position
            pos4 = atom4.position
    
            # Calculate the dihedral angle
            angle = calc_dihedrals(pos1, pos2, pos3, pos4, box=ts.dimensions)
        
            # Store the angle
            dihedral_angles.append(angle)  # angles is an array, get the first element
    
    with open('pippo.dat','w') as ofile:
        for i,dih in enumerate(dihedral_angles):
            ofile.write(f"{i} {dih*57.03}\n")
            
    
    return dihedral_angles
    
    


def calc_dihedral_statistics(file):
        dihedrals = []
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    try:
                        dihedrals.append(float(line.split()[1]))
                    except:
                        pass
        dihedrals = np.array(dihedrals)
        counts,bins,bars=plt.hist(dihedrals, bins=180,density=True)
        maxhist=np.max(counts)
        with open(os.path.join(os.path.dirname(file),"statistics.out"),"a") as ofile:
            ofile.write(f'Statistics for {file}\n')
            ofile.write(f'Mean: {np.mean(dihedrals)}\n')
            ofile.write(f'Std: {np.std(dihedrals)}\n')
            ofile.write(f'Min: {np.min(dihedrals)}\n')
            ofile.write(f'Max: {np.max(dihedrals)}\n')
            ofile.write(f'Count: {len(dihedrals)}\n\n')
        return(dihedrals, np.mean(dihedrals), np.std(dihedrals), np.min(dihedrals), np.max(dihedrals), len(dihedrals),maxhist)


def calc_dihedral_statistics_blocks(path,n):
    dihedrals = []
    with open(path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    try:
                        dihedrals.append(float(line.split()[1]))
                    except:
                        pass
    dihedrals = [dihedrals[i:i + len(dihedrals)//n] for i in range(0, len(dihedrals), len(dihedrals)//n)]
    dihedrals = np.array(dihedrals)
    
    #print('Statistics for', path)
    for i in range(len(dihedrals)):
        print('Block:', i)
        print('Mean:', np.mean(dihedrals[i]))
        print('Std:', np.std(dihedrals[i]))
        print('Min:', np.min(dihedrals[i]))
        print('Max:', np.max(dihedrals[i]))
        print('Count:', len(dihedrals[i]),'\n\n')
        
    return (dihedrals)

def plot_dihedrals(paths,means,stds,labels,title,fig_hist=None,fig_tseries=None,index=0,groups=1,ymax=0):
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
    print(len(paths),rows)
    if np.sqrt(len(paths))-rows > 0:
        cols=rows + 1
    else:
        cols=rows
    outer = gridspec.GridSpec(rows, cols, wspace=0.2, hspace=0.2)
    plt.figure(fig_hist.number)
    plt.suptitle(title)
    for k in range(len(paths)):
        inner = gridspec.GridSpecFromSubplotSpec( groups, 1,
                    subplot_spec=outer[k], wspace=0.2, hspace=0.2)
        
        abspath=os.path.abspath(paths[k])
        dihedrals = np.loadtxt(paths[k],unpack=True,comments="#")   
        plt.subplot(inner[index])
        counts,bins,bars=plt.hist(dihedrals[1], bins=180,range=(-180,180), alpha=0.8,color=clr[index],density=True,label=labels[index])
        if np.max(counts) > ymax:
            ymax=np.max(counts)
        plt.xlim(-180,180)
        plt.ylim(0,ymax)
        plt.legend(loc='upper right')
        if index == 0:
            plt.title(os.path.basename(paths[k]))
    
    plt.tight_layout()
    plt.figure(fig_tseries.number)
    for i,path in enumerate(paths):
        abspath=os.path.abspath(path)
        dihedrals = np.loadtxt(path,unpack=True,comments="#")#[float(line.split()[1]) for line in f if not line.startswith('#')]
        plt.subplot(rows,cols,i+1)
        plt.scatter(dihedrals[0],dihedrals[1], label=path, s=3,color=clr[index])
        if "RoG" in os.path.basename(path):
            plt.ylim(1,4)
        else:
            plt.ylim(-180,180)
        plt.title(os.path.basename(path))
        #t=plt.text(0.05,0.05+index/10,f"Mean: {means[abspath]:.2f} Std: {stds[abspath]:.2f}",transform=plt.gca().transAxes)
        #t.set_bbox(dict(facecolor='white', alpha=1,edgecolor=clr[index])) 
        #plt.fill_between(dihedrals[0],dihedrals[1]-stds[i],dihedrals[1]+stds[i],alpha=0.5)
        plt.legend(labels,loc='upper right')
        plt.tight_layout()
    return(fig_hist,fig_tseries)


def plot_dihedrals_same(paths,means,stds,labels,title,fig_hist=None,fig_tseries=None,index=0,groups=1,ymax=0):
    font = {"weight": "normal", "size": 24}
    mpl.rc("font", **font)
    mpl.rc("lines", linewidth=4, marker="o", markersize=6)
    mpl.rcParams["axes.linewidth"] = 3
    clr=['g','k','m','c','r','y']
    if fig_hist is None:
        fig_hist=plt.figure(figsize=(16,9))
    if fig_tseries is None:
        fig_tseries=plt.figure(figsize=(16,9))
    
    rows=len(paths)
    plt.figure(fig_hist.number)
    plt.suptitle(title)
    for k in range(len(paths)):
        
        abspath=os.path.abspath(paths[k])
        dihedrals = np.loadtxt(paths[k],unpack=True,comments="#")   
        plt.subplot(rows,1,k+1)
        counts,bins,bars=plt.hist(dihedrals[1], bins=180, alpha=0.8,color=clr[index],density=True,label=labels[index])
        if np.max(counts) > ymax:
            ymax=np.max(counts)
        plt.xlim(-180,180)
        plt.ylim(0,ymax)
        plt.legend(loc='upper right')
        if index == 0:
            plt.title(os.path.basename(paths[k]))
    
    plt.tight_layout()
    plt.figure(fig_tseries.number)
    for i,path in enumerate(paths):
        abspath=os.path.abspath(path)
        dihedrals = np.loadtxt(path,unpack=True,comments="#")#[float(line.split()[1]) for line in f if not line.startswith('#')]
        plt.subplot(rows,1,i+1)
        plt.scatter(dihedrals[0],dihedrals[1], label=path, s=3,color=clr[index])
        if "RoG" in os.path.basename(path):
            plt.ylim(1,4)
        else:
            plt.ylim(-180,180)
        plt.title(os.path.basename(path))
        #t=plt.text(0.05,0.05+index/10,f"Mean: {means[abspath]:.2f} Std: {stds[abspath]:.2f}",transform=plt.gca().transAxes)
        #t.set_bbox(dict(facecolor='white', alpha=1,edgecolor=clr[index])) 
        #plt.fill_between(dihedrals[0],dihedrals[1]-stds[i],dihedrals[1]+stds[i],alpha=0.5)
        plt.legend(labels,loc='upper right')
        plt.tight_layout()
    return(fig_hist,fig_tseries)
    

def main():
    parser = argparse.ArgumentParser(description='Calculate dihedral statistics')
    
    parser.add_argument('dihedral_files',  nargs='+', help='Folder containing the dihedral files to analyze')
    parser.add_argument('--pattern',  help='Pattern of data file names', default=None, type=str)
    parser.add_argument('--blocks', type=int, help='Number of blocks to split the dihedral data', default=1)
    parser.add_argument('--labels', type=str, help='Legend labels', nargs='+', default=None)
    parser.add_argument('--title', type=str, help='Title of the plot', default=None)
    parser.add_argument('--same', action='store_true', help='Plot all dihedrals in the same plot')
    
    
    args=parser.parse_args()
    
    paths = args.dihedral_files
    means={}
    stds={}
    fig_hist=None
    fig_tseries=None
    limy=[]
    for i,path in enumerate(paths):
        files=glob.glob(os.path.join(path,args.pattern))
        print(f"Found {len(files)} files in {path}")
        for file in files:
        
            dihedrals,mean,std,min_d,max_d,len_d,maxhist=calc_dihedral_statistics(file)
            #if args.blocks > 1:
            #    calc_dihedral_statistics_blocks(file,args.blocks)
            means[os.path.abspath(file)]=mean
            stds[os.path.abspath(file)]=std
            limy.append(maxhist)
    
    for i,path in enumerate(paths):
        files=glob.glob(os.path.join(path,args.pattern))
        if args.same:
            fig_hist,fig_tseries=plot_dihedrals_same(files,means,stds,args.labels,args.title,fig_hist,fig_tseries,index=i,
                                             groups=len(paths),ymax=np.max(limy))
        else:  
            fig_hist,fig_tseries=plot_dihedrals(files,means,stds,args.labels,args.title,fig_hist,fig_tseries,index=i,
                                                groups=len(paths),ymax=np.max(limy))
            
    print(means)
    
    
    
    plt.figure(fig_hist.number)
    plt.savefig('dihedral_histogram.png',dpi=300)
    plt.figure(fig_tseries.number)
    plt.savefig('dihedral_timeseries.png',dpi=300)
  
    
        
if __name__ == '__main__':
    main()
        