import numpy as np
import glob
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

def extract_data(input):
    data={}
    #print(glob.glob('*.metadynLog'))
    for i in glob.glob('{}'.format(input)):
        with open(i,'r') as ifile:
            lines=ifile.readlines()
        labx=lines[7].split(',')[0][0]
        laby=lines[7].split(',')[0][1]
        c1=[float(x) for x in lines[7].split(',')[1:-1]]
        c2=c1.copy()
        grid=[[float(x) for x in y.split(',')[1:-1]] for y in lines[8:]]
        grid=np.transpose(grid)
    return(c1,c2,grid,labx,laby)

def plot2d(x,y,value,file,labx,laby,cmap):
    fig=plt.figure(figsize=(16,10),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    mpl.rc('font', **font)
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #kjmol/plot
    lev=range(-20,22,1)
    plt.contourf(x, y,value,lev,cmap=cmap)
    ###
    
    #kcal/mol plot
    #lev=[x/2 for x in range(0,21,1)]
    #val_kcal=[x/4.184 for x in value]
    #plt.contourf(x, y,val_kcal,lev,vmin=0,vmax=10,cmap=cmap)
    
    plt.xlabel(labx)
    plt.ylabel(laby)
    #bounds=[1,2,3,4]
    #cbarkcal
    #cbar=plt.colorbar(label="$\Delta A\ (kcal\ mol^{-1})$",ticks=range(0,11,1))
    #cbar.ax.set_ylim(0,10)
    #cbar kj/mol
    cbar=plt.colorbar(label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(-20,22,2))
    cbar.ax.set_ylim(-20,20)
    #for i in minpath:
    #    plt.scatter(x[i[0]],y[i[1]],color='black')
    #plt.scatter(x[minima[0]],y[minima[1]])
    #plt.scatter(x[maxima[0]],y[maxima[1]],color='red')
    #plt.xlim([0.75,3.25])
    #plt.ylim([0.75,3.25])
    plt.savefig('{}.png'.format(file),format='png')
            

c1,c2,grid,labx,laby=extract_data(sys.argv[1])
#print(len(c1),len(c2),len(grid[0][0]))
file="NICS_plot"
cmap_active='jet'
plot2d(c1,c2,grid,file,labx,laby,cmap_active)
#print(grid[0])
            
"""
            
        for line in lines:
            data[i]['cv'].append(float(line.split()[1]))
        pos=i.split(sep)[field]
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
"""