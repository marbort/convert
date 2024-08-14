import plot2d_reweight as rwg
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np
from scipy.optimize import minimize
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage 
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def plot_joint(x,y,value,file,labx,laby,cmap,minima,min_pt,cv1,cv2,y_max_1d=120):
    fig=plt.figure(figsize=(16,10),dpi=150)
    #gs = fig.add_gridspec(1, 2, width_ratios=[1, 10])
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 32}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX=int(y_max_1d)
    
    #plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    #kjmol/plot
    lev=range(0,MAX+5,5)
    custom_params = {
                    "axes.facecolor":"white",
                    "axes.edgecolor":"black",
                    "axes.spines.right": True, "axes.spines.top": True,
                    "axes.spines.left": True,"axes.spines.bottom": True,
                    "xtick.bottom":True,"xtick.major.size":10,"xtick.major.width":3,
                    "ytick.left":True,"ytick.major.size":10,"ytick.major.width":3,
                    "axes.linewidth" : 5, 
                    "grid.color":"black","grid.alpha":0.0,"grid.linestyle":":"}
    sns.set_theme(font='Formular',font_scale=4,context="paper",rc=custom_params)
    #sns.set_style("ticks",font="Formular",font_scale=5)
    
    
    custom_headers = ['C1', '$\Delta$A ($kJ\ mol^{-1})$', 'C3']
    try:
        dataset_CV1=pd.read_csv(cv1,delim_whitespace=True,comment="#",names=custom_headers)
        dataset_CV2=pd.read_csv(cv2,delim_whitespace=True,comment="#",names=custom_headers)
    except:
        dataset_CV1=pd.read_csv(cv1,delim_whitespace=True,comment="#",names=custom_headers[:-1])
        dataset_CV2=pd.read_csv(cv2,delim_whitespace=True,comment="#",names=custom_headers[:-1])
        
    #Create a seaborn joint grid
    g = sns.JointGrid(height=20,ratio=3,marginal_ticks=True,space=0.02)

    # Plot your seaborn joint plot
    #g.plot_joint(sns.scatterplot, color='blue')

    # Get the matplotlib axes
    ax = g.ax_joint
    

    # Now, you can plot your matplotlib plot using the ax
    #ax.plot(x_data, y_data, color='red')  # Example matplotlib plot
    CLines=ax.contour(x,y,value,levels=range(0,80,20),vmin=0,vmax=MAX,linewidths=2,colors='white')
    ax.clabel(CLines,levels=range(0,80,20), inline=True, fontsize=24,colors='white')
    contour=ax.contourf(x, y,value,lev,vmin=0,vmax=MAX,cmap=cmap)
    ax.set_xlabel(labx)
    ax.set_ylabel(laby)
    ax.set_ylim(bottom=min(y),top=max(y))
    ax.set_xlim(left=min(x),right=max(x))
    #ax.set_frame_on(True)
    ###
    
    #kcal/mol plot
    #lev=[x/2 for x in range(0,21,1)]
    #val_kcal=[x/4.184 for x in value]
    #plt.contourf(x, y,val_kcal,lev,vmin=0,vmax=10,cmap=cmap)
    
    #plt.xlabel(labx)
    #plt.ylabel(laby)
    plt.xticks(np.arange(min(x),max(x)+0.5,0.5))
    #bounds=[1,2,3,4]
    #cbarkcal
    #cbar=plt.colorbar(label="$\Delta A\ (kcal\ mol^{-1})$",ticks=range(0,11,1))
    #cbar.ax.set_ylim(0,10)
    #cbar kj/mol
    #cax = fig.add_subplot(gs[1])
    cbaxes = inset_axes(ax, width="30%", height="3%",  bbox_to_anchor=(0.01, 0.1, 1, 1),bbox_transform=ax.transAxes,loc=4)
    plt.setp(cbaxes.get_xticklabels(), backgroundcolor="white")
    #cbaxes.set_facecolor([1,1,1,1])
    cbar=fig.colorbar(contour,label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(0,MAX+MAX//10,MAX//2),cax=cbaxes, orientation='horizontal')
    cbar.ax.set_xlim(0,MAX)
    
    #for i in minpath:
    #    plt.scatter(x[i[0]],y[i[1]],color='black')
    #plt.scatter(x[minima[0]],y[minima[1]])
    #plt.scatter(x[maxima[0]],y[maxima[1]],color='red')
    #plt.xlim([0.75,3.25])
    #plt.ylim([0.75,3.25])
    if len(minima) > 0:
        min_crd=[]
        with open(minima,'r') as ifile:
            lines=ifile.readlines()
        for line in lines:
            min_crd.append((int(line.split()[1]),float(line.split()[5]),float(line.split()[-1])))
        for pt in min_crd:
            ax.scatter(pt[1],pt[2],color='white')
    good_minima=[]
    try:
        for i,item in enumerate(min_pt[0]):
            if value[item,min_pt[1][i]] > 60:
                pass
            else:
                good_minima.append((x[min_pt[1][i]],y[item],value[item][min_pt[1][i]]))
                #plt.scatter(x[min_pt[1][i]],y[item],color='red')
        with open('minima.dat','w') as ofile:
            for i in good_minima:
                ofile.write(" ".join([f"{j:10.4f}" for j in i])+"\n")
    except:
        print("No minima found. Skipping")
        pass
    plt.tight_layout()
    
    sns_cmap=plt.get_cmap(cmap)
    
    #ax=g.ax_marg_x
    #ax.plot(dataset_CV1['C1'],dataset_CV1['C2'],label="CV1")
    sns.set_style("ticks")
    sns.lineplot(x=dataset_CV1['C1'],y=dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$'],ax=g.ax_marg_x,color='black',linewidth=4)#,hue=dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$'],
    sns.scatterplot(x=dataset_CV1['C1'],y=dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$'],ax=g.ax_marg_x,hue=dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$'],
                palette=sns_cmap,hue_norm=(0,120),legend=False,s=120)
    # Create a set of line segments so that we can color them individually
    # This creates the points as an N x 1 x 2 array so that we can stack points
    # together easily to get the segments. The segments array for line collection
    # needs to be (numlines) x (points per line) x 2 (for x and y)
    #points = np.array([dataset_CV1['C1'], dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$']]).T.reshape(-1, 1, 2)
    #segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # Create a continuous norm to map from data points to colors
    #norm = plt.Normalize(0, 200,clip=True)
    #lc = mpl.collections.LineCollection(segments, cmap=cmap, norm=norm)
    # Set the values used for colormapping
    #lc.set_array(dataset_CV1['$\Delta$A ($kJ\ mol^{-1})$'])
    #lc.set_linewidth(5)
    #lc.set_label(laby)
    #line = g.ax_marg_x.add_collection(lc)
    #plt.xlabel("CV")
    #plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.xlim([0.5,4.0])
    limy=int(y_max_1d)
    g.ax_marg_x.set_ylim([-10,limy])
    g.ax_marg_x.set_xlabel(labx)
    g.ax_marg_x.set_yticks(range(0,limy+limy//2,limy//2))
    g.ax_marg_x.xaxis.set_ticks_position('none')
    g.ax_marg_x.set_frame_on(True)
    g.ax_marg_x.set_title(labx)
    #g.ax_marg_x.set_ylabel("Ciao")
    #plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    #plt.tight_layout()
    
    #ax=g.ax_marg_y
    #ax.plot(dataset_CV2['C1'],dataset_CV2['C2'],label="CV1")
    sns.lineplot(x=dataset_CV2['$\Delta$A ($kJ\ mol^{-1})$'],y=dataset_CV2['C1'],ax=g.ax_marg_y,orient="y",linewidth=4,color='black')
    sns.scatterplot(x=dataset_CV2['$\Delta$A ($kJ\ mol^{-1})$'],y=dataset_CV2['C1'],ax=g.ax_marg_y,hue=dataset_CV2['$\Delta$A ($kJ\ mol^{-1})$'],
                palette=sns_cmap,hue_norm=(0,120),legend=False,s=120)
    #plt.xlabel("CV")
    #plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    g.ax_marg_y.set_xlim([-10,limy])
    g.ax_marg_y.set_ylabel(laby,x=-0.5,y=-0.5)
    g.ax_marg_y.set_xticks(range(0,limy+limy//2,limy//2))
    g.ax_marg_y.yaxis.set_ticks_position('none')
    g.ax_marg_y.set_frame_on(True)
    g.ax_marg_y.set_title(laby)
    #g.ax_marg_y.set_xlabel("Free Energy ($kJ\ mol^{-1})$")
    #plt.ylim([-1.,40.0])
    #plt.legend()
    #plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    
    #g.ax_marg_x.annotate('X Marginal Label', xy=(0.5, 0), xytext=(0, -20),
    #                  xycoords='axes fraction', textcoords='offset points',
    #                  ha='center', va='top')
    #g.ax_marg_x.annotate('Y Marginal Label', xy=(0, 0.5), xytext=(-40, 0),
    #                  xycoords='axes fraction', textcoords='offset points',
    #                  ha='right', va='center', rotation=90)
    
    
    #plt.tight_layout()
    plt.box(on=True)
    plt.savefig('{}_joint.png'.format(file),format='png')

    #plt.show()
    
    #print(dataset)
    #sns.kdeplot(data=dataset,x='Column1', y='Column2', cmap="Reds", fill=True, bw_adjust=.5,clip=(0,200),levels=20)
    #sns.jointplot(data=dataset, x="Column1", y="Column2", hue="Column3", kind="kde")
    #sns.scatterplot(x='Column1', y='Column2', hue='Column3', data=dataset)
    #plt.show()
    
    
        
    #plt.savefig('{}.png'.format(file),format='png')


def main():

    cv1,cv2,free_grid,min_pt=rwg.extract_data(sys.argv[1])
    #print(cv1[min_pt[1][0]],cv2[min_pt[0][0]],free_grid[min_pt[0][0],min_pt[1][0]])
    file=os.path.splitext(sys.argv[1])[0]
    labx=sys.argv[2]
    laby=sys.argv[3]
    file_cv1=sys.argv[4]
    file_cv2=sys.argv[5]
    try:
        y_max_1d=sys.argv[6]
    except:
        y_max_1d=120
    try:
        minima=sys.argv[7]
    except:
        minima=""
    cmap_active='rainbow'
    print(y_max_1d)
    
    #saddle=allSaddles(free_grid)
    #print(saddle)
    #cmap_active=ListedColormap(np.linspace([0.16862745098, 0.219607843137,1,1],[1, 1,1,1],12)) #blue
    #cmap_active=ListedColormap(np.linspace([1.0, 0.40784313725490196,0,1],[1, 1,1,1],12)) #orange
    #cmap_active=ListedColormap(np.linspace([0.1960784313725490, 0.7686274509803922,0.4980392156862745,1],[1, 1,1,1],12)) #green
    0.1960784313725490
    
    #0.16862745098, 0.219607843137,1
    #cmap_active=LinearSegmentedColormap.from_list("mycmap",["#2b38ff","#FFFFFF"])#"#24BC99","#D6BB61","#E18F26","#FF6800"])
                                                  #oldcmap"#000000","#2b38ff","#17d9ff","#f7059b"])
                                                  
    #cmap_active='Blues_r'
    #path=minpath(free_grid,1,1,2,2)
    #min_path_cv1,min_path_cv2=get_minimum_path(cv1,cv2,free_grid)
    #reduced_fes_1,reduced_fes_2=dim_red(cv1,cv2,free_grid,300)
    #minima=detect_local_minima(free_grid)
    #maxima=detect_local_minima(-free_grid)
    #plot2d(cv1,cv2,free_grid,file,labx,laby,cmap_active,minima,min_pt)
    print(max(cv2),min(cv2))
    plot_joint(cv1,cv2,free_grid,file,labx,laby,cmap_active,minima,min_pt,file_cv1,file_cv2,y_max_1d)
    #plot3d(cv1,cv2,free_grid,file+'.dat',labx,laby)
    #plotminpath(min_path_cv1,min_path_cv2,file)
    #plot_reduced(cv1,cv2,reduced_fes_1,reduced_fes_2,file)
    #print(cv1)
 

if __name__ == "__main__":
    main()
