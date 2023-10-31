#%%
import MDAnalysis as mda
import numpy as np
import matplotlib as mpl
#matplotlib.use("pgf")
#matplotlib.rcParams.update({
#    "pgf.texsystem": "pdflatex",
#    'font.family': 'serif',
#    'text.usetex': True,
#    'pgf.rcfonts': False,
#})
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
from MDAnalysis.analysis import *
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.density import DensityAnalysis
from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis import distances
import glob
from itertools import combinations
import regex as re
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
from PIL import Image
from natsort import natsorted
import subprocess
#%%

def angle_with_axes(groupA,groupB,Universe):
        atom_A=Universe.select_atoms(groupA)
        atom_B=Universe.select_atoms(groupB)
        groupA_pos=[atom_A.positions for ts in Universe.trajectory]
        groupB_pos=[atom_B.positions for ts in Universe.trajectory]
        vec_pos_traj=[x-groupA_pos[i][0] for i,x in enumerate(groupB_pos)]
        dists=[np.linalg.norm(x) for x in vec_pos_traj]
        cos_angle=[[k[0][j]/dists[i] for j in range(3)]  for i,k in enumerate(vec_pos_traj)] 
        angle=[[np.arccos(x[j])*57.30 for j in range(3)] for x in cos_angle ]
        return(cos_angle,angle)

def set_size(width_pt, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)
def coord_Number(x,y,dist_min,dist_max,mols,dens):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*dens*np.trapz(prod,val_x)
    CN_one=CN_all/mols
    return(CN_one)
    
def RDF_COM_volume(gA,gB,Universe,binsize,max):
    groupA=Universe.select_atoms(gA)
    groupB=Universe.select_atoms(gB)
    gA_coms=[groupA.center_of_mass(unwrap=True,compound='residues') for ts in Universe.trajectory]
    gB_coms=[groupB.center_of_mass(unwrap=True,compound='residues') for ts in Universe.trajectory]
    boxes=[Universe.dimensions for ts in Universe.trajectory]
    if groupA == groupB:
        com_arr  = [distances.self_distance_array(x,box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist = [np.tile(x,2) for x in com_arr]
    else:
        com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist=[x.reshape(x.shape[0]*x.shape[1]) for x in com_arr]
    bin_size = binsize
    # Create an array for the radial bins
    bins = np.arange(0, max, bin_size)
    hists=[np.histogram(x, bins=bins) for x in com_dist]
    vols=[u.dimensions[0]*u.dimensions[1]*u.dimensions[2] for ts in u.trajectory]
    avg_vol=np.mean(vols)
    dens=len(com_dist[0])/avg_vol
    sum_tmp=np.zeros(len(bins)-1)
    for i in hists:
        sum_tmp=sum_tmp+i[0] 
    sum_norm=(sum_tmp/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)*(len(gA_coms)-1)))
    
    return(bins,sum_tmp,sum_norm,avg_vol,dens)


def RDF_POS_volume(gA,gB,Universe,binsize,max):
    groupA=Universe.select_atoms(gA)
    groupB=Universe.select_atoms(gB)
    gA_coms=[groupA.positions for ts in Universe.trajectory]
    gB_coms=[groupB.positions for ts in Universe.trajectory]
    boxes=[Universe.dimensions for ts in Universe.trajectory]
    if groupA == groupB:
        com_arr  = [distances.self_distance_array(x,box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist = [np.tile(x,2) for x in com_arr]
    else:
        com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist=[x.reshape(x.shape[0]*x.shape[1]) for x in com_arr]
    bin_size = binsize
    # Create an array for the radial bins
    bins = np.arange(0, max, bin_size)
    hists=[np.histogram(x, bins=bins) for x in com_dist]
    vols=[u.dimensions[0]*u.dimensions[1]*u.dimensions[2] for ts in u.trajectory]
    avg_vol=np.mean(vols)
    dens=len(groupB)/avg_vol
    sum_tmp=np.zeros(len(bins)-1)
    for i in hists:
        sum_tmp=sum_tmp+i[0] 
    sum_norm=(sum_tmp/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)*(len(gA_coms))))
    
    return(bins,sum_tmp,sum_norm,avg_vol)
        
def RDF_density(gA,gB,Universe,binsize,max,density):
    groupA=Universe.select_atoms(gA)
    groupB=Universe.select_atoms(gB)
    gA_coms=[groupA.positions for ts in Universe.trajectory]
    gB_coms=[groupB.positions for ts in Universe.trajectory]
    boxes=[Universe.dimensions for ts in Universe.trajectory]
    if groupA == groupB:
        com_arr  = [distances.self_distance_array(x,box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist = [np.tile(x,2) for x in com_arr]
    else:
        com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
        com_dist=[x.reshape(x.shape[0]*x.shape[1]) for x in com_arr]
    bin_size = binsize
    # Create an array for the radial bins
    bins = np.arange(0, max, bin_size)
    hists=[np.histogram(x, bins=bins) for x in com_dist]
    vols=[u.dimensions[0]*u.dimensions[1]*u.dimensions[2] for ts in u.trajectory]
    avg_vol=np.mean(vols)
    dens=density*len(groupB)
    sum_tmp=np.zeros(len(bins)-1)
    for i in hists:
        sum_tmp=sum_tmp+i[0] 
    sum_norm=(sum_tmp/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)*(len(gA_coms))))
    
    return(bins,sum_tmp,sum_norm,avg_vol)

def RDF_gmx_number_density(gA,gB,Universe,binsize,max):
    groupA=Universe.select_atoms(gA)
    groupB=Universe.select_atoms(gB)
    gA_coms=[groupA.positions for ts in Universe.trajectory]
    gB_coms=[groupB.positions for ts in Universe.trajectory]
    boxes=[Universe.dimensions for ts in Universe.trajectory]
    if groupA == groupB:
        com_arr  = [distances.self_distance_array(x,box=boxes[i]) for i,x in enumerate(gA_coms)]
        com_dist = [np.tile(x,2) for x in com_arr]
        dens=1
    else:
        com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i]) for i,x in enumerate(gA_coms)]
        com_dist=[x.reshape(x.shape[0]*x.shape[1]) for x in com_arr]
        dens=1
    bin_size = binsize
    # Create an array for the radial bins
    bins = np.arange(0, max, bin_size)
    hists=[np.histogram(x, bins=bins) for x in com_dist]
    vols=[u.dimensions[0]*u.dimensions[1]*u.dimensions[2] for ts in u.trajectory]
    avg_vol=np.mean(vols)
    sum_tmp=np.zeros(len(bins)-1)
    for i in hists:
        sum_tmp=sum_tmp+i[0]/len(groupA)
    sum_tmp_avg=sum_tmp/len(gA_coms)
    sum_norm=(sum_tmp_avg/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)))
    
    return(bins,sum_tmp,sum_norm,avg_vol)

def plot_RDF_MDanal(RDF,r,nmols):
    plt.plot(RDF.results.bins,RDF.results.rdf)
    cn=coord_Number(RDF.results.bins,RDF.results.rdf,0,r,nmols)
    CN = 4*np.pi*np.trapz(RDF.rdf[RDF.bins <=r]*RDF.bins[RDF.bins <=r]**2, RDF.bins[RDF.bins <=r])
    print(CN/nmols,cn)

def plot_hbonds(folder:str,textsize:int,fraction:float,labels:list,title:str=""):
    HB={}
    hbonds=glob.glob(os.path.join(folder,"*hbond*"))
    for x in hbonds:
        frames=[]
        hb=[]
        with open(x,'r') as ifile:
            lines=ifile.readlines()
        title=lines[0].split()[1].split('[')[0]
        for j in lines:
            try:
                frames.append(int(j.split()[0]))
                hb.append(int(j.split()[1].rstrip()))
            except:
                continue
        HB[title]={'frame':frames,'HB':hb}
    HB_tot=sum([x for k in HB for x in HB[k]['HB'] ])/len(frames)
    wd,hg=set_size(textsize,fraction)
    lbls={x:labels[i] for i,x in enumerate(HB)}
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    clrs=['blue','cyan','tan','violet']+[x for x in mcolors.CSS4_COLORS]
    for i,x in enumerate(HB):
        plt.plot(HB[x]['frame'],[j/HB_tot for j in HB[x]['HB']],label=lbls[x],color=clrs[i])
        print("Average HB percentage of {} = {}".format(x,np.mean([j/HB_tot for j in HB[x]['HB']])))
    plt.legend()
    plt.ylim([0,1])
    plt.xlabel("Frame")
    plt.ylabel("H-Bond fraction")
    plt.savefig(os.path.join(folder,'HB.pgf'),format='pgf',bbox_inches = "tight")
#%%    
def plot_COM_RDF(root:str,textsize:int,fraction:float,labels:list,title:str,resnumb:dict):
    with open(os.path.join(root,'cpptraj.out')) as ifile:
        lines=ifile.readlines()
        names=[]
        Dd=[]
        numb=[]
        for i,line in enumerate(lines):
            if "Calculating RDF" in line:
                names.append(line.split()[7]+line.split()[-1])
            if "in common" in line:
                numb.append(max(int(line.split()[6].rstrip(',')),int(line.split()[11].rstrip(','))))
            if "Average density" in line:
                Dd.append(float(line.split()[3]))
    print(Dd)        
    D2={}

    id2="_rdf.dat"
    coms2=glob.glob(os.path.join(root,"*com*"+id2))
    RDFs2={}
    names2=[]
    ords2=[]
    print(coms2)
    for i in coms2:
        with open(i,'r') as ifile:
            name=os.path.basename(i).replace(id2,"")
            names2.append(name)
            lines=ifile.readlines()
            res=name.split('_')[2]
            ords2.append(name.split('_')[0])
            x=[]
            y=[]
            y_min=[]
            for line in lines[1:]:
                    x.append(float(line.split()[0]))
                    y.append(float(line.split()[1]))
            y_max_index=find_peaks(y,width=1)
            x_max=x[y_max_index[0][0]]
            y_min=find_peaks([-x for x in y],width=2)
            x_min=x[y_min[0][0]]
        RDFs2[name.split('_')[0]]={'name':'_'.join(name.split('_')[1:]),'RDF':[x,y],'r_max':x_max,'r_min':x_min,'minima':y_min,'maxima':y_max_index,'#_part':resnumb[res]}
    for i,x in enumerate(ords2):
        D2[x]={'name':names[i],'density':round(Dd[i],6)}

    for i in RDFs2:
        RDFs2[i]['CN']=coord_Number(RDFs2[i]['RDF'][0],RDFs2[i]['RDF'][1],0,RDFs2[i]['r_min'],
                                    RDFs2[i]['#_part'],D2[i]['density'])
    #print(RDFs['GCL_C2_GCL_C2'])

    with open (os.path.join(root,'results_COM.dat'),'w') as ofile:
        ofile.write("{:35s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN'))
        for i in RDFs2:
            ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f}\n".format(RDFs2[i]['name'],
                        RDFs2[i]['r_max'],RDFs2[i]['r_min'],RDFs2[i]['CN']))
    wd,hg=set_size(textsize,fraction)
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    clrs=['blue','orange','green','red','purple','black']+[x for x in mcolors.TABLEAU_COLORS]
    legd=labels
    for i,item in enumerate(RDFs2):
        plt.plot(RDFs2[item]['RDF'][0], RDFs2[item]['RDF'][1],color=clrs[i],label=item)
    plt.xlabel("Dist / \AA")
    plt.ylabel("$g(r)$")
    plt.legend(legd)    
    for i,item in enumerate(RDFs2):
        plt.scatter(RDFs2[item]['r_max'],RDFs2[item]['RDF'][1][RDFs2[item]['maxima'][0][0]],color=clrs[i])
        plt.scatter(RDFs2[item]['r_min'],RDFs2[item]['RDF'][1][RDFs2[item]['minima'][0][0]],color=clrs[i])
    tit=title
    plt.title(tit)
    #plt.savefig(os.path.join(root,tit+'.pgf'), format='pgf',bbox_inches = "tight")
    plt.savefig(os.path.join(root,tit+'.png'), format='png',bbox_inches = "tight")
#%%
def plot_DENS_cpptraj(root:str,textsize:int,fraction:float,labels:list,title:str="",ymax:float=0,shade:dict={}):
    files=glob.glob(os.path.join(root,"*DENS*"))
    dens={}
    for file in files:
        with open(file,'r') as ifile:
            lines=ifile.readlines()
            vars=len(lines[0].split())
            dens[ifile.name]={lines[0].split()[0]:[float(x.split()[0]) for x in lines[1:]]}
            for i in range(1,vars,2):
                vardata={lines[0].split()[i]:[float(x.split()[i]) for x in lines[1:]]}
                dens[ifile.name].update(vardata)
    wd,hg=set_size(textsize,fraction)
    for file in dens:
        max_dens=0
        fig=plt.figure(figsize=(wd,hg),dpi=150)
        vars=list(dens[file].keys())
        for var in range(1,len(dens[file])):
            if max(dens[file][vars[var]]) > max_dens:
                max_dens=max(dens[file][vars[var]])
            if ymax < max(dens[file][vars[var]]):
                ymax = max(dens[file][vars[var]])*0.1+max(dens[file][vars[var]])
            plt.plot(dens[file][vars[0]],dens[file][vars[var]],label=''.join(e for e in vars[var] if e.isalnum()))
        plt.ylim((0,ymax))
        plt.xlabel(" Z / \AA")
        plt.ylabel("$ Density / Kg m^{-3}$ ")
        plt.legend()
        for i in shade:
            plt.axvspan(shade[i]["min"],shade[i]["max"], color='#989898', alpha=0.5, lw=0)
        plt.title(title)
        plt.savefig(os.path.join(root,title+'.pgf'), format='pgf',bbox_inches = "tight")
        plt.close()
    
        
    return(dens,vars,fig)
def plot_DENS_gmx(root:str,files:dict,textsize:int,fraction:float,labels:list,title:str="",ymax:float=0,shade:dict={},phase:dict={}):
    dens={}
    arr_starts=[]
    arr_ends=[]
    for file in files:
        with open(os.path.join(root,files[file]),'r') as ifile:
            lines=ifile.readlines()
            dens[file]={"rho":[float(x.split()[1]) for x in lines if "@" not in x if "#" not in x],
                        "Z":[float(x.split()[0])*10 for x in lines if "@" not in x if "#" not in x]}
    wd,hg=set_size(textsize,fraction)
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    for file in dens:
        max_dens=0
        if max(dens[file]["rho"]) > max_dens:
            max_dens=max(dens[file]["rho"])
        if ymax < max(dens[file]["rho"]):
            ymax = max(dens[file]["rho"])*0.1+max(dens[file]["rho"])
        plt.plot(dens[file]["Z"],dens[file]["rho"],label=file)
        plt.ylim((0,ymax))
        plt.xlabel(" Z / \AA")
        plt.ylabel("$ Density / Kg m^{-3}$ ")
        plt.legend()
        for i in shade:
            plt.axvspan(shade[i]["min"],shade[i]["max"], color='#989898', alpha=0.5, lw=0)
        for i in phase:
            plt.annotate(text='', xy=(phase[i]['start'],phase[i]['arry']), xytext=(phase[i]['end'],phase[i]['arry']), arrowprops=dict(arrowstyle='<->')) 
            plt.text((phase[i]['start']+phase[i]['end'])/2,phase[i]['arry']+phase[i]['arry']*0.05,phase[i]['name'],ha='center')

        plt.title(title)
    plt.savefig(os.path.join(root,title+'.pgf'), format='pgf',bbox_inches = "tight")
    plt.close()
    
        
    return(dens,fig)

def extract_umbrella_profile(root:str,profile:str,histo:str):
    data=[]
    data_hist=[]
    chars=["#","@"]
    names=[]
    rho=[]
    with open(os.path.join(root,profile),'r') as ifile:
            lines=ifile.readlines()
            x=[float(x.split()[0])*10 for x in lines if x[0] not in chars ]
            y=[float(x.split()[1]) for x in lines if x[0] not in chars ]
            x_corr=[j+57.6 for j in x]
            y_scaled=[i-min(y) for i in y]
            data.append([x_corr,y_scaled])
    with open(os.path.join(root,histo),'r') as ifile:
            lines=[x for x in ifile.readlines() if "#" not in x if "@" not in x]
            data_hist.append([float(x.split()[0])*10+57.6 for x in lines])
            for i in range(1,len(lines[0].split())):
                data_hist.append([float(x.split()[i]) for x in lines])
    return(data,data_hist)

def extract_umbrella_profile_bootstrap(root:str,profile:str):
    data=np.loadtxt(os.path.join(root,profile),unpack=True,comments=["#","@"])
    return(data)
            
    

def plot_image(file:str,textsize:int,fraction:float):
    dir=os.path.dirname(file)
    wd,hg=set_size(textsize,fraction)
    img = np.asarray(Image.open(file))
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    plt.imshow(img)
    plt.axis('off')
    plt.savefig(file+'.pgf')
def calculate_orientation(root,top,sel1,sel2):
    us=[]
    paths_sorted=natsorted(glob.glob(os.path.join(root,"window*")))
    for path in paths_sorted:
            if os.path.exists(os.path.join(path,'ORIENTATION_OK')):
                    print("{} already done. Skipping...".format(path))
            else:
                    files=glob.glob(os.path.join(path,'*.xtc'))
                    dirs=[os.path.dirname(x) for x in files] 
                    fnames=[os.path.basename(x).replace('.xtc','') for x in files]
                    dumps=[os.path.basename(x).replace('_whole.xtc','.dump') for x in files]
                    us=[mda.Universe(top,x) for x in files]
                    orientation=[angle_with_axes(sel1,sel2,x) for x in us]
                    ang_avg=[]
                    std_ang_avg=[]
                    cos_avg=[]
                    std_cos_avg=[]
                    for index,i in enumerate(orientation):
                            angx=[]
                            angy=[]
                            angz=[]
                            cosx=[]
                            cosy=[]
                            cosz=[]
                            for cos in i[0]:
                                    cosx.append(cos[0])    
                                    cosy.append(cos[1])
                                    cosz.append(cos[2])
                            for ang in i[1]:
                                    angx.append(ang[0])    
                                    angy.append(ang[1])
                                    angz.append(ang[2])
                            ang_avg.append([np.mean(angx),np.mean(angy),np.mean(angz)])
                            std_ang_avg.append([np.std(angx),np.std(angy),np.std(angz)])
                            cos_avg.append([np.mean(cosx),np.mean(cosy),np.mean(cosz)])
                            std_cos_avg.append([np.std(cosx),np.std(cosy),np.std(cosz)])
                            with open(os.path.join(dirs[index],dumps[index]),'r') as ifile:
                                            lines=ifile.readlines()
                                            posu=[float(x.split()[-1]) for x in lines if " init " in x]
                            with open(os.path.join(dirs[index],fnames[index]+'_orientation.dat'),'w') as ofile:
                                    ofile.write("{:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n".format("cosX","cosY","cosZ","angX","angY","angZ"))
                                    for idx,x in enumerate(cosx):
                                            ofile.write("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(x,cosy[idx],cosz[idx],angx[idx],angy[idx],angz[idx]))
                                    ofile.write("{:8s} {:8.3f}\n".format("##POS UMBRELLA##",posu[0]))
                                    ofile.write("{:20s} {:8.3f} {:8.3f} {:8.3f} \n".format(
                                                    "##COS MEANS##", cos_avg[index][0],cos_avg[index][1],cos_avg[index][2]))
                                    ofile.write("{:20s} {:8.3f} {:8.3f} {:8.3f} \n".format(
                                                    "##COS STD##", std_cos_avg[index][0],std_cos_avg[index][1],std_cos_avg[index][2]))
                                    ofile.write("{:20s} {:8.3f} {:8.3f} {:8.3f} \n".format(
                                                    "##ANG MEANS##", ang_avg[index][0],ang_avg[index][1],ang_avg[index][2]))
                                    ofile.write("{:20s} {:8.3f} {:8.3f} {:8.3f} \n".format(
                                                    "##ANG STD##", std_ang_avg[index][0],std_ang_avg[index][1],std_ang_avg[index][2]))
                    subprocess.call(['touch',os.path.join(path,'ORIENTATION_OK')])
                    print("FINISHED {}".format(path))
        

def select_orientation(inputfile:str,outfile:str):
    root=os.path.dirname(inputfile)
    ang_avg_all=[]
    std_ang_all=[]
    pos_umbrella=[]
    dens_CHL=[]
    dens_GCL=[]
    dens_Clm=[]
    with open(inputfile,'r') as ifile:
            lines=ifile.readlines()
            fnames=[x.replace('.tpr','_whole_orientation.dat').rstrip() for x in lines]
    ang_mean_wd={}
    ang_std_wd={}
    nframes=[]
    for file in fnames:
            with open(os.path.join(root,file),'r') as ifile:
                    lines=ifile.readlines()
                    nframes=[(len(lines)-6)]
                    pos=[float(line.split()[2]) for line in lines if "##POS UMBRELLA##" in line]
                    #print(pos)
                    means=[[float(line.split()[2]),float(line.split()[3]),float(line.split()[4])] for line in lines if "##ANG MEANS##" in line]
                    if pos[0] in ang_mean_wd:
                            ang_mean_wd[pos[0]]['mean'].append(means[0])
                            ang_mean_wd[pos[0]]['frames'].append(nframes[0])
                            
                    else:
                            ang_mean_wd[pos[0]]={'mean':means}
                            ang_mean_wd[pos[0]]['frames']=nframes
                            
                    std_means=[[float(line.split()[2]),float(line.split()[3]),float(line.split()[4])] for line in lines if "##ANG STD##" in line ]
                    if pos[0] in ang_std_wd:
                            ang_std_wd[pos[0]].append(std_means[0])
                    else:
                            ang_std_wd[pos[0]]=std_means
    for k in ang_mean_wd:
            pos_umbrella.append(k)
            
            ang_avg_all.append([np.average([x[0] for x in ang_mean_wd[k]['mean']],weights=ang_mean_wd[k]['frames']),
                                    np.average([x[1] for x in ang_mean_wd[k]['mean']],weights=ang_mean_wd[k]['frames']),
                                    np.average([x[2] for x in ang_mean_wd[k]['mean']],weights=ang_mean_wd[k]['frames'])])
            std_ang_all.append([np.sqrt(sum([x[0]**2 for x in ang_std_wd[k] ])),np.sqrt(sum([x[1]**2 for x in ang_std_wd[k] ])),np.sqrt(sum([x[2]**2 for x in ang_std_wd[k] ]))])
    tuls=[(pos_umbrella[i],ang_avg_all[i],std_ang_all[i]) for i in range(len(pos_umbrella))]
    tuls.sort()
    with open(outfile,'w') as ofile:
            ofile.write("{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}\n".format("Pos","angX","stdX","angY","stdY","angZ","stdZ"))
            for x in tuls:
                    ofile.write("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(x[0],x[1][0],x[2][0],
                                                                                                x[1][1],x[2][1],
                                                                                                x[1][2],x[2][2]))
    return(tuls)
#%%
def extract_RDF_gmx(root:str,mols:list):                                                                                     
    paths_sorted=natsorted(glob.glob(os.path.join(root,'window*')))

    RDFS={x:{} for x in mols}
    COORD={x:{} for x in mols}

    for traj in paths_sorted:
            dir=os.path.dirname(traj)
            window=os.path.basename(traj)
            pairs=[]
            for j,mol in enumerate(RDFS):
                RDFS[mol][window]={}
                COORD[mol][window]={}
                pairs.append([])
                #RDFS
                with open(os.path.join(traj,'RDF_'+mol+'.xvg'),'r') as ifile:
                        lines=ifile.readlines()
                        for line in lines:
                                if re.findall('@ s\d',line) :
                                        pairs[j].append(line.split()[-1])
                        x=[float(x.split()[0])*10 for x in lines if "#" not in x if "@" not in x]
                        RDFS[mol]['dist']=x
                        for i,t in enumerate(pairs[j]):
                                RDFS[mol][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
                #COORDINATION NUMBER                
                with open(os.path.join(traj,'coord_'+mol+'.xvg'),'r') as ifile:
                        lines=ifile.readlines()
                        x=[float(x.split()[0])*10 for x in lines if "#" not in x if "@" not in x]
                        COORD[mol]['dist']=x
                        for i,t in enumerate(pairs[j]):
                                COORD[mol][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
    return(mols,pairs,RDFS,COORD)


def extract_coordination_Number(tuls:list,RDFS:list):
    posz_all=[x[0]*10+57.6 for x in tuls]
    axis=["X","Y","Z"]
    CN={x:{} for x in RDFS}
    data={x:{} for x in RDFS}
    for k in RDFS:
        for numb,i in enumerate(RDFS[k]):
            idx=0
            if i == 'dist':
                idx=numb
            else:
                if numb > idx:
                    CN[k][i]={'pos':posz_all[numb-1]}
                    data[k][i]={'pos':posz_all[numb-1]}
                else:
                    CN[k][i]={'pos':posz_all[numb]}
                    data[k][i]={'pos':posz_all[numb]}
                for j in RDFS[k][i]:
                        try:
                                y_min=find_peaks([-x for x in RDFS[k][i][j]],width=1)
                                y_max=find_peaks([x for x in RDFS[k][i][j]],width=1)
                                x_min=RDFS[k]['dist'][y_min[0][0]]# if RDFS[k]['dist'][y_min[0][0]] < 4 else 0]
                                x_min_data=[RDFS[k]['dist'][y_min[0][0] if RDFS[k]['dist'][y_min[0][0]] < 3.7 else 0]]
                                x_max_data=[RDFS[k]['dist'][y_max[0][0] if RDFS[k]['dist'][y_max[0][0]] < 3.7 else 0]]
                                data[k][i][j]={"val":COORD[k][i][j][RDFS[k]['dist'].index(x_min)],"rmax":x_max_data[0],"rmin":x_min_data[0]}
                                CN[k][i][j]=COORD[k][i][j][RDFS[k]['dist'].index(x_min)]
                        except:
                                CN[k][i][j]=0
                                data[k][i][j]={"val":0,"rmax":0,"rmin":0}
    return(CN,data)
#%%
def calc_avg_distance(universe,traj,sel1,sel2,max):
    temp=[]
    Universe = mda.Universe(universe,traj)
    groupA=Universe.select_atoms(sel1)
    groupB=Universe.select_atoms(sel2)
    dist=[mda.analysis.distances.distance_array(groupA,groupB) for ts in Universe.trajectory]
    avg_dist=[x[(0<=x)&(x<=7.0)].mean() for x in dist if len(x[(0<=x)&(x<=5.0)])>=1 ]
    return(dist,avg_dist)
    

un="/run/media/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/ANNEALED_BOX/UMBRELLA_BIG_IMC_MOD_ANNEALED.parm7"
trj="/run/media/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/ANNEALED_BOX/ANNEALED_BOX/window3.9/umbrella_whole.xtc"
dst,av_dst=calc_avg_distance(un,trj,'type mg','type Clm',3.5)
#%%
print(dst[0][(dst[0]<=10)])
#print(dst[0][(0<=dst[0])&(dst[0]<=5.0)])
#print(dst[0][(0<=dst[0])&(dst[0]<=5.0)].mean())
#print(av_dist)
                      
                      
                      
                      
                      
                                      
#%%
############################################
### ANALYSIS
############################################


root_HB="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/AcPh/ANNEALING/1.0mmol_DESonly"
plot_hbonds(root_HB,345,1,["Choline-Cl","Gly-Gly","Gly-Cl","choline-Gly","Gly-AcPh","Chl-AcPh"])

# %%
root_COM="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/AcPh/ANNEALING/0.2mmol_DESonly"
#root_COM="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/DES/ANNEALED_BOX"
resnumb={"CHL":400,"Clm":400,"GCL":800,"ACP":26}
#resnumb={"CHL":400,"Clm":400,"GCL":800}
plot_COM_RDF(root_COM,500,1,["Choline-Choline","Choline-chloride","Choline-glycerol","Glycerol-chloride","Glycerol-glycerol",'Chloride-Chloride',
                            'AcPh-Choline','AcPh-Chloride','AcPh-glycerol','AcPh-AcPh'],"0.90",resnumb)
#plot_COM_RDF(root_COM,500,1,["Choline-Choline","Choline-chloride","Choline-glycerol","Glycerol-chloride","Glycerol-glycerol",'Chloride-Chloride'
#                            ],"0.90",resnumb)

# %%
root_DENS_IMC="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/ANNEALED"
root_DENS_ACP="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/AcPh/ANNEALING"
#D,v,fig=plot_DENS_cpptraj(root_DENS,345,1,["Acetophenone"],"rho-ACP",150,{1:{"min":94,"max":107},2:{"min":204,"max":206},3:{"min":-1,"max":11.7}})
D_IMC,fig_IMC=plot_DENS_gmx(root_DENS_IMC,{"ACT":"ACT_density.xvg","IMC":"IMC_density.xvg"},345,1,["Acetophenone"],"rho-ACT-IMC",60,
                            {1:{"min":94,"max":107},2:{"min":204,"max":206},3:{"min":-1,"max":11.7}},
                            {1:{"name":"DES","start":11.7,"end":94,"arry":35},2:{"name":"THF","start":107,"end":204,"arry":35}})
D_ACP,fig_ACP=plot_DENS_gmx(root_DENS_ACP,{"ACP":"ACP_density.xvg"},345,1,["Acetophenone"],"rho-ACP",60,
                            {1:{"min":94,"max":107},2:{"min":204,"max":206},3:{"min":-1,"max":11.7}},
                            {1:{"name":"DES","start":11.7,"end":94,"arry":35},2:{"name":"THF","start":107,"end":204,"arry":35}})
# %% PLOT PMFS
root_PMF="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
#mols=["IMC","ACT","ACP"]
font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 32}
mpl.rc('font', **font)
mols=["ACP"]
dir="umbrella_30_200/ANNEALED_BOX/ANNEALED_BOX"
textsize=800
fraction=1.5
a=10.1
b=11.5
limx=[9.5,12.0]
limy=[-1,12]
arry=10
bootstrap=True
graph_data=[extract_umbrella_profile(os.path.join(root_PMF,x,dir),'profile.xvg','histo.xvg') for x in mols]
if bootstrap:
    bs=[extract_umbrella_profile_bootstrap(os.path.join(root_PMF,x,dir),'bootstrap_avg.xvg') for x in mols]
wd,hg=set_size(textsize,fraction)
#wd,hg=15,10

fig=plt.figure(figsize=(wd,hg),dpi=150)
for i,x in enumerate(graph_data):
    plt.subplot(len(graph_data),1,i+1)
    #plt.plot([k/10 for k in x[0][0][0]],x[0][0][1],label=mols[i])
    if bootstrap:
        print(min(bs[i][1]))
        plt.plot([k+5.7 for k in bs[i][0]],[k-min(bs[i][1]) for k in bs[i][1]] ,label=mols[i])
        plt.fill_between([x+5.7 for x in bs[i][0]],[k-min(bs[i][1]) for k in bs[i][1]]+bs[i][2],
                         [k-min(bs[i][1]) for k in bs[i][1]]-bs[i][2],alpha=0.3)
        
    plt.xlim(limx)
    plt.ylim(limy)
    #plt.legend()
    #plt.annotate(text='', xy=(limx[0],arry), xytext=(a,arry), arrowprops=dict(arrowstyle='<->')) 
    plt.text((limx[0]+a)/2,arry+arry*0.05,"DES",ha='center')
    #plt.annotate(text='', xy=(b,arry), xytext=(limx[1],arry), arrowprops=dict(arrowstyle='<->'))
    plt.text((b+limx[1])/2,arry+arry*0.05,"THF",ha='center')
    plt.axvspan(a, b, color='#989898', alpha=0.5, lw=0)
    plt.ylabel("Free Energy $(Kj\ mol^{-1})$")
    plt.tight_layout()
plt.xlabel("Z coordinate (nm)")


plt.savefig(os.path.join(root_PMF,'_'.join(mols)+'_PMF.png'), format='png',bbox_inches = "tight")

#PLOT HISTOS
fig=plt.figure(figsize=(wd,hg*2),dpi=150)
for i,x in enumerate(graph_data):
    plt.subplot(3,1,i+1)
    for j in range(1,len(x[1])):
        plt.plot(x[1][0],x[1][j])
    plt.xlim(limx)
    #plt.legend(mols[i])
    #plt.ylim(limy)
    #plt.legend()
    plt.axvspan(a, b, color='#989898', alpha=0.5, lw=0)
    plt.ylabel("Free Energy / $Kjmol^{-1}$")
    plt.tight_layout
plt.xlabel("Z / $\AA$")


plt.savefig(os.path.join(root_PMF,'_'.join(mols)+'_HIST.pgf'), format='pgf',bbox_inches = "tight")



# %%
plot_image('/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/img/interface_struct.pov.png',426,1)




#%%
# calc ORIENTATION
ORIENT_root="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/ACP/umbrella_30_200/ANNEALED_BOX/ANNEALED_BOX"
TOP="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/ACP/umbrella_30_200/ANNEALED_BOX/UMBRELLA_BIG_ACP_MOD_ANNEALED.parm7"
calculate_orientation(ORIENT_root,TOP,"resname ACP and type c","resname ACP and type o")

#%% extract ORIENTATION

tuls=select_orientation(os.path.join(ORIENT_root,'tprfiles.dat'),os.path.join(ORIENT_root,'orientation_results.dat'))
mols,pairs,RDFS,COORD=extract_RDF_gmx(ORIENT_root,["O"])
orientation_PROFILE=extract_umbrella_profile(ORIENT_root,'profile.xvg','histo.xvg')
CN,dst=extract_coordination_Number(tuls,RDFS)
# %% PLOT ORIENTATION

wd,hg=set_size(426,1)
fig=plt.figure(figsize=(wd,hg*2),dpi=150)

root_ORIENT='/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/ACP/umbrella_30_200/ANNEALED_BOX/ANNEALED_BOX'
limx_all=[95,140]
limy_all=[0,40]
int_start=105
int_end=115
arry_orient=2.5
limx=[95,140]
limy=[-1,15]
a=105
b=115
ltr=['a','b','c','d','e']
#labels=[["THF-O","CHL-O","GCL-O","Clm"],["CHL-HO","CHL-N","GCL-HO"]]
labels=[["CHL-HO","GCL-HO","CHL-N"]]
ax=plt.subplot(4,1,1)
plt.errorbar([x[0]*10+57.6 for x in tuls],[x[1][2] for x in tuls],[x[2][2] for x in tuls],linestyle=None)
plt.plot([x[0]*10+57.6 for x in tuls],[x[1][2] for x in tuls])
plt.text(.9,.9,'a',transform=ax.transAxes)
plt.axvline(x=int_start, color = '#989898')
plt.axvline(x=int_end, color = '#989898')
plt.axhline(y=90, color = '#989898')
plt.xlim(limx_all)
plt.ylabel("Z Angle / deg")
for i,mol in enumerate(mols):
    ax=plt.subplot(4,1,2+i)
    plt.text(.9,.9,ltr[1+i],transform=ax.transAxes)
    for k,pair in enumerate(pairs[i]):
            #plt.plot([CN[mol][x]['pos'] for x in CN[mol]],[CN[mol][x][pair]['val'] for x in CN[mol] if x != 'dist'],'o-',label=labels[i][k])
            plt.plot([CN[mol][x]['pos'] for x in CN[mol]],[CN[mol][x][pair] for x in CN[mol] if x != 'dist'],'o-',label=labels[i][k])
            plt.ylabel("Coord Number")
            plt.axvline(x=int_start, color = '#989898')
            plt.axvline(x=int_end, color = '#989898')
            plt.legend(loc='lower right')
            plt.xlim(limx_all)
#plt.plot([CN['Mg'][x]['pos'] for x in CN['Mg']],[CN['Mg'][x][pairsMg[0]] + CN['Mg'][x][pairsMg[1]]+CN['Mg'][x][pairsMg[2]]+CN['Mg'][x][pairsMg[3]] for x in CN['Mg']],'-o',label="Sum")
            plt.ylim([-0.5,3.5])
plt.annotate(text='', xy=(limx[0],arry_orient), xytext=(a,arry_orient), arrowprops=dict(arrowstyle='<->')) 
plt.text((limx[0]+a)/2,arry_orient+arry_orient*0.05,"DES",ha='center')
plt.annotate(text='', xy=(b,arry_orient), xytext=(limx[1],arry_orient), arrowprops=dict(arrowstyle='<->'))
plt.text((b+limx[1])/2,arry_orient+arry_orient*0.05,"THF",ha='center')


ax=plt.subplot(4,1,2+len(mols))
plt.text(.9,.9,ltr[1+len(mols)],transform=ax.transAxes)
plt.plot(orientation_PROFILE[0][0][0],orientation_PROFILE[0][0][1])
plt.axvline(x=int_start, color = '#989898')
plt.axvline(x=int_end, color = '#989898')
plt.ylabel("Free Energy / Kj/mol")
plt.xlabel("Z / \AA")
plt.xlim(limx_all)
plt.tight_layout()
plt.savefig(os.path.join(root_ORIENT,'IMC_ORIENT.png'), format='png',bbox_inches = "tight")
plt.show()
plt.close()
# %%
# %% PLOT RDFS MAXIMA and MINIMA
arry_orient=1
limy_dst=[0,3.7]
wd,hg=set_size(426,1)
fig=plt.figure(figsize=(wd,hg),dpi=150)
#ax=plt.subplot(2,1,1)
plt.plot([dst['Mg'][x]['pos'] for x in dst['Mg']] ,[dst['Mg'][x]['"Clm"']['rmax'] for x in dst['Mg'] if x != 'dist'],'o-',label="rmax")
plt.plot([dst['Mg'][x]['pos'] for x in dst['Mg']] ,[dst['Mg'][x]['"Clm"']['rmin'] for x in dst['Mg'] if x != 'dist'],'o-',label="rmin")
plt.annotate(text='', xy=(limx[0],arry_orient), xytext=(a,arry_orient), arrowprops=dict(arrowstyle='<->')) 
plt.text((limx[0]+a)/2,arry_orient+arry_orient*0.05,"DES",ha='center')
plt.annotate(text='', xy=(b,arry_orient), xytext=(limx[1],arry_orient), arrowprops=dict(arrowstyle='<->'))
plt.text((b+limx[1])/2,arry_orient+arry_orient*0.05,"THF",ha='center')
plt.legend()
plt.axvline(x=int_start, color = '#989898')
plt.axvline(x=int_end, color = '#989898')
plt.ylabel("R / \AA")
plt.xlabel("Z / \AA")
plt.xlim(limx_all)
plt.ylim(limy_dst)
"""
ax=plt.subplot(2,1,2)
plt.plot([dst['Mg'][x]['pos'] for x in dst['Mg']] ,[dst['Mg'][x]['"Clm"']['rmin'] for x in dst['Mg'] if x != 'dist'],'o-',label="rmin")
plt.annotate(text='', xy=(limx[0],arry_orient), xytext=(a,arry_orient), arrowprops=dict(arrowstyle='<->')) 
plt.text((limx[0]+a)/2,arry_orient+arry_orient*0.05,"DES",ha='center')
plt.annotate(text='', xy=(b,arry_orient), xytext=(limx[1],arry_orient), arrowprops=dict(arrowstyle='<->'))
plt.text((b+limx[1])/2,arry_orient+arry_orient*0.05,"THF",ha='center')
plt.axvline(x=int_start, color = '#989898')
plt.axvline(x=int_end, color = '#989898')
plt.legend()
plt.ylabel("R / \AA")
plt.xlabel("Z / \AA")
plt.xlim(limx_all)
plt.ylim(limy_dst)
"""
plt.savefig(os.path.join(root_ORIENT,'IMC_Mg-Cl_dist.pgf'), format='pgf',bbox_inches = "tight")
plt.close()
# %%
