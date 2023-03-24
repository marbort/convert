#%%
import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
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
def coord_Number(x,y,dist_min,dist_max,mols):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4/3*np.pi*np.trapz(prod,val_x)
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
    clrs=['blue','cyan','tan','violet']
    for i,x in enumerate(HB):
        plt.plot(HB[x]['frame'],[j/HB_tot for j in HB[x]['HB']],label=lbls[x],color=clrs[i])
        print("Average HB percentage of {} = {}".format(x,np.mean([j/HB_tot for j in HB[x]['HB']])))
    plt.legend()
    plt.ylim([0,1])
    plt.xlabel("Frame")
    plt.ylabel("H-Bond fraction")
    plt.savefig(os.path.join(folder,'HB.pgf'),format='pgf',bbox_inches = "tight")
    
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
        RDFs2[i]['CN']=coord_Number(RDFs2[i]['RDF'][0],RDFs2[i]['RDF'][1],0,RDFs2[i]['r_min'],RDFs2[i]['#_part'])*D2[i]['density']
    #print(RDFs['GCL_C2_GCL_C2'])

    with open (os.path.join(root,'results_COM.dat'),'w') as ofile:
        ofile.write("{:35s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN'))
        for i in RDFs2:
            ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f}\n".format(RDFs2[i]['name'],RDFs2[i]['r_max'],RDFs2[i]['r_min'],RDFs2[i]['CN']))
    wd,hg=set_size(textsize,fraction)
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    clrs=['blue','orange','green','red','purple','black']
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
    plt.savefig(os.path.join(root,tit+'.pgf'), format='pgf',bbox_inches = "tight")
#%%
def plot_DENS_cpptraj(root:str,textsize:int,fraction:float,labels:list,title:str="",a:float=0,b:float=0):
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
            plt.plot(dens[file][vars[0]],dens[file][vars[var]],label=''.join(e for e in vars[var] if e.isalnum()))
        plt.xlabel(" Z / \AA")
        plt.ylabel("$ Density / Kg m^{-3}$ ")
        plt.legend()
        plt.axvspan(a, b, color='#989898', alpha=0.5, lw=0)
        plt.title(title)
        plt.savefig(os.path.join(root,title+'.pgf'), format='pgf',bbox_inches = "tight")
    
        
    return(dens,vars,fig)
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
            x_corr=[j+51.4 for j in x]
            y_scaled=[i-min(y) for i in y]
            data.append([x_corr,y_scaled])
    with open(os.path.join(root,histo),'r') as ifile:
            lines=[x for x in ifile.readlines() if "#" not in x if "@" not in x]
            data_hist.append([float(x.split()[0])*10+51.4 for x in lines])
            for i in range(1,len(lines[0].split())):
                data_hist.append([float(x.split()[i]) for x in lines])
            
            
    return(data,data_hist)

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
    for path in paths_sorted:
            if os.path.exists(os.path.join(path,'ORIENTATION_OK')):
                    print("{} already done. Skipping...".format(path))
            else:
                    files=glob.glob(os.path.join(path,'*.xtc'))
                    dirs=[os.path.dirname(x) for x in files] 
                    fnames=[os.path.basename(x).replace('.xtc','') for x in files]
                    dumps=[os.path.basename(x).replace('_whole.xtc','.tpr.dump') for x in files]
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
    posz_all=[x[0]*10+51.4 for x in tuls]
    axis=["X","Y","Z"]
    CN={x:{} for x in RDFS}
    for k in RDFS:
        for numb,i in enumerate(RDFS[k]):
            if i != 'dist':
                try:
                    CN[k][i]={'pos':posz_all[-1-numb]}
                except:
                    continue
                for j in RDFS[k][i]:
                        try:
                                y_min=find_peaks([-x for x in RDFS[k][i][j]],width=1)
                                x_min=RDFS[k]['dist'][y_min[0][0]]
                                CN[k][i][j]=COORD[k][i][j][RDFS[k]['dist'].index(x_min)]
                        except:
                               CN[k][i][j]=0
    return(CN)
                                        
#%%

root_HB="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/DES"
plot_hbonds(root_HB,345,1,["Choline-Cl","Gly-Gly","Gly-Cl","choline-Gly"])
# %%
root_COM="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/DES"
resnumb={"CHL":400,"Clm":400,"GCL":800}
plot_COM_RDF(root_COM,345,1,["Choline-Choline","Choline-chloride","Choline-glycerol","Glycerol-chloride","Glycerol-glycerol",'Chloride-Chloride','Chloride-Chloride'],"0.90",resnumb)

# %%
root_DENS="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/AcPh"
D,v,fig=plot_DENS_cpptraj(root_DENS,345,1,["Acetophenone"],"rho-ACP",84,106)
# %% PLOT PMFS
root_PMF="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
mols=["IMC","ACT","ACP"]
dir="umbrella_30_200"
textsize=426
fraction=1
a=84
b=106
limx=[65,120]
limy=[-1,40]
graph_data=[extract_umbrella_profile(os.path.join(root_PMF,x,dir),'profile_all.xvg','hist_all.xvg') for x in mols]
wd,hg=set_size(textsize,fraction)
fig=plt.figure(figsize=(wd,hg*2),dpi=150)
for i,x in enumerate(graph_data):
    plt.subplot(3,1,i+1)
    plt.plot(x[0][0][0],x[0][0][1],label=mols[i])
    plt.xlim(limx)
    plt.ylim(limy)
    plt.legend()
    plt.axvspan(a, b, color='#989898', alpha=0.5, lw=0)
    plt.ylabel("Free Energy / $Kjmol^{-1}$")
    plt.tight_layout
plt.xlabel("Z / \AA")


plt.savefig(os.path.join(root_PMF,'_'.join(mols)+'_PMF.pgf'), format='pgf',bbox_inches = "tight")

#PLOT HISTOS
fig=plt.figure(figsize=(wd,hg*2),dpi=150)
for i,x in enumerate(graph_data):
    plt.subplot(3,1,i+1)
    for j in range(1,len(x[1])):
        plt.plot(x[1][0],x[1][j])
    plt.xlim(limx)
    plt.legend(mols[i])
    #plt.ylim(limy)
    #plt.legend()
    plt.axvspan(a, b, color='#989898', alpha=0.5, lw=0)
    plt.ylabel("Free Energy / $Kjmol^{-1}$")
    plt.tight_layout
plt.xlabel("Z / \AA")


plt.savefig(os.path.join(root_PMF,'_'.join(mols)+'_HIST.pgf'), format='pgf',bbox_inches = "tight")



# %%
plot_image('/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/img/interface_struct.pov.png',426,1)
# %%extract ORIENTATION
tuls=select_orientation('/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/tprfiles_DES.dat','/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/orientation_results.dat')
mols,pairs,RDFS,COORD=extract_RDF_gmx('/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200',["Mg","Cl"])
orientation_PROFILE=extract_umbrella_profile('/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200','profile_all.xvg','hist_all.xvg')
CN=extract_coordination_Number(tuls,RDFS)
# %% PLOT ORIENTATION

wd,hg=set_size(426,1)
fig=plt.figure(figsize=(wd,hg*2),dpi=150)

root_ORIENT='/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200'
limx_all=[65,125]
limy_all=[0,40]
ltr=['a','b','c','d','e']
ax=plt.subplot(4,1,1)
plt.errorbar([x[0]*10+51.4 for x in tuls],[x[1][2] for x in tuls],[x[2][2] for x in tuls],linestyle=None)
plt.plot([x[0]*10+51.4 for x in tuls],[x[1][2] for x in tuls])
plt.text(.9,.9,'a',transform=ax.transAxes)
plt.axvline(x=106, color = '#989898')
plt.axvline(x=84, color = '#989898')
plt.axhline(y=90, color = '#989898')
plt.xlim(limx_all)
plt.ylabel("Z Angle / deg")
for i,mol in enumerate(mols):
    ax=plt.subplot(4,1,2+i)
    plt.text(.9,.9,ltr[1+i],transform=ax.transAxes)
    for pair in pairs[i]:
            plt.plot([CN[mol][x]['pos'] for x in CN[mol]],[CN[mol][x][pair] for x in CN[mol] if x != 'dist'],'o-',label=''.join(e for e in pair if e.isalnum()))
            plt.ylabel("Coord Number")
            plt.axvline(x=106, color = '#989898')
            plt.axvline(x=84, color = '#989898')
            plt.legend(loc='lower right')
            plt.xlim(limx_all)
#plt.plot([CN['Mg'][x]['pos'] for x in CN['Mg']],[CN['Mg'][x][pairsMg[0]] + CN['Mg'][x][pairsMg[1]]+CN['Mg'][x][pairsMg[2]]+CN['Mg'][x][pairsMg[3]] for x in CN['Mg']],'-o',label="Sum")
            plt.ylim([-0.5,4.5])


ax=plt.subplot(4,1,2+len(mols))
plt.text(.9,.9,ltr[1+len(mols)],transform=ax.transAxes)
plt.plot(orientation_PROFILE[0][0][0],orientation_PROFILE[0][0][1])
plt.axvline(x=106, color = '#989898')
plt.axvline(x=84, color = '#989898')
plt.ylabel("Free Energy / Kj/mol")
plt.xlim(limx_all)
plt.tight_layout()
plt.savefig(os.path.join(root_ORIENT,'_IMC_ORIENT.pgf'), format='pgf',bbox_inches = "tight")
plt.close()
# %%
