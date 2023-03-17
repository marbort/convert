#%%
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
from MDAnalysis.analysis import *
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.lineardensity import LinearDensity
import glob
import numpy as np
from scipy.signal import argrelextrema
from natsort import natsorted, ns
import subprocess
import regex as re
from scipy.signal import find_peaks

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

#%%


root="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
mols="IMC"
sel1="resname IMC and name Mg1"
sel2="resname IMC and name Cl1"
paths=glob.glob(root+mols+'/umbrella_30_200/window*')
paths_sorted=natsorted(paths)


#%%
profile=root+mols+"/umbrella_30_200/profile_all.xvg"
histo=root+mols+"/umbrella_30_200/histo_all.xvg"
data=[]
data_hist=[]
chars=["#","@"]
names=[]
rho=[]

with open(profile,'r') as ifile:
        lines=ifile.readlines()
        #x=[]
        #y=[]
        #for i in lines:
            #if i[0] in chars:
             #   continue
            #else:
                #x.append(i.split()[0])
                #y.append(i.split()[0])
        x=[float(x.split()[0]) for x in lines if x[0] not in chars ]
        y=[float(x.split()[1]) for x in lines if x[0] not in chars ]
        
        x_corr=[j+5.14 for j in x]
        y_scaled=[i-min(y) for i in y]
        #frcs_tmp=[[float(x.split()[i]) for x in lines[1:]] for i in range(6)]
        data.append([x_corr,y_scaled])
        
#%% ORIENTATION
us=[] 
for path in paths_sorted:
        if os.path.exists(os.path.join(path,'ORIENTATION_OK')):
                print("{} already done. Skipping...".format(path))
        else:
                files=glob.glob(os.path.join(path,'*.xtc'))
                dirs=[os.path.dirname(x) for x in files] 
                fnames=[os.path.basename(x).replace('.xtc','') for x in files]
                dumps=[os.path.basename(x).replace('_whole.xtc','.tpr.dump') for x in files]
                us=[mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",x) for x in files]
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
#%%

#%% DENSITIES
U=mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",files[0])
GCL_all_1=U.select_atoms("resname GCL")
CHL_all_1=U.select_atoms("resname CHL")
THF_all_1=U.select_atoms("resname THF")
Clm_all_1=U.select_atoms("resname Clm")

GCL_dens=LinearDensity(GCL_all_1,binsize=1)
CHL_dens=LinearDensity(CHL_all_1,binsize=1)
THF_dens=LinearDensity(THF_all_1,binsize=1)
Clm_dens=LinearDensity(Clm_all_1,binsize=1)
#densities=[LinearDensity(x,binsize=1) for x in us]
GCL_dens_res=GCL_dens.run()
CHL_dens_res=CHL_dens.run()
THF_dens_res=THF_dens.run()
Clm_dens_res=Clm_dens.run()
#DESd=LinearDensity(DES,binsize=1)
#DESd.run()
#%%
#pos_umbrella=[9.435,9.178,8.888,8.65,8.372,8.11,7.85,7.572,7.307,7.038,6.771,
#6.507,6.238,5.973,5.711,5.448,5.168,4.9,4.647,4.391,4.105,3.839,3.583,3.312,
#3.01,2.789,2.518,2.266,1.981,1.716]l
ang_avg_all=[]
std_ang_all=[]
pos_umbrella=[]
dens_CHL=[]
dens_GCL=[]
dens_Clm=[]
for path in paths_sorted:
        fnames=glob.glob(os.path.join(path,'*_orientation.dat'))
        
        ang_mean_wd={}
        ang_std_wd={}
        nframes=[]
        for file in fnames:
            with open(file,'r') as ifile:
                lines=ifile.readlines()
                nframes=[(len(lines)-6)]
                pos=[float(line.split()[2]) for line in lines if "##POS UMBRELLA##" in line]
                with open(file.replace("_whole_orientation.dat",".tpr_CHL_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_CHL.append(rho[tgt])
                with open(file.replace("_whole_orientation.dat",".tpr_GCL_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_GCL.append(rho[tgt])
                with open(file.replace("_whole_orientation.dat",".tpr_Clm_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_Clm.append(rho[tgt])
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
tuls=[(pos_umbrella[i],ang_avg_all[i],std_ang_all[i],dens_CHL[i],dens_GCL[i],dens_Clm[i]) for i in range(len(pos_umbrella))]
tuls.sort()
with open(os.path.join(root,mols,'orientation_results.dat'),'w') as ofile:
        ofile.write("{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}\n".format("Pos","angX","stdX","angY","stdY","angZ","stdZ"))
        for x in tuls:
                ofile.write("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(x[0],x[1][0],x[2][0],
                                                                                            x[1][1],x[2][1],
                                                                                            x[1][2],x[2][2]))
        
        
#print(pos_umbrella)


#%% SELECT ONLY FILES IN tprfiles.dat

ang_avg_all=[]
std_ang_all=[]
pos_umbrella=[]
dens_CHL=[]
dens_GCL=[]
dens_Clm=[]
with open(os.path.join(root,mols,'umbrella_30_200','tprfiles_DES.dat'),'r') as ifile:
        lines=ifile.readlines()
        fnames=[os.path.join(root,mols,'umbrella_30_200',x.replace('.tpr','_whole_orientation.dat').rstrip()) for x in lines]
ang_mean_wd={}
ang_std_wd={}
nframes=[]
for file in fnames:
        with open(file,'r') as ifile:
                lines=ifile.readlines()
                nframes=[(len(lines)-6)]
                pos=[float(line.split()[2]) for line in lines if "##POS UMBRELLA##" in line]
                with open(file.replace("_whole_orientation.dat",".tpr_CHL_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_CHL.append(rho[tgt])
                with open(file.replace("_whole_orientation.dat",".tpr_GCL_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_GCL.append(rho[tgt])
                with open(file.replace("_whole_orientation.dat",".tpr_Clm_density.xvg"),'r') as dfile:
                        lines2=dfile.readlines()
                        z=[float(x.split()[0]) for x in lines2 if "#" not  in x if "@" not in x]
                        rho=[float(x.split()[1]) for x in lines2 if "#" not  in x if "@" not in x]
                        diff=[abs(x-pos[0]-5.14) for x in z]
                        tgt=diff.index(min(diff))
                        dens_Clm.append(rho[tgt])
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
tuls=[(pos_umbrella[i],ang_avg_all[i],std_ang_all[i],dens_CHL[i],dens_GCL[i],dens_Clm[i]) for i in range(len(pos_umbrella))]
tuls.sort()
with open(os.path.join(root,mols,'orientation_results.dat'),'w') as ofile:
        ofile.write("{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}{:>8s}\n".format("Pos","angX","stdX","angY","stdY","angZ","stdZ"))
        for x in tuls:
                ofile.write("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(x[0],x[1][0],x[2][0],
                                                                                            x[1][1],x[2][1],
                                                                                            x[1][2],x[2][2]))
        
        
#print(pos_umbrella)


#%%
RDFS={"Mg":{},"Cl":{}}
COORD={"Mg":{},"Cl":{}}

for traj in paths_sorted:
        print(traj)
        dir=os.path.dirname(traj)
        window=os.path.basename(traj)
        pairsMg=[]
        pairsCl=[]
        RDFS['Mg'][window]={}
        RDFS['Cl'][window]={}
        COORD['Mg'][window]={}
        COORD['Cl'][window]={}
        #RDFS
        with open(os.path.join(traj,'RDF_Mg.xvg'),'r') as ifile:
                lines=ifile.readlines()
                for line in lines:
                        if re.findall('s\d',line) :
                                pairsMg.append(line.split()[-1])
                x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
                RDFS['dist']=x
                for i,t in enumerate(pairsMg):
                        RDFS['Mg'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        with open(os.path.join(traj,'RDF_Cl.xvg'),'r') as ifile:
                lines=ifile.readlines()
                for line in lines:
                        if re.findall('s\d',line) :
                                pairsCl.append(line.split()[-1])
                for i,t in enumerate(pairsCl):
                        RDFS['Cl'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        #COORDINATION NUMBER                
        with open(os.path.join(traj,'coord_Mg.xvg'),'r') as ifile:
                lines=ifile.readlines()
                x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
                COORD['dist']=x
                for i,t in enumerate(pairsMg):
                        COORD['Mg'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        with open(os.path.join(traj,'coord_Cl.xvg'),'r') as ifile:
                lines=ifile.readlines()
                for i,t in enumerate(pairsCl):
                        COORD['Cl'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]


#%%
fig = plt.figure(dpi=150)
ax = fig.add_subplot(projection='3d')
for i,item in enumerate(RDFS['Mg']):
        print(item)
        ax.plot(RDFS['dist'], [(x+1)/(x+1)+i for x in RDFS['dist']], RDFS['Mg'][item][pairs[0]], zdir='z',color='b')
plt.ylabel("WINDOW")
plt.show()

#%%
fig = plt.figure(dpi=150)
for l,i in enumerate(RDFS['Mg']):
#%%PLOT RDF Mg
clrs={}
ref=["blue","red","green","purple"]
for i,item in enumerate(pairsMg):
        clrs[item]=ref[i]
for j in pairsMg:
        clr=clrs[j]
        fig=plt.figure(figsize=(60,30),dpi=150)
        for k,i in enumerate(RDFS['Mg']):
                Coord=True
                try:
                        y_min=find_peaks([-x for x in RDFS['Mg'][i][j]],width=1)
                        x_min=RDFS['dist'][y_min[0][0]]
                except:
                        Coord=False
                ax = plt.subplot(4,8,k+1)
                ax.set_ylim([0,30])
                ax.plot(RDFS['dist'],RDFS['Mg'][i][j],color=clr,label="{}".format(j))   
                ax.scatter(x_min,RDFS['Mg'][i][j][RDFS['dist'].index(x_min)],color=clr)
                ax.set_title("window{}".format(k+1))             
                ax.legend(loc='lower right')
                if Coord:
                        c=ax.inset_axes([.5,.5,.48,.48])
                        c.plot(COORD['dist'],COORD['Mg'][i][j],color=clr) 
                        c.scatter(x_min,COORD['Mg'][i][j][RDFS['dist'].index(x_min)],color=clr)
                        c.text(x_min,COORD['Mg'][i][j][RDFS['dist'].index(x_min)]+(COORD['Mg'][i][j][RDFS['dist'].index(x_min)]*0.5),"{:.3f}".format(COORD['Mg'][i][j][RDFS['dist'].index(x_min)]) )
        plt.suptitle(j,size=24)
        plt.show()
#print(COORD['Mg']['window1'][pairs[0]][RDFS['dist'].index(x_min)])
#%% CL RDFS PLOTS
clrs={}
ref=["blue","red","green","purple"]
for i,item in enumerate(pairsCl):
        clrs[item]=ref[i]
for j in pairsCl:
        clr=clrs[j]
        fig=plt.figure(figsize=(35,30),dpi=150)
        for k,i in enumerate(RDFS['Cl']):
                Coord=True
                try:
                        y_min=find_peaks([-x for x in RDFS['Cl'][i][j]],width=1)
                        x_min=RDFS['dist'][y_min[0][0]]
                except:
                        Coord=False
                ax = plt.subplot(5,6,k+1)
                ax.set_ylim([0,30])
                ax.plot(RDFS['dist'],RDFS['Cl'][i][j],color=clr,label="{}".format(j))   
                ax.scatter(x_min,RDFS['Cl'][i][j][RDFS['dist'].index(x_min)],color=clr)
                ax.set_title("window{}".format(k+1))             
                ax.legend(loc='lower right')
                if Coord:
                        c=ax.inset_axes([.5,.5,.48,.48])
                        c.plot(COORD['dist'],COORD['Cl'][i][j],color=clr) 
                        c.scatter(x_min,COORD['Cl'][i][j][RDFS['dist'].index(x_min)],color=clr)
                        c.text(x_min,COORD['Cl'][i][j][RDFS['dist'].index(x_min)]*2,"{:.3f}".format(COORD['Cl'][i][j][RDFS['dist'].index(x_min)]) )
        plt.suptitle(j,size=44)
        plt.show()
                                
#%%
posz_all=[x[0]+5.14 for x in tuls]
limy_all=[0,180]
limx_all=[6.8,12]
axis=["X","Y","Z"]

CN={'Mg':{},'Cl':{}}
for k in RDFS:
        if k == 'dist':
                continue
        else:
                for numb,i in enumerate(RDFS[k]):
                        CN[k][i]={'pos':posz_all[-1-numb]}
                        for j in RDFS[k][i]:
                                try:
                                        y_min=find_peaks([-x for x in RDFS[k][i][j]],width=1)
                                        x_min=RDFS['dist'][y_min[0][0]]
                                        CN[k][i][j]=COORD[k][i][j][RDFS['dist'].index(x_min)]
                                except:
                                        CN[k][i][j]=0


#%%


fig=plt.figure(figsize=(10,15),dpi=150)

for i in range(3):
        plt.subplot(7,1,i+1)
        plt.errorbar(posz_all,[x[1][i] for x in tuls],[x[2][i] for x in tuls],linestyle=None)
        plt.plot(posz_all,[x[1][i] for x in tuls])
        plt.axvline(x=10.6, color = '#989898')
        plt.axvline(x=8.4, color = '#989898')
        plt.axhline(y=90, color = '#989898')
        plt.ylim(limy_all)
        plt.xlim(limx_all)
        plt.ylabel("{} Angle / deg".format(axis[i]))
plt.subplot(7,1,4)
for pair in pairsMg:
        plt.plot([CN['Mg'][x]['pos'] for x in CN['Mg']],[CN['Mg'][x][pair] for x in CN['Mg']],'o-',label=pair)
        plt.ylabel("Coord Number")
plt.plot([CN['Mg'][x]['pos'] for x in CN['Mg']],[CN['Mg'][x][pairsMg[0]] + CN['Mg'][x][pairsMg[1]]+CN['Mg'][x][pairsMg[2]]+CN['Mg'][x][pairsMg[3]] for x in CN['Mg']],'-o',label="Sum")
plt.ylim([-0.5,4.5])
plt.xlim(limx_all)
plt.legend(location='lower right')
plt.subplot(7,1,5)
for pair in pairsCl:
        plt.plot([CN['Cl'][x]['pos'] for x in CN['Cl']],[CN['Cl'][x][pair] for x in CN['Cl']],'o-',label=pair)
        plt.ylabel("Coord Number")
plt.ylim([-0.5,3.5])
plt.xlim(limx_all)
plt.legend()
plt.subplot(7,1,7)
plt.plot(data[0][0],data[0][1])
plt.axvline(x=10.6, color = '#989898')
plt.axvline(x=8.4, color = '#989898')
plt.ylabel("Free Energy / Kj/mol")
plt.xlim(limx_all)
plt.subplot(7,1,6)
plt.plot(posz_all,[x[3] for x in tuls],'-o',label="CHL density")
plt.plot(posz_all,[x[4] for x in tuls],'-o',label="GCL density")
plt.plot(posz_all,[x[5] for x in tuls],'-o',label="Clm density")
#plt.plot(GCL_dens_res.results.z.hist_bin_edges[:-1]/10 ,GCL_dens_res.results.z.mass_density)
#plt.plot(CHL_dens_res.results.z.hist_bin_edges[:-1]/10 ,CHL_dens_res.results.z.mass_density)
#plt.plot(THF_dens_res.results.z.hist_bin_edges[:-1]/10 ,THF_dens_res.results.z.mass_density)
#plt.plot(Clm_dens_res.results.z.hist_bin_edges[:-1]/10 ,Clm_dens_res.results.z.mass_density)
plt.axvline(x=10.6, color = '#989898')
plt.axvline(x=8.4, color = '#989898')
plt.xlim(limx_all)
plt.ylabel("Density / kg/m3")
plt.legend()
#plt.legend(["GCL","CHL","THF","Clm"])
plt.xlabel('$z\ /\ nm$')
plt.suptitle(mols)
plt.tight_layout()
# %%
with open('/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/window29/angles_cpptraj.dat','r') as ifile:
        lines=ifile.readlines()
        angls=[float(x.split()[1]) for x in lines[1:]]
ag=np.mean(angls)
std_ag=np.std(angls)
print(ag,std_ag)
# %%DENS ACCURATE WINDOW
posz=9.611+5.14
with open(os.path.join(root,mols,"umbrella_30_200/window1/density.xvg"),'r') as ifile:
        lines=ifile.readlines()
        z=[float(x.split()[0]) for x in lines if "#" not  in x if "@" not in x]
        rho=[float(x.split()[1]) for x in lines if "#" not  in x if "@" not in x]
diff=[abs(x-posz) for x in z]
tgt=diff.index(min(diff))
print(rho[tgt])
        



# %%
