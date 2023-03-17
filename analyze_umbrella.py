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
import regex as re
from scipy.signal import find_peaks
#%%
def coord_number(x,y,dist_min,dist_max,mols):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*np.trapz(prod,val_x)
    CN_one=CN_all/mols
    return(CN_one)

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
w_THF=range(1,17)
w_INT=range(15,22)
w_DES=range(17,31)
w_all=range(1,31)
root="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
mols="IMC"
paths_THF=[root+mols+"/umbrella_30_200/window{}/umbrella_whole.xtc".format(i) for i in w_THF]
paths_INT=[root+mols+"/umbrella_30_200/window{}/umbrella_whole.xtc".format(i) for i in w_INT]
paths_DES=[root+mols+"/umbrella_30_200/window{}/umbrella_DES_whole.xtc".format(i) for i in w_DES]
paths_ALL=[root+mols+"/umbrella_30_200/window{}/umbrella_whole.xtc".format(i) for i in w_all]

paths_ALL_mix=paths_THF+paths_DES

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



print(paths_ALL_mix)
#%%
# Load the LAMMPS trajectory
#for i in paths:
u_THF = mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",paths_THF)
u_INT = mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",paths_INT)
u_DES = mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",paths_DES)
u_INT_sigle=[mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",i) for i in paths_INT]
#%%
u_single=[mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",i) for i in paths_ALL_mix]
#%%

GCL_all_1=u_single[0].select_atoms('resname GCL')
CHL_all_1=u_single[0].select_atoms('resname CHL')
THF_all_1=u_single[0].select_atoms('resname THF')
Clm_all_1=u_single[0].select_atoms('resname Clm')

# Select the atoms to be used in the calculation THF
CHL = u_THF.select_atoms('resname CHL Clm')
GCL = u_THF.select_atoms('resname GCL')
THF = u_THF.select_atoms('resname THF')
ACT = u_THF.select_atoms('resname ACT')
IPM = u_THF.select_atoms('resname IPM')
IMC = u_THF.select_atoms('resname IMC')
MGCL = u_THF.select_atoms('resname MG Cl-')
us=[u_THF.select_atoms('resname {}'.format(i)) for i in mols]
DES = u_DES.select_atoms('resname CHL GCL Clm')

CHLd = u_DES.select_atoms('resname CHL Clm')
GCLd = u_DES.select_atoms('resname GCL')
THFd = u_DES.select_atoms('resname THF')
ACTd = u_DES.select_atoms('resname ACT')
IMCd = u_DES.select_atoms('resname IMC')
MGCLd = u_DES.select_atoms('resname MG Cl-')
usd=[u_DES.select_atoms('resname {}'.format(i)) for i in mols]


#DESd = u.select_atoms('resname CHL GCL Clm')



#%% RDF THF
bins=80
CHL1_ho=CHL.select_atoms('resid 1 and type ho')
CHL1_oh=CHL.select_atoms('resid 1 and type oh')
CHL1_c=CHL.select_atoms('resid 1 and type c3')
GCL_oh=GCL.select_atoms('type oh')
GCL_ho=GCL.select_atoms('type ho')
THF_o =THF.select_atoms('type os')
THF_c = THF.select_atoms('type c3')
IMC_Cl=IMC.select_atoms('type cl')
IMC_Mg=IMC.select_atoms('type mg')


#%% DENSITIES
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

#%%ORIENTATION THF
vec=[IMC_Mg.positions - IMC_Cl.positions[i] for i,x in enumerate(IMC_Mg.positions)]
#%%
#print(vec)
IMC_Mg_traj_pos=[IMC_Mg.positions for ts in u_THF.trajectory]
IMC_Cl_traj_pos=[IMC_Cl.positions for ts in u_THF.trajectory]
#%%
dist_trj=[x-IMC_Mg_traj_pos[i][0] for i,x in enumerate(IMC_Cl_traj_pos)]
#print(dist_trj[1][0])
norm=[np.linalg.norm(x) for x in dist_trj]
print(norm[0])
cos=[[k[0][j]/norm[i] for j in range(3) ] for i,k in enumerate(dist_trj)]
ang=[[np.arccos(k[j])*57.30 for j in range(3) ] for k in cos]
#print(ang[2])
print(np.mean([x[1] for x in cos]))
#%%
cos_func,ang_func=angle_with_axes("resname IMC and type mg","resname IMC and type cl",u_THF)
print(np.mean([x[1] for x in cos_func]))


#%%
distances = []
IMC2_cl=IMC.select_atoms('type cl')
IMC2_mg=IMC.select_atoms('type mg')
print(IMC2_mg.atoms.n_atoms)
print(IMC2_cl.atoms.n_atoms)
#%%
for ts in u_THF.trajectory:
    dist = mda.analysis.distances.dist(IMC2_mg, IMC2_cl)
    distances.append(dist)
print(distances[0])
#%%

#print(np.arccos(np.mean(cos[:][0]))*57.29)
#print(np.mean(cos[:][1]))
#print(np.mean(cos[:][2]))
#print(cos[1][1])
#print(np.mean([x[2] for x in cos]))
#plt.plot([i for i in range(6370)],[x[2] for x in ang])
#plt.show()
#print(IMC_Cl_traj_pos[1][0]-IMC_Mg_traj_pos[1][0])

#print(len(u_THF.trajectory))

#%%
cos_angle_THF=angle_with_axes(IMC_Mg,IMC_Cl,u_THF)
print(np.arccos(np.mean(cos_angle_THF[0])),np.arccos(np.mean(cos_angle_THF[1])),np.arccos(np.mean(cos_angle_THF[2])))


#%%
RDF_CHL1_ho_GCL_oh=mda.analysis.rdf.InterRDF(CHL1_ho, GCL_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1_oh_GCL_ho=mda.analysis.rdf.InterRDF(CHL1_oh, GCL_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1_oh_THF_o=mda.analysis.rdf.InterRDF(CHL1_oh, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1_c_THF_c=mda.analysis.rdf.InterRDF(CHL1_c, THF_c, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)



RDF_CHL1_ho_GCL_oh.run()
RDF_CHL1_oh_GCL_ho.run()
RDF_CHL1_oh_THF_o.run()
RDF_CHL1_c_THF_c.run()

#A1a=mda.analysis.bat.BAT(A1)
#%% RDF IMC THF
#RDF_IMC_Cl_THF_h=mda.analysis.rdf.InterRDF(IMC_Mg, THF_h, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_Mg_THF_o=mda.analysis.rdf.InterRDF(IMC_Mg, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm='density')
#RDF_IMC_Cl_IMC_Mg=mda.analysis.rdf.InterRDF(IMC_Cl, IMC_Mg, nbins=bins, range=(1.0, 8.0), exclusion_block=(2,2),norm=normp)
RDF_IMC_Mg_THF_o.run()
fig=plt.figure(dpi=150)
#plt.plot(RDF_IMC_Cl_IMC_Mg.results.bins,RDF_IMC_Cl_IMC_Mg.results.rdf)
plt.plot(RDF_IMC_Mg_THF_o.results.bins,RDF_IMC_Mg_THF_o.results.rdf)
#plt.plot(RDF_IMC_Cl_THF_h.results.bins,[x/4 for x in RDF_IMC_Cl_THF_h.results.rdf])
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_cl...IMC_mg","IMC_mg...THF_o","IMC_Cl...THF_h"],prop={'size': 14})
#%%
minx=[x for x in RDF_IMC_Mg_THF_o.results.bins if 2 < x < 3]
miny=[RDF_IMC_Mg_THF_o.results.rdf[i] for i,item in enumerate(RDF_IMC_Mg_THF_o.results.bins) if 2 < item < 3]
ing=np.trapz(miny,x=minx)
print(ing)
#RDF_IMC_Cl_IMC_Mg.run()
#RDF_IMC_Cl_THF_h.run()
#%% RDF DES
bins=80
normp="density"
IMC1d_ch=u_DES.select_atoms('resid 1 and name C1 ')
IMC1d_cl=u_DES.select_atoms('resid 1 and type cl ')
IMC1d_mg=u_DES.select_atoms('resid 1 and type mg ')
CHL1d_ho=CHLd.select_atoms('resid 1 and type ho')
CHL1d_oh=CHLd.select_atoms('resid 1 and type oh')
CHL1d_c=CHLd.select_atoms('resid 1 and type c3')
GCLd_oh=GCLd.select_atoms('type oh')
GCLd_ho=GCLd.select_atoms('type ho')
THFd_o =THFd.select_atoms('type os')
THFd_c = THFd.select_atoms('type c3')


#%% ORIENTATION DES
cos_func_DES,ang_func_DES=angle_with_axes("resname IMC and type mg","resname IMC and type cl",u_DES)
print(np.mean([x[2] for x in ang_func_DES]))



#%%

RDF_CHL1d_ho_GCLd_oh=mda.analysis.rdf.InterRDF(CHL1d_ho, GCLd_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1d_oh_GCLd_ho=mda.analysis.rdf.InterRDF(CHL1d_oh, GCLd_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1d_oh_THFd_o=mda.analysis.rdf.InterRDF(CHL1d_oh, THFd_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL1d_c_THFd_c=mda.analysis.rdf.InterRDF(CHL1d_c, THFd_c, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1d_ch_GCL_ho=mda.analysis.rdf.InterRDF(IMC1d_ch, GCLd_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1d_cl_GCL_ho=mda.analysis.rdf.InterRDF(IMC1d_cl, GCLd_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1d_mg_GCL_oh=mda.analysis.rdf.InterRDF(IMC1d_mg, GCLd_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)



RDF_CHL1d_ho_GCLd_oh.run()
RDF_CHL1d_oh_GCLd_ho.run()
RDF_CHL1d_oh_THFd_o.run()
RDF_CHL1d_c_THFd_c.run()
RDF_IMC1d_ch_GCL_ho.run()
RDF_IMC1d_cl_GCL_ho.run()
RDF_IMC1d_mg_GCL_oh.run()

#%% RDF INT
bins=80
normp="density"
CHLi=u_INT.select_atoms("resname CHL")
GCLi=u_INT.select_atoms('resname GCL')
THFi=u_INT.select_atoms('resname THF')
IMC1i=u_INT.select_atoms('resname IMC')
IMC1i_ch=u_INT.select_atoms('resid 1 and name C1 ')
IMC1i_cl=u_INT.select_atoms('resid 1 and type cl ')
IMC1i_mg=u_INT.select_atoms('resid 1 and type mg ')
#CHL1i_ho=CHLi.select_atoms('resid 1 and type ho')
#CHL1i_oh=CHLi.select_atoms('resid 1 and type oh')
#CHL1i_c=CHLi.select_atoms('resid 1 and type c3')
GCLi_oh=GCLi.select_atoms('type oh')
GCLi_ho=GCLi.select_atoms('type ho')
THFi_o =THFi.select_atoms('type os')
THFi_c = THFi.select_atoms('type c3')
Clm_cl = u_INT.select_atoms('type Clm')

#%%ORIENTATION INT

cos_func_INT,ang_func_INT=angle_with_axes("resname IMC and type mg","resname IMC and type cl",u_INT)
print(np.mean([x[2] for x in ang_func_INT]))




#%% RDFs

#RDF_CHL1i_ho_GCLi_oh=mda.analysis.rdf.InterRDF(CHL1i_ho, GCLi_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
#RDF_CHL1i_oh_GCLi_ho=mda.analysis.rdf.InterRDF(CHL1i_oh, GCLi_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
#RDF_CHL1i_oh_THFi_o=mda.analysis.rdf.InterRDF(CHL1i_oh, THFi_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
#RDF_CHL1i_c_THFi_c=mda.analysis.rdf.InterRDF(CHL1i_c, THFi_c, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1i_ch_GCL_ho=mda.analysis.rdf.InterRDF(IMC1i_ch, GCLi_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1i_cl_GCL_ho=mda.analysis.rdf.InterRDF(IMC1i_cl, GCLi_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC1i_mg_GCL_oh=mda.analysis.rdf.InterRDF(IMC1i_mg, GCLi_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)



#RDF_CHL1i_ho_GCLi_oh.run()
#RDF_CHL1i_oh_GCLi_ho.run()
#RDF_CHL1i_oh_THFi_o.run()
#RDF_CHL1i_c_THFi_c.run()
RDF_IMC1i_ch_GCL_ho.run()
RDF_IMC1i_cl_GCL_ho.run()
RDF_IMC1i_mg_GCL_oh.run()
#A1a.run()

#%% RDF INT SINGLE
bins=80
normp="density"
print(len(u_INT_sigle))
CHLi_s=[x.select_atoms('resname CHL') for x in u_INT_sigle]
GCLi_s=[x.select_atoms('resname GCL') for x in u_INT_sigle]
THFi_s=[x.select_atoms('resname THF') for x in u_INT_sigle]
Clm_s=[x.select_atoms('type Clm') for x in u_INT_sigle]


IMC1i_ch_s=[x.select_atoms('resid 1 and name C1 ') for x in u_INT_sigle]
IMC1i_cl_s=[x.select_atoms('resid 1 and type cl ') for x in u_INT_sigle]
IMC1i_mg_s=[x.select_atoms('resid 1 and type mg ') for x in u_INT_sigle]


CHLi_oh_s=[x.select_atoms('type oh ') for x in CHLi_s]
CHLi_ho_s=[x.select_atoms('type ho ') for x in CHLi_s]

GCLi_oh_s=[x.select_atoms('type oh ') for x in GCLi_s]
GCLi_ho_s=[x.select_atoms('type ho ') for x in GCLi_s]
GCLi_c_s=[x.select_atoms('type c3 ') for x in GCLi_s]

THFi_o_s =[x.select_atoms('type os ') for x in THFi_s]
THFi_c_s =[x.select_atoms('type c3 ') for x in THFi_s] 



#%%orientation INT single

#print(u_INT_sigle[0])
orientation_INT=[angle_with_axes("resname IMC and type mg","resname IMC and type cl",x) for x in u_INT_sigle]
#cos_func_INT_s,ang_func_INT_s=[angle_with_axes("resname IMC and type mg","resname IMC and type cl",x) for x in u_INT_sigle]
#%%
#print(np.array(orientation_INT).shape)
print(orientation_INT[0][1][0])
np.mean([x[0]  for x in orientation_INT[0][0]] )
#%%
ang_avg=[]
std_ang_avg=[]
cos_avg=[]
std_cos_avg=[]
for i in orientation_INT:
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
print(ang_avg[0],std_ang_avg[0])
       
#%%
cos_avg=[[np.mean(k[0]),np.mean(k[1]),np.mean(k[2])] for k in x[0] for x in orientation_INT]
print(cos_avg)
len(np.mean(orientation_INT,axis=0))
#%% plot orientation
fig=plt.figure(figsize=(10,10),dpi=150)
rows=(len(u_INT_sigle)//3+1)
limxs=[1,8]
limys=[0,0.035]

dists=[5.94778,5.57157,5.32892,5.09535,4.93143,4.67588,4.37779]
posz=[x + 5.14 for x in dists]
for i in range(3):
        plt.subplot(3,1,i+1)
        plt.errorbar(posz,[x[i] for x in ang_avg],[x[i] for x in std_ang_avg])
plt.show()
        
        
#print(np.mean([x[2] for x in ang_func_DES]))

#%% ORIENTATION ALL 
orientation_ALL=[angle_with_axes("resname IMC and type mg","resname IMC and type cl",x) for x in u_single]

#%%
ang_avg=[]
std_ang_avg=[]
cos_avg=[]
std_cos_avg=[]
for index,i in enumerate(orientation_ALL):
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
        with open(root+mols+'/umbrella_30_200/angles_window_{}.dat'.format(index+1),'w') as ofile:
                ofile.write("{:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n".format("cosX","cosY","cosZ","angX","angY","angZ"))
                for idx,x in enumerate(cosx):
                        ofile.write("{:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(x,cosy[idx],cosz[idx],angx[idx],angy[idx],angz[idx]))
                


#%%


RDF_IMC1i_ch_GCL_ho_s=[mda.analysis.rdf.InterRDF(x, GCLi_ho_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_ch_s)]
RDF_IMC1i_mg_GCL_oh_s=[mda.analysis.rdf.InterRDF(x, GCLi_oh_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_mg_s)]
RDF_IMC1i_cl_GCL_ho_s=[mda.analysis.rdf.InterRDF(x, GCLi_ho_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_cl_s)]
RDF_IMC1i_ch_GCL_c_s=[mda.analysis.rdf.InterRDF(x, GCLi_c_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_ch_s)]
RDF_IMC1i_cl_CHL_ho_s=[mda.analysis.rdf.InterRDF(x, CHLi_ho_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_cl_s)]
RDF_IMC1i_mg_CHL_oh_s=[mda.analysis.rdf.InterRDF(x, CHLi_oh_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_mg_s)]
RDF_IMC1i_mg_THF_o_s=[mda.analysis.rdf.InterRDF(x, THFi_o_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_mg_s)]
RDF_IMC1i_mg_Clm_s=[mda.analysis.rdf.InterRDF(x, Clm_s[i], nbins=bins, range=(1.0, 8.0), 
                        exclusion_block=None,norm=normp) for i,x in enumerate(IMC1i_mg_s)]







#%%
RDF_IMC1i_ch_GCL_ho_res=[x.run() for x in RDF_IMC1i_ch_GCL_ho_s]
RDF_IMC1i_mg_GCL_oh_res=[x.run() for x in RDF_IMC1i_mg_GCL_oh_s]
RDF_IMC1i_cl_GCL_ho_res=[x.run() for x in RDF_IMC1i_cl_GCL_ho_s]
RDF_IMC1i_cl_CHL_ho_res=[x.run() for x in RDF_IMC1i_cl_CHL_ho_s]
RDF_IMC1i_mg_CHL_oh_res=[x.run() for x in RDF_IMC1i_mg_CHL_oh_s]
RDF_IMC1i_ch_GCL_c_res=[x.run() for x in RDF_IMC1i_ch_GCL_c_s]
RDF_IMC1i_mg_Clm_res=[x.run() for x in RDF_IMC1i_mg_Clm_s]
#%%
RDF_IMC1i_mg_THF_o_res=[x.run() for x in RDF_IMC1i_mg_THF_o_s]

#%%
n_res=1
dist=3.5
CN_IMC_mg_GCL_oh =[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_mg_GCL_oh_res]
CN_IMC_mg_CHL_oh =[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_mg_CHL_oh_res]
CN_IMC_mg_THF_o=[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_mg_CHL_oh_res]
CN_IMC1i_cl_CHL_ho=[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_cl_CHL_ho_res]
CN_IMC1i_cl_GCL_ho=[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_cl_GCL_ho_res]
CN_IMC1i_ch_GCL_ho=[coord_number(x.results.bins,x.results.rdf,0,dist,n_res) for x in RDF_IMC1i_ch_GCL_ho_res]

min_ind=argrelextrema(RDF_IMC1i_mg_GCL_oh_res[0].results.rdf, np.less)
print(min_ind[0][0])
print(RDF_IMC1i_mg_GCL_oh_res[0].results.bins[min_ind[0][:5]])
#%%
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL1_oh_THF_o.results.bins,RDF_CHL1_oh_THF_o.results.rdf)
plt.plot(RDF_CHL1_c_THF_c.results.bins,RDF_CHL1_c_THF_c.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_oh...THF_o","CHL_c...THF_c"],prop={'size': 14})
#%%
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL1d_oh_THFd_o.results.bins,RDF_CHL1d_oh_THFd_o.results.rdf)
plt.plot(RDF_CHL1d_c_THFd_c.results.bins,RDF_CHL1d_c_THFd_c.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_oh...THF_o","CHL_c...THF_c"],prop={'size': 14})

#%% PLOT RDF IMC DES
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC1d_ch_GCL_ho.results.bins,[x/3 for x in RDF_IMC1d_ch_GCL_ho.results.rdf])
plt.plot(RDF_IMC1d_cl_GCL_ho.results.bins,[x/3 for x in RDF_IMC1d_cl_GCL_ho.results.rdf])
plt.plot(RDF_IMC1d_mg_GCL_oh.results.bins,[x/3 for x in RDF_IMC1d_mg_GCL_oh.results.rdf])
#plt.plot(RDF_CHL1d_c_THFd_c.results.bins,RDF_CHL1d_c_THFd_c.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_ch...GCL_ho","IMC_cl...GCL_ho","IMC_mg...GCL_oh"],prop={'size': 14})


#%% PLOT RDF IMC INT
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC1i_ch_GCL_ho.results.bins,[x/3 for x in RDF_IMC1i_ch_GCL_ho.results.rdf])
plt.plot(RDF_IMC1i_cl_GCL_ho.results.bins,[x/3 for x in RDF_IMC1i_cl_GCL_ho.results.rdf])
plt.plot(RDF_IMC1i_mg_GCL_oh.results.bins,[x/3 for x in RDF_IMC1i_mg_GCL_oh.results.rdf])
#plt.plot(RDF_CHL1d_c_THFd_c.results.bins,RDF_CHL1d_c_THFd_c.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_ch...GCL_ho","IMC_cl...GCL_ho","IMC_mg...GCL_oh"],prop={'size': 14})

#%%PLOT IMC RDF
fig=plt.figure(dpi=150)
#plt.plot(RDF_IMC_Cl_IMC_Mg.results.bins,RDF_IMC_Cl_IMC_Mg.results.rdf)
plt.plot(RDF_IMC_Mg_THF_o.results.bins,RDF_IMC_Mg_THF_o.results.rdf)
#plt.plot(RDF_IMC_Cl_THF_h.results.bins,[x/4 for x in RDF_IMC_Cl_THF_h.results.rdf])
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_cl...IMC_mg","IMC_mg...THF_o","IMC_Cl...THF_h"],prop={'size': 14})
print(IMC_Mg)
#plt.xlim(limx)
#plt.ylim(limy2)

        
#%%PLOT IMC RDF SINGLES
fig=plt.figure(figsize=(10,5),dpi=150)
rows=(len(u_INT_sigle)//3+1)
limxs=[1,8]
limys=[0,0.035]
fig=plt.figure(figsize=(15,10),dpi=150)
dists=[5.94778,5.57157,5.32892,5.09535,4.93143,4.67588,4.37779]
posz=[x + 5.14 for x in dists]
lgnds=["IMC_ch...GCL_ho","IMC_mg...GCL_oh",
       "IMC_cl...GCL_ho","IMC_ch...GCL_c",
       "IMC_cl...CHL_ho",
       "IMC_mg...CHL_oh","IMC_mg...THF_o","IMC_mg...Clm"]
for i in range(len(u_INT_sigle)):
    
    ax=plt.subplot(rows,3,i+1)
    plt.plot(RDF_IMC1i_ch_GCL_ho_res[i].results.bins,[x/3 for x in RDF_IMC1i_ch_GCL_ho_res[i].results.rdf])
    plt.plot(RDF_IMC1i_mg_GCL_oh_res[i].results.bins,[x/3 for x in RDF_IMC1i_mg_GCL_oh_res[i].results.rdf])
    plt.plot(RDF_IMC1i_cl_GCL_ho_res[i].results.bins,[x/3 for x in RDF_IMC1i_cl_GCL_ho_res[i].results.rdf])
    plt.plot(RDF_IMC1i_ch_GCL_c_res[i].results.bins,[x/3 for x in RDF_IMC1i_ch_GCL_c_res[i].results.rdf])
    plt.plot(RDF_IMC1i_cl_CHL_ho_res[i].results.bins,RDF_IMC1i_cl_CHL_ho_res[i].results.rdf)
    plt.plot(RDF_IMC1i_mg_CHL_oh_res[i].results.bins,RDF_IMC1i_mg_CHL_oh_res[i].results.rdf)
    plt.plot(RDF_IMC1i_mg_THF_o_res[i].results.bins,RDF_IMC1i_mg_CHL_oh_res[i].results.rdf)
    plt.plot(RDF_IMC1i_mg_Clm_res[i].results.bins,RDF_IMC1i_mg_Clm_res[i].results.rdf)
    ax.set_xlim(limxs)
    ax.set_ylim(limys) 
    ax.set_xlabel(r'Distance$\ /\ \AA$')
    ax.set_ylabel(r'$g(r)$')
    ax.text(x=0, y=1.05, s='Z coordinate {:.3f} nm'.format(posz[i]) , fontsize=12, transform=ax.transAxes)
    ax.legend(lgnds)    

#plt.legend(lgnds)
plt.tight_layout()

#%%PLOT IMC RDF SINGLES one RDF on each plot

rows=(len(lgnds)//3+1)
limxs=[0,8]
limys=[0,0.035]


fig=plt.figure(figsize=(10,7),dpi=150)
ax=plt.subplot(rows,3,1)
for i in range(len(u_INT_sigle)):
        plt.plot(RDF_IMC1i_ch_GCL_ho_res[i].results.bins,[x/3 for x in RDF_IMC1i_ch_GCL_ho_res[i].results.rdf])
        plt.legend(["{:.3f}".format(x) for x in dists])
        plt.title(lgnds[0])
        
ax=plt.subplot(rows,3,2)
for i in range(len(u_INT_sigle)):
       plt.plot(RDF_IMC1i_mg_GCL_oh_res[i].results.bins,[x/3 for x in RDF_IMC1i_mg_GCL_oh_res[i].results.rdf])
       plt.title(lgnds[1])
ax=plt.subplot(rows,3,3)
for i in range(len(u_INT_sigle)):
        plt.plot(RDF_IMC1i_cl_GCL_ho_res[i].results.bins,[x/3 for x in RDF_IMC1i_cl_GCL_ho_res[i].results.rdf])
        plt.title(lgnds[2])

ax=plt.subplot(rows,3,4)
for i in range(len(u_INT_sigle)):
        plt.plot(RDF_IMC1i_ch_GCL_c_res[i].results.bins,RDF_IMC1i_ch_GCL_c_res[i].results.rdf)
        plt.title(lgnds[3])

ax=plt.subplot(rows,3,5)
for i in range(len(u_INT_sigle)):
        plt.plot(RDF_IMC1i_cl_CHL_ho_res[i].results.bins,RDF_IMC1i_cl_CHL_ho_res[i].results.rdf)
        plt.title(lgnds[4])
ax=plt.subplot(rows,3,6)
for i in range(len(u_INT_sigle)):
        plt.plot(RDF_IMC1i_mg_CHL_oh_res[i].results.bins,RDF_IMC1i_mg_CHL_oh_res[i].results.rdf)
        plt.title(lgnds[5])



plt.tight_layout()

#%%PLOT IMC RDF SINGLES one max of each RDF on each plot
dists=[5.94778,5.57157,5.32892,5.09535,4.93143,4.67588,4.37779]

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
xlab=('$z\ /\ nm$')

rows=(len(lgnds)//3+1)
limxs=[1,8]
limys=[0,0.035]
fig=plt.figure(figsize=(10,7),dpi=150)
ax=plt.subplot(rows,3,1)
plt.plot(posz,[max(x/3 for x in RDF_IMC1i_ch_GCL_ho_res[i].results.rdf) 
        for i in range(len(u_INT_sigle))],color=colors[0])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')
plt.title(lgnds[0])

ax=plt.subplot(rows,3,2)
plt.plot(posz,[max(x/3 for x in RDF_IMC1i_mg_GCL_oh_res[i].results.rdf)
       for i in range(len(u_INT_sigle))],color=colors[1])
plt.title(lgnds[1])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')

ax=plt.subplot(rows,3,3)
plt.plot(posz,[max(x/3 for x in RDF_IMC1i_cl_GCL_ho_res[i].results.rdf)
        for i in range(len(u_INT_sigle))],color=colors[2])
plt.title(lgnds[2])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')

x=plt.subplot(rows,3,4)
plt.plot(posz,[max(x/3 for x in RDF_IMC1i_ch_GCL_c_res[i].results.rdf) 
        for i in range(len(u_INT_sigle))],color=colors[3]) 
plt.title(lgnds[3])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')


ax=plt.subplot(rows,3,5)
plt.plot(posz,[max(RDF_IMC1i_cl_CHL_ho_res[i].results.rdf) 
        for i in range(len(u_INT_sigle))],color=colors[4]) 
plt.title(lgnds[4])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')

ax=plt.subplot(rows,3,6)
plt.plot(posz,[max(RDF_IMC1i_mg_CHL_oh_res[i].results.rdf) 
        for i in range(len(u_INT_sigle))],color=colors[5])
plt.title(lgnds[5])
plt.xlabel(xlab)
plt.ylabel(r'$g(r)$')
plt.tight_layout()



#%%PLOT CN SINGLES vs z pos
dists=[5.94778,5.57157,5.32892,5.09535,4.93143,4.67588,4.37779]
lgnds2=["IMC_ch...GCL_ho","IMC_mg...GCL_oh",
       "IMC_cl...GCL_ho","IMC_mg...CHL_oh",
       "IMC1i_cl...CHL_ho",
       "IMC_mg...THF_o"]
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
xlab=('$z\ /\ nm$')
ylab=('$CN$')

limy3=[0,2.8]

rows=(len(lgnds)//3+1)
limxs=[1,8]
limys=[0,0.035]
fig=plt.figure(figsize=(10,7),dpi=150)
ax=plt.subplot(rows,3,1)
plt.plot(posz,CN_IMC1i_ch_GCL_ho,color=colors[0])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.title(lgnds2[0])
plt.ylim(limy3)

ax=plt.subplot(rows,3,2)
plt.plot(posz,CN_IMC_mg_GCL_oh,color=colors[1])
plt.title(lgnds2[1])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(limy3)

ax=plt.subplot(rows,3,3)
plt.plot(posz,CN_IMC1i_cl_GCL_ho,color=colors[2])
plt.title(lgnds2[2])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(limy3)

x=plt.subplot(rows,3,4)
plt.plot(posz,CN_IMC_mg_CHL_oh,color=colors[3]) 
plt.title(lgnds2[3])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(limy3)


ax=plt.subplot(rows,3,5)
plt.plot(posz,CN_IMC1i_cl_CHL_ho,color=colors[4]) 
plt.title(lgnds2[4])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(limy3)

ax=plt.subplot(rows,3,6)
plt.plot(posz,CN_IMC_mg_THF_o,color=colors[5])
plt.title(lgnds2[5])
plt.xlabel(xlab)
plt.ylabel(ylab)
plt.ylim(limy3)
plt.tight_layout()
#%%
tot=np.trapz(RDF2.rdf,RDF2.bins)
#print(max(RDF1.rdf))
rdf_norm=[x/tot for x in RDF1.rdf]
plt.plot(RDF1.bins,RDF1.rdf/max(RDF1.rdf))
plt.plot(RDF2.bins,RDF2.rdf/max(RDF2.rdf))
plt.plot(RDF3.bins,RDF3.rdf/max(RDF3.rdf))
plt.ylabel("RDF")
#plt.subplot(3,1,2)
#plt.plot(RDF2.bins,RDF2.rdf)
#plt.ylabel("Mg1-Cl RDF")
#plt.subplot(3,1,3)
#plt.plot(RDF3.bins,RDF3.rdf)
#plt.plot(RDF1.bins, RDF1.rdf,RDF2.bins, RDF2.rdf)
plt.xlabel("Distance (Å)")
#plt.ylabel("Mg2-Cl RDF")
plt.tight_layout()
lgnd_dist=["Mg-Cl","Mg-Mg","Cl-Cl"]
plt.legend(lgnd_dist)
#plt.show()




#%%
# Calculate the RDF using the "AtomNeighborSearch" method from the MDAnalysis library
ans1=mda.analysis.distances.distance_array(Mg.positions,Cl.positions)
#ans1 = mda.analysis.distances.AtomNeighborSearch(both,5)
print(ans1)
hist, edges = np.histogram(ans1, bins=bins)
rdf = hist / (4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3) * len(Mg) * len(Cl))

# Plot the RDF
plt.plot(edges[:-1], rdf)
plt.xlabel("Distance (A)")
plt.ylabel("RDF")
plt.show()


# %%
r = 3.5

# Define the size of the radial bins
bin_size = 0.1

# Create an array for the radial bins
bins = np.arange(0, u.dimensions[0]/2, bin_size)

# Calculate the RDF between the central atom and all other atoms
ans = mda.analysis.distances.distance_array(atom_A.positions, atoms_other.positions)
hist, edges = np.histogram(ans, bins=bins)
rdf = hist / (4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3) * len(atom_A) * len(atoms_other))

# Integrate the area under the RDF curve up to distance r to get the coordination number
coord_num = np.trapz(rdf[edges[:-1] <= r], edges[:-1][edges[:-1] <= r])

# Print the coordination number
print("Coordination number up to a distance of {0} nm: {1}".format(r, coord_num))
# %%PLOT RDF COORD FROM GMX
RDFS={"Mg":{},"Cl":{}}
COORD={"Mg":{},"Cl":{}}

for traj in paths_ALL_mix:
        dir=os.path.dirname(traj)
        window=os.path.basename(dir)
        pairsMg=[]
        pairsCl=[]
        RDFS['Mg'][window]={}
        RDFS['Cl'][window]={}
        COORD['Mg'][window]={}
        COORD['Cl'][window]={}
        #RDFS
        with open(os.path.join(dir,'RDF_Mg.xvg'),'r') as ifile:
                lines=ifile.readlines()
                for line in lines:
                        if re.findall('s\d',line) :
                                pairsMg.append(line.split()[-1])
                x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
                RDFS['dist']=x
                for i,t in enumerate(pairsMg):
                        RDFS['Mg'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        with open(os.path.join(dir,'RDF_Cl.xvg'),'r') as ifile:
                lines=ifile.readlines()
                for line in lines:
                        if re.findall('s\d',line) :
                                pairsCl.append(line.split()[-1])
                for i,t in enumerate(pairsCl):
                        RDFS['Cl'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        #COORDINATION NUMBER                
        with open(os.path.join(dir,'coord_Mg.xvg'),'r') as ifile:
                lines=ifile.readlines()
                x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
                COORD['dist']=x
                for i,t in enumerate(pairsMg):
                        COORD['Mg'][window][t]=[float(x.split()[i+1]) for x in lines if "#" not in x if "@" not in x]
        with open(os.path.join(dir,'coord_Cl.xvg'),'r') as ifile:
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
        fig=plt.figure(figsize=(30,35),dpi=150)
        for k,i in enumerate(RDFS['Mg']):
                Coord=True
                try:
                        y_min=find_peaks([-x for x in RDFS['Mg'][i][j]],width=1)
                        x_min=RDFS['dist'][y_min[0][0]]
                except:
                        Coord=False
                ax = plt.subplot(5,6,k+1)
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
        plt.suptitle(j)
        plt.show()
#print(COORD['Mg']['window1'][pairs[0]][RDFS['dist'].index(x_min)])
#%% CL RDFS
for i,item in enumerate(pairsCl):
        clrs[item]=ref[i]
for j in pairsCl:
        clr=clrs[j]
        fig=plt.figure(figsize=(30,35),dpi=150)
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
                        c.text(x_min,COORD['Cl'][i][j][RDFS['dist'].index(x_min)]+(COORD['Cl'][i][j][RDFS['dist'].index(x_min)]*0.5),"{:.3f}".format(COORD['Cl'][i][j][RDFS['dist'].index(x_min)]) )
        plt.suptitle(j)
        plt.show()
                                
#%%
CN={'Mg':{},'Cl':{}}
for k in RDFS:
        if k == 'dist':
                continue
        else:
                for i in RDFS[k]:
                        CN[k][i]={}
                        for j in RDFS[k][i]:
                                try:
                                        y_min=find_peaks([-x for x in RDFS[k][i][j]],width=1)
                                        x_min=RDFS['dist'][y_min[0][0]]
                                        CN[k][i][j]=COORD[k][i][j][RDFS['dist'].index(x_min)]
                                except:
                                        CN[k][i][j]=0
# %%
#%%
#pos_umbrella=[9.435,9.178,8.888,8.65,8.372,8.11,7.85,7.572,7.307,7.038,6.771,
#6.507,6.238,5.973,5.711,5.448,5.168,4.9,4.647,4.391,4.105,3.839,3.583,3.312,
#3.01,2.789,2.518,2.266,1.981,1.716]
pos_umbrella=[]
for i in paths_ALL_mix:
        with open( ''.join(os.path.dirname(i)+'/umbrella.dump'),'r') as ifile:
                lines=ifile.readlines()
                pos=[x.split()[-1] for x in lines if " init " in x]
                pos_umbrella.append(float(pos[0]))
print(pos_umbrella)
                
#%%
posz_all=[x+5.14 for x in pos_umbrella]
limy_all=[0,180]
limx_all=[6,1]
axis=["X","Y","Z"]

fig=plt.figure(figsize=(10,10),dpi=150)

for i in range(3):
        plt.subplot(5,1,i+1)
        plt.errorbar(posz_all,[x[i] for x in ang_avg],[x[i] for x in std_ang_avg],linestyle=None)
        plt.plot(posz_all,[x[i] for x in ang_avg])
        plt.axvline(x=10.6, color = '#989898')
        plt.axvline(x=8.4, color = '#989898')
        plt.axhline(y=90, color = '#989898')
        plt.ylim(limy_all)
        plt.xlim(limx_all)
        plt.ylabel("{} Angle / deg".format(axis[i]))
plt.subplot(5,1,4)
plt.plot(data[0][0],data[0][1])
plt.xlim(limx_all)
plt.subplot(5,1,5)
plt.plot(GCL_dens_res.results.z.hist_bin_edges[:-1]/10 ,GCL_dens_res.results.z.mass_density)
plt.plot(CHL_dens_res.results.z.hist_bin_edges[:-1]/10 ,CHL_dens_res.results.z.mass_density)
plt.plot(THF_dens_res.results.z.hist_bin_edges[:-1]/10 ,THF_dens_res.results.z.mass_density)
plt.plot(Clm_dens_res.results.z.hist_bin_edges[:-1]/10 ,Clm_dens_res.results.z.mass_density)
plt.axvline(x=10.6, color = '#989898')
plt.axvline(x=8.4, color = '#989898')
plt.xlim(limx_all)
plt.legend(["GCL","CHL","THF","Clm"])
plt.xlabel('$z\ /\ nm$')
plt.show()