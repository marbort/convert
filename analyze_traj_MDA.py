#%%
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
from MDAnalysis.analysis import *
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.density import DensityAnalysis
from MDAnalysis.analysis.lineardensity import LinearDensity
import glob

#%%
def coord_number(x,y,dist_min,dist_max,mols):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*np.trapz(prod,val_x)
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
    
    return(bins,sum_tmp,sum_norm,avg_vol)


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
    dens=len(com_dist[0])/avg_vol
    sum_tmp=np.zeros(len(bins)-1)
    for i in hists:
        sum_tmp=sum_tmp+i[0] 
    sum_norm=(sum_tmp/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)*(len(gA_coms))))
    
    return(bins,sum_tmp,sum_norm,avg_vol,dens)
        
def RDF_density(gA,gB,Universe,binsize,max,density):
    groupA=Universe.select_atoms(gA)
    groupB=Universe.select_atoms(gB)
    gA_coms=[groupA.positions for ts in Universe.trajectory]
    gB_coms=[groupB.positions for ts in Universe.trajectory]
    boxes=[Universe.dimensions for ts in Universe.trajectory]
    if groupA == groupB:
        com_arr  = [distances.self_distance_array(x,boxes[i]) for i,x in enumerate(gA_coms)]
        com_dist = [np.tile(x,2) for x in com_arr]
    else:
        com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i]) for i,x in enumerate(gA_coms)]
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
    
    return(bins,sum_tmp,sum_norm,avg_vol,dens,len(groupB))
        
    
    


root="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/"
coms=glob.glob(root+"*_com.xvg")
#paths=glob.glob(root+"Iteration*/*/*.dcd")
#print(paths)
chars=["#","@"]
coms_dict={}
# Load the LAMMPS trajectory
#for i in paths:
u = mda.Universe(root+"interface_BIG_MOD_BOTH.parm7",root+"npt_iso_C_rescale_whole.xtc",root+"npt_iso_C_rescale2.xtc")
print(len(u.trajectory))
#%%
mols=["CHL Clm","GCL","THF","IPM","IMC","MG Cl-","ACT","Clm"]
labs={"CHL":"CholCl","GCL":"Glycerol","THF":"THF","IPM":"iPr2Mg",
"IMC":"iPrMgCl","MG":"MgCl-","ACT":"Acetone"}
clrs=['#F22727','#12ab07','#072dab','#4b68c9','#8393c7','#ffa6a7','#ff8533']
clrs_dict={"CHL":'#F22727',"GCL":'#12ab07',"THF":'#072dab',
           "IPM":'#4b68c9',"IMC":'#8393c7',"MG":'#ffa6a7',"ACT":'#ff8533',"Clm":'#8A2BE2'}
# Select the atoms to be used in the calculation
CHL = u.select_atoms('resname CHL')
CHL_com=CHL.center_of_mass()
GCL = u.select_atoms('resname GCL')
GCL_com=GCL.center_of_mass()
THF = u.select_atoms('resname THF')
ACT = u.select_atoms('resname ACT')
IMC = u.select_atoms('resname IMC')
IPM = u.select_atoms('resname IPM')
MGCL = u.select_atoms('resname MG Cl-')
ACT = u.select_atoms('resname ACT')
CHL_cl=CHL.select_atoms('type Clm')
O = u.select_atoms('name O*')

print(len(GCL)/14)
print(len(O))
us=[u.select_atoms('resname {}'.format(i)) for i in mols]
DES = u.select_atoms('resname CHL GCL Clm')

#RDF_com=mda.analysis.rdf.InterRDF(CHL_com, GCL_com, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)

#%%
for file in coms:
    with open(file,'r') as comfile:
        lines=comfile.readlines()
        name=os.path.basename(comfile.name).replace("_com.xvg","")
        t=[float(x.split()[0])/1000 for x in lines if x[0] not in chars ]
        x=[float(x.split()[1]) for x in lines if x[0] not in chars ]
        y=[float(x.split()[2]) for x in lines if x[0] not in chars ]
        z=[float(x.split()[3]) for x in lines if x[0] not in chars ]
        coms_dict.update({name:[t,x,y,z]})
#print(coms_dict)

fig=plt.figure(dpi=150)
lgnd=[]
for i in coms_dict:
    plt.plot(coms_dict[i][3],coms_dict[i][0],color=clrs_dict[i])
    lgnd.append(labs[i])
plt.legend(lgnd)
plt.axvline(x=10.6, color = '#606060')
plt.axvline(x=8.37, color = '#606060')       

plt.axvline(x=20.2, color = '#606060')
plt.axvline(x=1.15, color = '#606060')
plt.xlabel(r'$z\ /\ nm$')
plt.ylabel(r'$t\ /\ ns$')
plt.title("Center of Mass Position Over Time")


#print(len(u.trajectory))
#%%

densities=[LinearDensity(x,binsize=1) for x in us]
dens_results=[x.run() for x in densities]
DESd=LinearDensity(DES,binsize=1)
DESd.run()

#%%
#print(THFd.results.z.dim)
#plt.plot(dens_results[0].results.z.hist_bin_edges[:-1],dens_results[0].results.z.mass_density)
fig=plt.figure(dpi=150)
for i,x in enumerate(dens_results):
    plt.plot(x.results.z.hist_bin_edges[:-1]/10 ,x.results.z.mass_density,color=clrs[i])
plt.plot(DESd.results.z.hist_bin_edges[:-1]/10 ,DESd.results.z.mass_density,'--k')


lgnd=[labs[i] for i in labs]+["DES"]
plt.xlabel(r'$z\ /\ nm$')
plt.ylabel(r'$\rho\ /\ g cm^{-3}$')
plt.legend(lgnd,loc='lower right')

plt.axvline(x=10.6, color = '#606060')
plt.axvline(x=8.37, color = '#606060')       

plt.axvline(x=20.2, color = '#606060')
plt.axvline(x=1.15, color = '#606060')
plt.title("Density")
#%%
for i,x in enumerate(dens_results):
    plt.plot(x.results.z.hist_bin_edges[:-1]/10 ,x.results.z.mass_density,color=clrs[i])

limyd=[0,0.2]
limxd=[7.5,12.5]
lgnd=mols
plt.xlabel(r'$z\ /\ nm$')
plt.ylabel(r'$\rho\ /\ g cm^{-3}$')
plt.legend(mols)

plt.axvline(x=10.6, color = '#606060')
plt.axvline(x=8.37, color = '#606060')       

plt.axvline(x=20.2, color = '#606060')
plt.axvline(x=1.15, color = '#606060')
plt.title("Density")
plt.ylim(limyd)
plt.xlim(limxd)
#plt.legend(lgnd)



#%% Plot RDFs
bins=80
normp="density"
CHL_ho=CHL.select_atoms('type ho')
CHL_oh=CHL.select_atoms('type oh')
CHL_c=CHL.select_atoms('type c3')
CHL_n=CHL.select_atoms('type n4')

GCL_oh=GCL.select_atoms('type oh')
GCL_oh1=GCL.select_atoms('name O1')
GCL_ho=GCL.select_atoms('type ho')
GCL_ho6=GCL.select_atoms('name H6')
GCL_c=GCL.select_atoms('type c3')
THF_o =THF.select_atoms('type os')
THF_c = THF.select_atoms('type c3')
THF_h = THF.select_atoms('name H1 H2 H3 H4 H5 H6 H7 H8')
IMC_Cl=IMC.select_atoms('name Cl1')
IMC_Mg=IMC.select_atoms('name Mg1')
IMC_c=IMC.select_atoms('name C1')
ACT_o=ACT.select_atoms('type o')
ACT_c=ACT.select_atoms('name C1')
IPM_Mg=IPM.select_atoms('name Mg1')
IPM_c=IPM.select_atoms('name C1 C4')
MGCL_Mg=MGCL.select_atoms('name MG')
MGCL_Cl=MGCL.select_atoms('name Cl-')

IMC_coms=[IMC.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]
GCL_coms=[GCL.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]
CHL_coms=[CHL.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]
ACT_coms=[ACT.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]
IPM_coms=[IPM.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]
THF_coms=[THF.center_of_mass(pbc=True,compound='residues') for ts in u.trajectory]

#%%
print(len(GCL.center_of_mass(pbc=True,compound='residues')))
print(len(IMC.center_of_mass(pbc=True,compound='residues')))
#%%
print(GCL_coms[0])
#print(CHL_coms)

IMC_3378=u.select_atoms('resid 3356')
IMC_3378_Mg=IMC_3378.select_atoms('name Mg1')
#print(CHL_cl)


#%% COMS
COMS={mols[i]:[x.center_of_mass for ts in u.trajectory] for i,x in enumerate(us)}
#%%
print(COMS['GCL'][0])
#%%
#RDF_CHL
RDF_CHL_ho_GCL_oh=mda.analysis.rdf.InterRDF(CHL_ho, GCL_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL_oh_GCL_ho=mda.analysis.rdf.InterRDF(CHL_oh, GCL_ho, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL_oh_THF_o=mda.analysis.rdf.InterRDF(CHL_oh, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL_c_THF_c=mda.analysis.rdf.InterRDF(CHL_c, THF_c, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL_ho_ACT_o=mda.analysis.rdf.InterRDF(CHL_ho, ACT_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_CHL_oh_CHL_ho=mda.analysis.rdf.InterRDF(CHL_ho, CHL_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=(1,1),norm=normp)
RDF_CHL_n_CHL_n=mda.analysis.rdf.InterRDF(CHL_n, CHL_n, nbins=bins, range=(1.0, 8.0), exclusion_block=(1,1),norm=normp)


RDF_CHL_ho_GCL_oh.run()
RDF_CHL_oh_GCL_ho.run()
RDF_CHL_oh_THF_o.run()
RDF_CHL_c_THF_c.run()
RDF_CHL_ho_ACT_o.run()
RDF_CHL_oh_CHL_ho.run()
RDF_CHL_n_CHL_n.run()

#%%
RDF_CHL_N_CHL_N=RDF_density("resname CHL and type n4","resname CHL and type n4",u,0.1,15,0.002255665)
#%%
r=8
print(RDF_CHL_N_CHL_N[-1])
print(RDF_CHL_N_CHL_N[-2])
plt.plot(RDF_CHL_N_CHL_N[0][:-1],RDF_CHL_N_CHL_N[2])
cn=np.trapz(RDF_CHL_N_CHL_N[2][RDF_CHL_N_CHL_N[0][:-1] <= r],RDF_CHL_N_CHL_N[0][:-1][RDF_CHL_N_CHL_N[0][:-1] <= r])
print(cn)
#%% CN CHL_ho_CHL_oh
CN_CHL_ho_CHL_oh =[]
n_res=834
rdfa=RDF_CHL_oh_CHL_ho
CN_CHL_ho_CHL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,0,2.5,n_res))
CN_CHL_ho_CHL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,2.6,6,n_res))
CN_CHL_ho_CHL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,6,8,n_res))
print(CN_CHL_ho_CHL_oh)

#%% CN CHL_ho_GCL_OH
CN_CHL_ho_GCL_oh =[]
n_res=834
rdfa=RDF_CHL_oh_GCL_ho
CN_CHL_ho_GCL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,0,2.5,n_res))
CN_CHL_ho_GCL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,2.6,6,n_res))
CN_CHL_ho_GCL_oh.append(coord_number(rdfa.results.bins,rdfa.results.rdf,6,8,n_res))
print(CN_CHL_ho_GCL_oh)

#%%RDF GCL
RDF_GCL_ho_THF_o=mda.analysis.rdf.InterRDF(GCL_ho, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_GCL_c_GCL_c=mda.analysis.rdf.InterRDF(GCL_c, GCL_c, nbins=bins, range=(1.0, 8.0), exclusion_block=(3,3),norm=normp)
RDF_GCL_ho_GCL_oh=mda.analysis.rdf.InterRDF(GCL_ho, GCL_oh, nbins=bins, range=(1.0, 8.0), exclusion_block=(3,3),norm=normp)
RDF_GCL_ho6_GCL_oh1=mda.analysis.rdf.InterRDF(GCL_ho6, GCL_oh1, nbins=bins, range=(1.0, 8.0), exclusion_block=(1,1),norm=normp)
RDF_GCL_ho_THF_o.run()
RDF_GCL_ho_GCL_oh.run()
RDF_GCL_c_GCL_c.run()
RDF_GCL_ho6_GCL_oh1.run()
#%% CN_GCL_ho_GCL_OH
CN_GCL_ho_GCL_OH =[]
CN_GCL_ho_GCL_OH.append(coord_number(RDF_GCL_ho_GCL_oh.results.bins,RDF_GCL_ho_GCL_oh.results.rdf,0,2.5,1668))
CN_GCL_ho_GCL_OH.append(coord_number(RDF_GCL_ho_GCL_oh.results.bins,RDF_GCL_ho_GCL_oh.results.rdf,2.6,4,1668))
CN_GCL_ho_GCL_OH.append(coord_number(RDF_GCL_ho_GCL_oh.results.bins,RDF_GCL_ho_GCL_oh.results.rdf,4.1,6.5,1668))
print(CN_GCL_ho_GCL_OH)
#%%

RDF_GCL_ho_GCL_oh_1st_peakx=[x for x in RDF_GCL_ho_GCL_oh.results.bins if 2.5<=x<=4]
RDF_GCL_ho_GCL_oh_1st_peaky=[RDF_GCL_ho_GCL_oh.results.rdf[i] for i,x in enumerate(RDF_GCL_ho_GCL_oh.results.bins) if 2.5<=x<=4]
y2=[x**2*RDF_GCL_ho_GCL_oh_1st_peaky[i] for i,x in enumerate(RDF_GCL_ho_GCL_oh_1st_peakx)]
tot=4*np.pi*np.trapz(y2,RDF_GCL_ho_GCL_oh_1st_peakx)
print(tot/1668)

#%%
RDF_IMC_Cl_THF_h=mda.analysis.rdf.InterRDF(IMC_Mg, THF_h, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_Mg_THF_o=mda.analysis.rdf.InterRDF(IMC_Mg, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_Cl_IMC_Mg=mda.analysis.rdf.InterRDF(IMC_Cl, IMC_Mg, nbins=bins, range=(1.0, 8.0), exclusion_block=(2,2),norm=normp)
RDF_IMC_Mg_ACT_o=mda.analysis.rdf.InterRDF(IMC_Mg, ACT_o ,nbins=bins, range=(1.0, 8.0), exclusion_block=(2,2),norm=normp)
RDF_IMC_c_ACT_c=mda.analysis.rdf.InterRDF(IMC_c, ACT_c ,nbins=bins, range=(1.0, 8.0), exclusion_block=(2,2),norm=normp)
RDF_IMC_Mg_MGCL_Cl=mda.analysis.rdf.InterRDF(IMC_Mg, MGCL_Cl, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_Cl_MGCL_Mg=mda.analysis.rdf.InterRDF(IMC_Cl, MGCL_Mg, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_3378_Mg_THF_o=mda.analysis.rdf.InterRDF(IMC_3378_Mg, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)

RDF_IMC_Mg_O=mda.analysis.rdf.InterRDF(IMC_Mg, O ,nbins=bins, range=(1.0, 8.0), exclusion_block=(2,2),norm=normp)

RDF_IMC_Mg_THF_o.run()
RDF_IMC_Cl_IMC_Mg.run()
RDF_IMC_Cl_THF_h.run()
RDF_IMC_Mg_ACT_o.run()
RDF_IMC_c_ACT_c.run()
RDF_IMC_Mg_MGCL_Cl.run()
RDF_IMC_Cl_MGCL_Mg.run()
RDF_IMC_Mg_O.run()

RDF_IMC_3378_Mg_THF_o.run()



RDF_IMC_3378_Mg_THF_o_1st_peakx=[x**2 for x in RDF_IMC_3378_Mg_THF_o.results.bins if x<=3]
RDF_IMC_3378_Mg_THF_o_1st_peaky=[RDF_IMC_3378_Mg_THF_o.results.rdf[i] for i,x in enumerate(RDF_IMC_3378_Mg_THF_o.results.bins) if x<=3]

#%%
CN_IMC_mg_THF_o =[]
CN_IMC_mg_THF_o.append(coord_number(RDF_IMC_Mg_THF_o.results.bins,RDF_IMC_Mg_THF_o.results.rdf,0,3,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_oresults.rdf,2.6,4,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_o.results.rdf,4.1,6.5,142))
print(CN_IMC_mg_THF_o)

CN_IMC_mg_ACT_o =[]
CN_IMC_mg_ACT_o.append(coord_number(RDF_IMC_Mg_ACT_o.results.bins,RDF_IMC_Mg_ACT_o.results.rdf,0,3,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_oresults.rdf,2.6,4,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_o.results.rdf,4.1,6.5,142))
print(CN_IMC_mg_ACT_o)

plt.plot(RDF_IMC_Mg_O.results.bins,RDF_IMC_Mg_O.results.rdf)
plt.show()
CN_IMC_mg_O=coord_number(RDF_IMC_Mg_O.results.bins,RDF_IMC_Mg_O.results.rdf,0,3,42)
print(CN_IMC_mg_O)


RDF_IMC_mg_THF_o_1st_peakx=[x for x in RDF_IMC_Mg_THF_o.results.bins if 2<x<3]
RDF_IMC_mg_THF_o_1st_peaky=[RDF_IMC_Mg_THF_o.results.rdf[i] for i,x in enumerate(RDF_IMC_Mg_THF_o.results.bins) if 2<x<3]
y2=[x**2*RDF_IMC_mg_THF_o_1st_peaky[i] for i,x in enumerate(RDF_IMC_mg_THF_o_1st_peakx)]
tot=4*np.pi*np.trapz(y2,RDF_IMC_mg_THF_o_1st_peakx)
#print(tot/42)

#%%
RDF_IMC_Cl_IPM_Mg=mda.analysis.rdf.InterRDF(IMC_Cl, IPM_Mg, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IMC_Cl_IPM_Mg.run()


#%% IPM RDF
RDF_IPM_mg_THF_o=mda.analysis.rdf.InterRDF(IPM_Mg, THF_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IPM_mg_ACT_o=mda.analysis.rdf.InterRDF(IPM_Mg, ACT_o, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IPM_c_ACT_c=mda.analysis.rdf.InterRDF(IPM_c, ACT_c, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IPM_mg_MGCL_Cl=mda.analysis.rdf.InterRDF(IPM_Mg, MGCL_Cl, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_IPM_mg_THF_o.run()
RDF_IPM_mg_ACT_o.run()
RDF_IPM_c_ACT_c.run()
RDF_IPM_mg_MGCL_Cl.run()
RDF_IPM_mg_THF_o_1st_peakx=[x for x in RDF_IPM_mg_THF_o.results.bins if x<=3]

#%%
CN_IPM_mg_THF_o =[]
CN_IPM_mg_THF_o.append(coord_number(RDF_IPM_mg_THF_o.results.bins,RDF_IPM_mg_THF_o.results.rdf,0,3,24))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_oresults.rdf,2.6,4,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_o.results.rdf,4.1,6.5,142))
print(CN_IPM_mg_THF_o)

CN_IPM_mg_ACT_o =[]
CN_IPM_mg_ACT_o.append(coord_number(RDF_IPM_mg_ACT_o.results.bins,RDF_IPM_mg_ACT_o.results.rdf,0,3,24))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_oresults.rdf,2.6,4,42))
#CN_IMC_mg_THF_o.append(coord_number(CN_IMC_mg_THF_o.results.bins,CN_IMC_mg_THF_o.results.rdf,4.1,6.5,142))
print(CN_IPM_mg_ACT_o)

#%%
RDF_IPM_mg_THF_o_1st_peaky=[RDF_IPM_mg_THF_o.results.rdf[i] for i,x in enumerate(RDF_IPM_mg_THF_o.results.bins) if x<=3]
y2=[x**2*RDF_IPM_mg_THF_o_1st_peaky[i] for i,x in enumerate(RDF_IPM_mg_THF_o_1st_peakx)]
tot=4*np.pi*np.trapz(y2,RDF_IPM_mg_THF_o_1st_peakx)
print(tot/24)

#%%RDF MGCL-MGCL

RDF_MGCL_mg_MGCL_cl=mda.analysis.rdf.InterRDF(MGCL_Mg, MGCL_Cl, nbins=bins, range=(1.0, 8.0), exclusion_block=None,norm=normp)
RDF_MGCL_mg_MGCL_cl.run()
RDF_MGCL_mg_MGCL_cl_site=mda.analysis.rdf.InterRDF_s(u,[[MGCL_Mg, MGCL_Cl]])
RDF_MGCL_mg_MGCL_cl_site.run()
cdf=RDF_MGCL_mg_MGCL_cl_site.get_cdf()

RDF_MGCL_mg_MGCL_cl_1st_peakx=[x for x in RDF_MGCL_mg_MGCL_cl.results.bins if 2<x<3]
RDF_MGCL_mg_MGCL_cl_1st_peaky=[RDF_MGCL_mg_MGCL_cl.results.rdf[i] for i,x in enumerate(RDF_MGCL_mg_MGCL_cl.results.bins) if 2<x<3]
y2=[x**2*RDF_MGCL_mg_MGCL_cl_1st_peaky[i] for i,x in enumerate(RDF_MGCL_mg_MGCL_cl_1st_peakx)]
tot=4*np.pi*np.trapz(y2,RDF_MGCL_mg_MGCL_cl_1st_peakx)
tot_all=np.trapz(RDF_MGCL_mg_MGCL_cl.results.rdf,RDF_MGCL_mg_MGCL_cl.results.bins)
print(tot,tot_all)
#%% LIMITS
limx=[1,8]
limy=[0,9]
limy2=[0,0.015]
limy3=[0.0,0.25]
#%% RDF CHL-GCL PLOT
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL_ho_GCL_oh.results.bins,[x/3 for x in RDF_CHL_ho_GCL_oh.results.rdf])
plt.plot(RDF_CHL_oh_GCL_ho.results.bins,[x/3 for x in RDF_CHL_oh_GCL_ho.results.rdf])
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_oh...GCL_ho","CHL_ho...GCL_oh"],prop={'size': 14})
#plt.xlim(limx)
plt.ylim(limy)
#RDF2.run()
#RDF3.run()


#%% RDF CHL-THF
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL_oh_THF_o.results.bins,RDF_CHL_oh_THF_o.results.rdf)
plt.plot(RDF_CHL_c_THF_c.results.bins,RDF_CHL_c_THF_c.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_oh...THF_o","CHL_c...THF_c"],prop={'size': 14})
plt.xlim(limx)
plt.ylim(limy)
#%% RDF CHL-ACT
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL_ho_ACT_o.results.bins,RDF_CHL_ho_ACT_o.results.rdf)
#plt.plot(RDF_CHL_oh_ACT_o.results.bins,RDF_CHL_oh_ACT_o.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_oh...ACT_o","CHL_c...THF_c"],prop={'size': 14})
plt.xlim(limx)
plt.ylim(limy)

#%% RDF CHL-CHL
fig=plt.figure(dpi=150)
plt.plot(RDF_CHL_n_CHL_n.results.bins,RDF_CHL_n_CHL_n.results.rdf)
plt.plot(RDF_CHL_oh_CHL_ho.results.bins,RDF_CHL_oh_CHL_ho.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["CHL_n...CHL_n","CHL_oh...CHL_ho"],prop={'size': 14})
plt.xlim(limx)
plt.ylim(limy)

#%% RDF GCL-GCL
fig=plt.figure(dpi=150)
div=[x/9 for x in RDF_GCL_ho_GCL_oh.results.rdf]
div2=[x/9 for x in RDF_GCL_c_GCL_c.results.rdf]
plt.plot(RDF_GCL_ho_GCL_oh.results.bins,div)
plt.plot(RDF_GCL_c_GCL_c.results.bins,div2)
plt.plot(RDF_GCL_ho_THF_o.results.bins,RDF_GCL_ho_THF_o.results.rdf)


plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["GCL_ho...GCL_oh","GCL_c...GCL_c","GCL_ho...THF_o"],prop={'size': 14})
#plt.xlim(limx)
plt.ylim(limy)

#%%PLOT IMC RDF
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC_Cl_IMC_Mg.results.bins,[x/42 for x in RDF_IMC_Cl_IMC_Mg.results.rdf])
plt.plot(RDF_IMC_Mg_THF_o.results.bins,[x/42 for x in RDF_IMC_Mg_THF_o.results.rdf])
plt.plot(RDF_IMC_Mg_ACT_o.results.bins,[x/42 for x in RDF_IMC_Mg_ACT_o.results.rdf])
plt.plot(RDF_IMC_c_ACT_c.results.bins,[x/42 for x in RDF_IMC_c_ACT_c.results.rdf])
plt.plot(RDF_IMC_Cl_THF_h.results.bins,[x/42 for x in RDF_IMC_Cl_THF_h.results.rdf])





plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_cl...IMC_mg","IMC_mg...THF_o","IMC_mg...ACT_o",
            "IMC_c...ACT_c","IMC_Cl...THF_h","IMC_cl...MGCL_mg","IMC_mg...MGCL...cl"],
            prop={'size': 14})
plt.xlim(limx)
plt.ylim(limy2)


#%%RDF IMC MGCL
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC_Cl_MGCL_Mg.results.bins,[x/42 for x in RDF_IMC_Cl_MGCL_Mg.results.rdf])
plt.plot(RDF_IMC_Mg_MGCL_Cl.results.bins,[x/42 for x in RDF_IMC_Mg_MGCL_Cl.results.rdf])
plt.plot(RDF_IMC_Mg_MGCL_Cl_site.re)

plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_cl...MGCL_mg","IMC_mg...MGCL...cl"],
            prop={'size': 14})
#plt.ylim(limy3)
plt.ylim(limy2)
#%% RDF IPM
fig=plt.figure(dpi=150)
plt.plot(RDF_IPM_mg_THF_o.results.bins,[x/24 for x in RDF_IPM_mg_THF_o.results.rdf],color='orange')
plt.plot(RDF_IPM_mg_ACT_o.results.bins,[x/24 for x in RDF_IPM_mg_ACT_o.results.rdf],color='green')
plt.plot(RDF_IPM_c_ACT_c.results.bins,[x/24 for x in RDF_IPM_c_ACT_c.results.rdf],color='red')
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IPM_mg...THF_o","IPM_mg...ACT_o","IPM_c...ACT_c"],prop={'size': 14})
plt.xlim(limx)
plt.ylim(limy2)

#%%RDF IPM MGCL
fig=plt.figure(dpi=150)
plt.plot(RDF_IPM_mg_MGCL_Cl.results.bins,[x/24 for x in RDF_IPM_mg_MGCL_Cl.results.rdf],color='orange')
plt.legend(["IPM_mg...MGCL_Cl"],prop={'size': 14})
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.ylim(limy2)
#%% RDF IPM-IMC
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC_Cl_IPM_Mg.results.bins,RDF_IMC_Cl_IPM_Mg.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_Cl...IPM_Mg"],prop={'size': 14})
plt.xlim(limx)
#plt.ylim(limy2)
#print(RDF_IMC_Cl_IPM_Mg.results.rdf)

#%% PLOT RDF MGCL MGCL
fig=plt.figure(dpi=150)
plt.plot(RDF_MGCL_mg_MGCL_cl.results.bins,[x/24 for x in RDF_MGCL_mg_MGCL_cl.results.rdf])
#plt.plot(RDF_MGCL_mg_MGCL_cl_site.results.bins,RDF_MGCL_mg_MGCL_cl_site.results.rdf[0][0,0])
#plt.plot(RDF_MGCL_mg_MGCL_cl_site.results.bins,RDF_MGCL_mg_MGCL_cl_site.results.cdf[0][0,0])
plt.legend(["MGCL_mg...MGCL_cl"],prop={'size': 14})
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.ylim(limy2)

#%% RDF IPM-IMC
fig=plt.figure(dpi=150)
plt.plot(RDF_IMC_Cl_IPM_Mg.results.bins,RDF_IMC_Cl_IPM_Mg.results.rdf)
plt.xlabel(r'Distance$\ /\ \AA$')
plt.ylabel(r'$g(r)$')
plt.legend(["IMC_Cl...IPM_Mg"],prop={'size': 14})
plt.xlim(limx)

#A1a=mda.analysis.bat.BAT(A1)
#A1a.run()
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
# Define the size of the radial bins
bin_size = 0.1

# Create an array for the radial bins
bins = np.arange(0, u.dimensions[0]/2, bin_size)
ans1=mda.analysis.distances.distance_array(IMC_Mg.positions,THF_o.positions)
#ans1 = mda.analysis.distances.AtomNeighborSearch(both,5)
print(ans1)
hist, edges = np.histogram(ans1, bins=bins)
rdf = hist / (4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3) * len(IMC_Mg) * len(THF_o))

# Plot the RDF
plt.plot(edges[:-1], rdf)
plt.xlabel("Distance (A)")
plt.ylabel("RDF")
plt.show()

#%%


# %%
r = 7.0

# Define the size of the radial bins
bin_size = 0.1

# Create an array for the radial bins
bins = np.arange(0, 15, bin_size)

# Calculate the RDF between the central atom and all other atoms


#%%
dists=[mda.analysis.distances.distance_array(x, CHL_coms[i]) for i,x in enumerate(GCL_coms)]

#%%
#ans=mda.analysis.distances.AtomNeighborSearch(GCL_coms)
#print(min(ans[-1]))
hists=[np.histogram(x, bins=bins) for x in dists]
#%%
ans = mda.analysis.distances.distance_array(GCL, CHL)
hist, edges = np.histogram(ans, bins=bins)
#%%
#print(hists[0][1])
#print(hists[0])
print(len(GCL_coms[0]))
print(len(hists[0][0]/(4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3) * len(GCL_coms)*len(CHL_coms))))
sum_tmp=np.zeros(149)
for i in hists:
    sum_tmp=sum_tmp+i[0]
sum_norm=(sum_tmp/(4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3) * len(GCL_coms[0])))
#print(sum_tmp/(4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3) * len(GCL_coms)*len(CHL_coms)))
    



#%%
#rdfa = hist / (4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3) * len(GCL_coms)*len(CHL_coms))
#print(edges[50])

#plt.plot(edges[:-1], rdfa)
plt.plot(hists[25][1][:-1], sum_norm)
plt.xlabel("Distance (A)")
plt.ylabel("RDF")
#plt.ylim([0,0.01])
plt.show()
#print(rdfa)
#Integrate the area under the RDF curve up to distance r to get the coordination number
coord_num = np.trapz(sum_norm[hists[0][1][:-1] <= r], hists[0][1][:-1][hists[0][1][:-1] <= r])

# Print the coordination number
#print(coord_num)
print("Coordination number up to a distance of {0} Å: {1}".format(r, coord_num))
# %%
with open('/run/media/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/GCL_CHL_rdf_vol.dat','r') as ifile:
    x=[]
    y=[]
    lines=ifile.readlines()
    for line in lines:
        if "#" in line:
            continue
        else:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))
x_tomin=[k  for k in x if k <= r]
y_tomin=[y[j] for j,k in enumerate(x) if k <= r ]

CN=np.trapz(y_tomin,x_tomin)
print(CN)
            

# %%
ress=["Clm","GCL","CHL","THF","ACT","IPM","IMC"]
rdfs=[]
bins_all=[]
names=[]
for i,x in enumerate(ress):
    for j in ress:
        if ress.index(j)<=i:  
            #print(x,j)        
            bs,rf,av=RDF_COM_volume("resname "+x,"resname "+j,u,0.1,15)
            rdfs.append(rf)
            bins_all.append(bs)
            names.append("{} {}".format(x,j))
            
for i,item in enumerate(rdfs):
    plt.plot(bins_all[i][:-1],item,color=clrs[i])
plt.legend(names)
# %%
#%%
RDF_IMC_Mg_Clm=RDF_POS_volume("resname IMC and name Mg1","type Clm",u,0.1,15)
#%%
RDF_GCL_C_GCL_C=RDF_POS_volume("resname GCL and name C2","resname GCL and name C2",u,0.1,15)
#%%
plt.plot(RDF_IMC_Mg_Clm[0][:-1],RDF_IMC_Mg_Clm[2])
# %%
