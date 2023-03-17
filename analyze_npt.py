#%%
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
from MDAnalysis.analysis import *
from MDAnalysis.analysis.bat import BAT
import glob

w_THF=range(1,15)
w_INT=range(15,22)
w_DES=range(22,27)
root="/run/media/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
mols="ACT"
#paths_THF=[root+mols+"/umbrella_30_200/window{}/umbrella.xtc".format(i) for i in w_THF]
#paths_INT=[root+mols+"/umbrella_30_200/window{}/umbrella.xtc".format(i) for i in w_INT]
#paths_DES=[root+mols+"/umbrella_30_200/window{}/umbrella.xtc".format(i) for i in w_DES]
#print(paths_THF)

# Load the LAMMPS trajectory
#for i in paths:
u = mda.Universe(root+"UMBRELLA_BIG_"+mols+"_MOD.parm7",root+mols+"/npt_iso_C_rescale_whole.xtc")

#%%

# Select the atoms to be used in the calculation
GCL_HO = u.select_atoms('resname GCL and type ho')
print(R1)
#%%
Mg2 = u.select_atoms('index 6')
Cl1 = u.select_atoms('index 3')
Cl2 = u.select_atoms('index 7')
Mg  = u.select_atoms('type 4')
Cl  = u.select_atoms('type 5')
both=u.select_atoms('type 5 or type 4')
print(len(u.trajectory))
#%%

A1 = u.select_atoms('index 3 2 7')
A2 = u.select_atoms('index 3 6 7')
angles1=[mda.core.topologyobjects.Angle([3,2,7],u).value() for ts in u.trajectory]
angles2=[mda.core.topologyobjects.Angle([3,6,7],u).value() for ts in u.trajectory]
angles3=[mda.core.topologyobjects.Angle([0,2,7],u).value() for ts in u.trajectory]
angles4=[mda.core.topologyobjects.Angle([0,2,3],u).value() for ts in u.trajectory]

lgnd=["Cl-Mg1-Cl","Cl-Mg2-Cl","C1-Mg1-Cl1","C1-Mg1-Cl2"]
#ang=mda.core.topologyobjects.Angle([3,2,7],u)
print(ang.value())
#print(len(u.trajectory))
print(angles1[1],angles1[2])
# Define the size of the radial bins
bin_size = 0.05
bins = 180
hist, edges = np.histogram(angles1, bins=bins)
#a1_prob = hist / (4/3 * np.pi * (edges[1:]**3 - edges[:-1]**3))
plt.plot(edges[:-1], hist)
hist, edges = np.histogram(angles2, bins=bins)
plt.plot(edges[:-1], hist)
hist, edges = np.histogram(angles3, bins=bins)
plt.plot(edges[:-1], hist)
hist, edges = np.histogram(angles4, bins=bins)
plt.plot(edges[:-1], hist)
# Create an array for the radial bins
plt.legend(lgnd)
plt.xlabel("Amplitude (deg)")

#%%
RDF1=mda.analysis.rdf.InterRDF(Mg, Cl, nbins=bins, range=(1.0, 5.0), exclusion_block=None)
RDF2=mda.analysis.rdf.InterRDF(Mg1, Mg2, nbins=bins, range=(1.0, 5.0), exclusion_block=None)
RDF3=mda.analysis.rdf.InterRDF(Cl1, Cl2, nbins=bins, range=(1.0, 5.0), exclusion_block=None)
RDF1.run()
RDF2.run()
RDF3.run()


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
