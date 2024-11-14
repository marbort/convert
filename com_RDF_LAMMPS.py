import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rdf import InterRDF_s
from MDAnalysis.analysis import *


class replace_with_COM:
    """Replace special atom `atomname` in each fragment with COM of the fragment."""
    def __init__(self, polymer, atomname):
        self.polymer = polymer
        self.com_atoms = polymer.select_atoms(f"name {atomname}")
        
        # sanity check
        assert self.get_com().shape == self.com_atoms.positions.shape
        
    def get_com(self):
        return self.polymer.center_of_mass(compound="residues")
    def get_cog(self):
        return self.polymer.center_of_geometry(compound="residues")
    
    def __call__(self, ts):
        self.com_atoms.positions = self.get_cog()
        return ts

class replace_with_COM_atoms:
    """Replace special atom `atomname` in each fragment with COM of the fragment."""
    def __init__(self, polymer, atomname):
        self.polymer = polymer
        self.com_atoms = polymer.select_atoms(f"name {atomname}")
        
        # sanity check
        #assert self.get_com().shape == self.com_atoms.positions.shape
        
    def get_com(self):
        return self.polymer.center_of_mass()
    
    def get_cog(self):
        return self.polymer.center_of_geometry()
    
    def __call__(self, ts):
        self.com_atoms.positions = self.get_cog()
        return ts



# Load your LAMMPS trajectory and topology files
u = mda.Universe('test.pdb', 'traj.dcd')
#polymer = ucom.select_atoms("resname AAA BBB CCC DDD")
#ucom.trajectory.add_transformations(replace_with_COM(polymer, "T1"))

# Define atom groups based on atom selections, e.g., by residue or atom types
group1 = u.select_atoms("bynum 1 9 13")  # replace 'RES1' with your group 1 selection
#group1 = u.select_atoms("bynum 1 8")  # replace 'RES1' with your group 1 selection
#group1 = u.select_atoms("bynum 1 ")  # replace 'RES1' with your group 1 selection
#group2 = u.select_atoms("(around 3 name Mg) and resname THF",updating=True)  # replace 'RES2' with your group 2 selection
group3 = u.select_atoms("resname THF")  # replace 'RES2' with your group 2 selection
u.trajectory.add_transformations(replace_with_COM(group3,"O"),replace_with_COM_atoms(group1,"C0"))

com_f1=u.select_atoms("name C0")
com_f2=u.select_atoms("name O and not around 3.5 name Mg", updating=True)
#com_f2=u.select_atoms("name O")








# Ensure the groups are non-empty
#if len(group1) == 0 or len(group2) == 0:
#    raise ValueError("One of the atom groups is empty. Check your selection criteria.")


# Calculate the RDF between the centers of mass of group1 and group2
rdf_com = mda.analysis.rdf.InterRDF(com_f1, com_f2, nbins=75, range=(0.0, 15.0), exclusion_block=(1,1))
rdf_com.run()

np.savetxt('rdf.txt',list(zip(np.transpose(rdf_com.results.bins),np.transpose(rdf_com.results.rdf))),fmt="%.4f")
# Plotting the RDF
#plt.figure(figsize=(8, 5))
#plt.plot(np.transpose(rdf_com.results.bins), np.transpose(rdf_com.results.rdf)], label="COM RDF (Group1 vs Group2)")
#plt.xlabel("Distance (Å)")
#plt.ylabel("g(r)")
#plt.title("Radial Distribution Function between Group1 and Group2")
#plt.legend()
#plt.grid()
"""
gA_coms=[group1.center_of_mass() for ts in u.trajectory]
gB_coms=[group2.center_of_mass(compound='residues') for ts in u.trajectory]
boxes=[u.dimensions for ts in u.trajectory]
if group1 == group2:
    com_arr  = [distances.self_distance_array(x,box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
    com_dist = [np.tile(x,2) for x in com_arr]
else:
    com_arr=[distances.distance_array(x,gB_coms[i],box=boxes[i],backend='OpenMP') for i,x in enumerate(gA_coms)]
    com_dist=[x.reshape(x.shape[0]*x.shape[1]) for x in com_arr]
print(com_dist[0])
bin_size = 100
max_val=10
# Create an array for the radial bins
bns = np.arange(0, max_val, bin_size)
hists=[np.histogram(x, bins=bns) for x in com_dist]
print(np.histogram(com_dist[0],bins=bns))
vols=[u.dimensions[0]*u.dimensions[1]*u.dimensions[2] for ts in u.trajectory]
avg_vol=np.mean(vols)
dens=len(com_dist[0])/avg_vol
sum_tmp=np.zeros(len(bns)-1)
for i in hists:
    sum_tmp=sum_tmp+i[0] 
sum_norm=(sum_tmp/(dens*4/3 * np.pi * (hists[0][1][1:]**3 - hists[0][1][:-1]**3)*(len(gA_coms)-1)))
print(sum_norm)
np.savetxt('rdf_mio.txt',list(zip(bns,sum_norm)),fmt='%.4f')
"""
