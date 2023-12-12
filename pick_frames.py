
import MDAnalysis as mda
import os

#from MDAnalysis.tests.datafiles import PDB, GRO, XTC

#os.chdir('/Users/Kiki/clusters/saga/Validation/LiX-opes2d/LiI/runB')

input_file='lmp.lammpstrj'
output_file='lmp_selected.lammpstrj'

frame_interval = 100
# Load the input trajectory
u = mda.Universe('lmp.lammpstrj', format='LAMMPSDUMP')
all = u.select_atoms("all")


for ts in u.trajectory:
    if ts.frame % frame_interval == 0:
        # Write the current frame to the output file
        all.write(output_file, format='LAMMPSDUMP')
