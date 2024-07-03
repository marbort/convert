import glob
import os
import dpdata as dp
import subprocess
from Create_plumed_groups_xyz import create_plumed_input

def plumed_driver(input):
    create_plumed_input(input)
    subprocess.run(["plumed", "driver", "--plumed", "plumed.dat", "--ixyz", input, "--timestep", "1", "--trajectory-stride", "1", "--length-units", "A"])



files=glob.glob("Iteration*")
for file in files:
    plumed_driver(file)





