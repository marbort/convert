import MDAnalysis as mda
#from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import contacts

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



u = mda.Universe('DES_BIG_MOD.parm7', 'npt_iso_C_rescale_fromannealing_090_100ns_okparams_whole.xtc')

sel_ref = "byres resid 1380"
sel_group = "byres resname CHL"
ref = u.select_atoms(sel_ref).center_of_mass(compound='residues')
group = u.select_atoms(sel_group).center_of_mass(compound='residues')


def contacts_within_cutoff(u, group_a, group_b, radius=7.7):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a, group_b)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

ca = contacts_within_cutoff(u, ref, group, radius=7.7)
print(ref,group[2])
print(ca[23])