import numpy as np
import itertools

def extract_QM_atoms(input):
    labels=np.loadtxt(input,skiprows=2,usecols=(1,2,3))
    index=[j for j,x in enumerate(labels) if x[0] > 0]
    return(index)

def intervals_extract(iterable):
     
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


input="Dimer_neutral_AcPh_QMMM_free_v3_adaptive-fmlabels-1_0.xyz"
index=extract_QM_atoms(input)
intervals=list(intervals_extract(index))
with open(f"vmd_index_adapt.dat",'w') as ofile:
        for interval in intervals:
            try:
                ofile.write(f"{interval[0]} to {interval[1]} ")
            except:
                ofile.write(f"{interval[0]} ")

