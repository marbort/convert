import itertools
import glob
import sys

def intervals_extract(iterable):
     
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


def extract_QMMM(file):
    with open(file,'r') as ifile:
        name=ifile.name.replace(".inp","")
        lines=ifile.readlines()
    start=False
    atoms=""
    for line in lines:
        if "&END QM_KIND" in line:
            start=False
        if "MM_INDEX" in line:
            start=True
        if start:
            if "MM_INDEX" in line:
                newline=line.replace("MM_INDEX","")
                atoms+=newline
            else:
                atoms+=line
    atoms_list=atoms.split()
    print(len(atoms_list))
    atom_idx=sorted([int(x)-1 for x in atoms_list])
    intervals=list(intervals_extract(atom_idx))
    with open(f"vmd_index_{name}.dat",'w') as ofile:
        for interval in intervals:
            try:
                ofile.write(f"{interval[0]} to {interval[1]} ")
            except:
                ofile.write(f"{interval[0]} ")
                

extract_QMMM(sys.argv[1])

