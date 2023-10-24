import MDAnalysis as mda
import itertools
import glob
import sys

def intervals_extract(iterable):
     
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def lastdcd_to_qmmm(topo,trj,sel):
    u = mda.Universe(topo,trj,in_memory=True)
    val=trj.split('-')[0].split('_')[-1]
    u.trajectory[-1]
    ag1 = u.select_atoms("all")
    ag1.write(f"QMMM_last_{val}.dcd")
    ag2 = u.select_atoms(sel)
    atom_idx=[x.index for x in ag2]
    atom_elm=[x.element.upper() for x in ag2 ]
    elms=[]
    for elm in atom_elm:
        if elm not in elms:
            elms.append(elm)
    QMmm_atoms={}
    for elm in elms:
        QMmm_atoms[elm]=[]

    for i,idx in enumerate(atom_idx):
        QMmm_atoms[atom_elm[i]].append(idx+1)
    intervals=list(intervals_extract(atom_idx))
    with open(f"vmd_index_{val}.dat",'w') as ofile:
        for interval in intervals:
            try:
                ofile.write(f"{interval[0]} to {interval[1]} ")
            except:
                ofile.write(f"{interval[0]} ")

    with open(f'QMMM_atoms_{val}.dat','w') as ofile:
        for i in QMmm_atoms:
            ofile.write("&QM_KIND {}\n ".format(i))
            ofile.write("MM_INDEX ")
            for j in QMmm_atoms[i]:
                ofile.write(f"{j} ")
            ofile.write("\n")
            ofile.write("&END QM_KIND\n")
    return()



topo=sys.argv[1]
inputs=glob.glob(sys.argv[2])
print(inputs)
sel='byres around 8 index 1636'

for input in inputs:
    print(input)
    lastdcd_to_qmmm(topo,input,sel)