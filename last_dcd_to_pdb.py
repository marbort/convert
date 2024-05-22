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
    #ag1.write(f"QMMM_last_{val}.dcd")
    #ag1.write(f"QMMM_last_{val}.nc")
    ag2 = u.select_atoms(sel)
    atom_idx=[x.index for x in ag2]
    QMmm_residues=sorted(list(set([x.resid for x in ag2])))
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
    if 'N' not in list(QMmm_atoms.keys()):
        chrg=len(QMmm_atoms['MG'])-len(QMmm_atoms['CL'])
    else:
        chrg=len(QMmm_atoms['MG'])+len(QMmm_atoms['N'])-len(QMmm_atoms['CL'])
    with open(f"vmd_index_{val}_{chrg}.dat",'w') as ofile:
        for interval in intervals:
            try:
                ofile.write(f"{interval[0]} to {interval[1]} ")
            except:
                ofile.write(f"{interval[0]} ")
    with open(f"residues_{val}_{chrg}.dat",'w') as ofile:
        for resid in QMmm_residues:
                ofile.write(f"{resid} ")


    with open(f'QMMM_atoms_{val}_{chrg}.dat','w') as ofile:
        for i in QMmm_atoms:
            ofile.write("&QM_KIND {}\n ".format(i))
            ofile.write("MM_INDEX ")
            for j in QMmm_atoms[i]:
                ofile.write(f"{j} ")
            ofile.write("\n")
            ofile.write("&END QM_KIND\n")
    return()

def get_avg_Cl(Cl1,Cl2,topo,trj):
    u = mda.Universe(topo,trj,in_memory=True)
    val=trj.split('-')[0].split('_')[-1]
    u.trajectory[-1]
    Cl1_sel=u.select_atoms(f"name Cl1 and resid {Cl1}")
    Cl2_sel=u.select_atoms(f"name Cl1 and resid {Cl2}")
    Mg1_sel=u.select_atoms(f"name Mg1 and resid {Cl1}")
    Mg2_sel=u.select_atoms(f"name Mg1 and resid {Cl2}")
    Cl1_pos=Cl1_sel.positions[0]
    Cl2_pos=Cl2_sel.positions[0]
    Mg1_pos=Mg1_sel.positions[0]
    Mg2_pos=Mg2_sel.positions[0]
    avg=(Cl1_pos+Cl2_pos)/2
    avg_Mg=(Mg1_pos+Mg2_pos)/2
    return(Cl1_pos,Cl2_pos,Mg1_pos,Mg2_pos,avg,avg_Mg)
    
    
    


topo=sys.argv[1]
inputs=glob.glob(sys.argv[2])
radius=sys.argv[3]
print(inputs)


#sel='resid 4 44 or (not resname IM2 and byres around 8 index 1636)'


for input in inputs:
    Cl1_pos,Cl2_pos,Mg1_pos,Mg2_pos,avg,avg_Mg=get_avg_Cl(4,44,topo,input)
    Cl_avg_pos=" ".join([str(x) for x in avg])
    Mg_avg_pos=" ".join([str(x) for x in avg_Mg])
    #sel=f"byres point {Cl_avg_pos} 8.0"
    #sel=f"resid 4 44 or (not resname IM2 and (byres point {Cl_avg_pos} {radius} and not resid 6223))"
    sel=f"resid 4 44 or (not resname IM2 and (byres point {Mg_avg_pos} {radius} and not resid 6223))"
    print(input)
    print(Cl1_pos,Cl2_pos,Cl_avg_pos)
    lastdcd_to_qmmm(topo,input,sel)