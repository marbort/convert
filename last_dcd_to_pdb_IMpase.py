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

def lastdcd_to_qmmm(topo,trj,sel,check,test=None):
    u = mda.Universe(topo,trj,in_memory=True)
    check_atms=u.select_atoms(check)
    if len(check_atms) == 0:
        print("#"*80)
        print(f"ERROR in {trj} due to no {check} atoms")
        print("#"*80)
        sys.exit()
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
    chrg=0
    with open(f"vmd_index_{val}_{chrg}{test}.dat",'w') as ofile:
        for interval in intervals:
            try:
                ofile.write(f"{interval[0]} to {interval[1]} ")
            except:
                ofile.write(f"{interval[0]} ")
    with open(f"residues_{val}_{chrg}{test}.dat",'w') as ofile:
        for resid in QMmm_residues:
                ofile.write(f"{resid} ")


    with open(f'QMMM_atoms_{val}_{chrg}{test}.dat','w') as ofile:
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
radius=sys.argv[3]
print(inputs)


#sel='resid 4 44 or (not resname IM2 and byres around 8 index 1636)'


for input in inputs:
    #Cl1_pos,Cl2_pos,Mg1_pos,Mg2_pos,avg,avg_Mg=get_avg_Cl(4,46,topo,input)
    #Cl_avg_pos=" ".join([str(x) for x in avg])
    #Mg_avg_pos=" ".join([str(x) for x in avg_Mg])
    #sel=f"byres point {Cl_avg_pos} 8.0"
    #sel=f"resid 4 44 or (not resname IM2 and (byres point {Cl_avg_pos} {radius} and not resid 6223))"
    check="resname L1A"
    sel_test="resname L1A"
    sel=f"not resname WAT and (resname L1A or byres around {radius} resname L1A)"
    #sel=f"bynum 115 to 126 or bynum 1635 to 1646 or (not resname IM2 and (byres point {Mg_avg_pos} {radius}))"
    print(input)
    #print(Cl1_pos,Cl2_pos,Cl_avg_pos)
    print(f"Extracting atoms from sel: {sel}")
    lastdcd_to_qmmm(topo,input,sel,check)
    lastdcd_to_qmmm(topo,input,sel_test,check,"test")