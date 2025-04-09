import numpy as np
import MDAnalysis as mda

def subst_resid(topo,coords):
    IM2THFS = {"THF1Mg1":(127,139),"THF2Mg1":(140,152),"THF1Mg2":(1647,1659),"THF2Mg2":(1660,1672)}
    u=mda.Universe(topo,coords)
    
    sel="not resid 4 44 and name O1 and resname THF and around 10 bynum 117 1637"
    
    THF=u.select_atoms(sel)
    Mg1=u.select_atoms("bynum 117")
    Mg2=u.select_atoms("bynum 1637")
    
    
    IM2THFS_orig_positions = [ u.atoms[IM2THFS[item][0]:IM2THFS[item][1]+1].positions for item in IM2THFS]
  
    
    dist_Mg1=[np.linalg.norm(x-Mg1.positions[0]) for x in THF.positions]
    dist_Mg2=[np.linalg.norm(x-Mg2.positions[0]) for x in THF.positions]
    
    close_THFS=[
        
        THF.residues[np.argsort(dist_Mg1)[0]].ix,THF.residues[np.argsort(dist_Mg1)[1]].ix, 
        THF.residues[np.argsort(dist_Mg2)[0]].ix,THF.residues[np.argsort(dist_Mg2)[1]].ix 
                
                
                ]
    close_THFS_indexs=[ u.residues[x].atoms.ix for x in close_THFS]
    
    print(close_THFS_indexs)
    print(close_THFS)

    for i,item in enumerate(IM2THFS):
        #print(f"Before {IM2THFS_orig_positions[i].positions}")
        u.atoms[IM2THFS[item][0]:IM2THFS[item][1]+1].positions = u.residues[close_THFS[i]].atoms.positions
        u.residues[close_THFS[i]].atoms.positions = IM2THFS_orig_positions[i]
        #print(f"After {IM2THFS_orig_positions[i].positions}")
    
    u.atoms.write("test.pdb")
    return(close_THFS_indexs,IM2THFS)
"""
    print(f"Closest to Mg1: {THF.residues[np.argmin(dist_Mg1)].ix}")
    print(f"Closest to Mg1: {THF.residues[np.argsort(dist_Mg1)[0]]}")
    print(f"2nd Closest to Mg1: {THF.residues[np.argsort(dist_Mg1)[1]]}")
    print(f"Closest to Mg2: {THF.residues[np.argsort(dist_Mg2)[0]]}")
    print(f"2nd Closest to Mg2: {THF.residues[np.argsort(dist_Mg2)[1]]}")
    
    #for item in IM2THFS:
        
    
    print(u.atoms[1647:1657].positions)
    u.atoms[1647:1659].positions=u.atoms[140:152].positions
    print(u.atoms[1647:1659].positions)
"""
subst_resid('../QMMM_Reactivity_Charged_AcPh_OKLJ.parm7','Dimer_charged_AcPh_QMMM_Belim_v6_1.1-pos-1.dcd.pdb')
        
    