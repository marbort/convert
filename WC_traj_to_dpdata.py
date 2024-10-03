import numpy as np
import sys
import os



def read_traj_with_wannier(input):
    with open(sys.argv[1],'r') as ifile:
        lines=ifile.readlines()
    natoms=int(lines[0])
    frames=len(lines)//(natoms+2)
    atoms=[]
    elms=[]
    WC=[]
    WC_cut=[]
    for frame in range(frames):
        wc_tmp=[]
        atom_tmp=[]
        elm_tmp=[]
        for line in range(frame*(natoms+2)+2,frame*(natoms+2)+2+natoms):  
            if lines[line].split()[0] == "X":
                for x in range(1,4):
                    wc_tmp.append(lines[line].split()[x])
            else:
                elm_tmp.append(lines[line].split()[0])
                for x in range(1,4):
                    atom_tmp.append(lines[line].split()[x])
        WC_cut.append([wc_tmp[x] for x in range(len(atom_tmp))])
        WC.append(wc_tmp)
        atoms.append(atom_tmp)
        elms.append(elm_tmp)
    return(WC,WC_cut,atoms,elms)
                
def convert_WCtraj_to_dpdata(elms,coord,box,WC,type_map):
    
    os.mkdir('dpdata')
    os.mkdir('dpdata/set.000')
    with open('dpdata/type_map.raw','w') as ofile:
        for item in list(type_map.keys()):
            ofile.write(f"{item}\n")
    with open('dpdata/type.raw','w') as ofile:
        for elm in elms[0]:
            ofile.write(f"{type_map[elm]}\n")
    BOX=[]
    for i in range(len(coord)):
        BOX.append(box)
    BOX=np.array(BOX)
    
    np.save('dpdata/set.000/atomic_dipole.npy',WC)
    np.save('dpdata/set.000/coord.npy',coord)
    np.save('dpdata/set.000/box.npy',BOX)

def main():
    type_map={"C":0,"Cl":1,"H":2,"N":3,"O":4}
    dim=[sys.argv[2], sys.argv[3], sys.argv[4]]
    box=[dim[0], 0, 0, 0, dim[1], 0, 0, 0, dim[2]]
    WC,WC_cut,atoms,elms=read_traj_with_wannier(sys.argv[1])
    print(len(WC[0]),len(WC_cut[0]),len(atoms[0]))
    convert_WCtraj_to_dpdata(elms,atoms,box,WC_cut,type_map)

if __name__=="__main__":
    main()
    




            