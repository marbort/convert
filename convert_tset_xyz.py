import glob
import os
import dpdata as dp

def convert_tset_xyz(file,path):
    xyzfolder=os.path.join(path,f"{file}_xyz")
    os.makedirs(xyzfolder,exist_ok=True)
    box=[]
    data=dp.LabeledSystem(file,'deepmd/npy')
    box.append(data['cells'])
    xyzname=file+".xyz"
    name=os.path.join(xyzfolder,xyzname)
    print(name)
    data.to('xyz',name)
    return(box,name)

def add_box_toxyz(file,box):
    with open(file,'r') as ifile:
        lines=ifile.readlines()
    natoms=int(lines[0].rstrip())
    
    for i,line in enumerate(lines[1::natoms+2]):
        cell=f"{box[0][i][0][0]} {box[0][i][1][1]} {box[0][i][2][2]}\n"
        lines[i*(natoms+2)+1]=line.replace("\n",cell)
    with open(file,'w') as ofile:
        for line in lines:
            ofile.write(line)

    
    
    





path=("./")
files=glob.glob(path+"tset*")
for file in files:
    box,name=convert_tset_xyz(file,path)
    add_box_toxyz(name,box)
    