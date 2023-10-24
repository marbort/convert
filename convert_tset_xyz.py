import glob
import os
import dpdata as dp

def convert_tset_xyz(root,iteration="Iteration_",tset="tset_"):
    files=glob.glob(os.path.join(root,iteration+"*",tset+"*"))
    xyzfolder=os.path.join(root,"tset_xyz")
    os.makedirs(xyzfolder,exist_ok=True)
    for file in files:
        data=dp.LabeledSystem(file,'deepmd/npy')
        
        xyzname=os.path.dirname(file).split('/')[-1]+"_"+os.path.basename(file)+".xyz"
        data.to('xyz',os.path.join(xyzfolder,xyzname))




convert_tset_xyz("/home/marco/SHARED/RATIO/WP1/test_check")
    