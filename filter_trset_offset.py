import dpdata as dp
import glob
import sys
import os
from pathlib import Path



files=glob.glob(f"{sys.argv[1]}/**/tset*",recursive=True)
dest=sys.argv[2]

for file in files:
    tset=os.path.basename(file)
    Iteration=file.split('/')[-2]
    Residue=file.split('/')[-4]

    name=f"{tset}_{Residue}_{Iteration}"
    print(name)

    data=dp.LabeledSystem(file,'deepmd/npy')
    data2=data.sub_system(range(0,len(data),10))
    data2.to('deepmd/npy',f"{dest}/{name}")


#print(Iteration,Residue)



#for file in files
