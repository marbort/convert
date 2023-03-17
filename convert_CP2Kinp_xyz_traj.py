
import os

root=os.getcwd()

crds=[]
for file in os.listdir(root):
    if file.endswith('.inp'):
        crds_file=[]
        with open(file,'r') as ifile:
            lines=ifile.readlines()
            start=False
            for line in lines:
                    if "&END COORD" in line:
                        start=False
                    if start:
                        crds_file.append(line)
                    if "&COORD" in line:
                        start=True
            crds.append(crds_file)
print(len(crds[0]))
with open(os.path.join(root,'input_traj.xyz'),'w') as ofile:
    for i,x in enumerate(crds):
        ofile.write("{}\n Input {}\n".format(len(x),i))
        for line in x:
            ofile.write("{}".format(line))

                    
                        
                        