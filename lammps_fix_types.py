import dpdata as dp
import numpy as np
import argparse
import re
import os

# Parse the command line arguments
parser = argparse.ArgumentParser(description='Fix a  lammps data file\
     with new type map')

parser.add_argument('--i', default='conf.data', type=str, help='File name \
    storing deepmd/npy files.')
parser.add_argument('--o', default='conf.data', type=str, help='Name of \
    lammps data file.')
parser.add_argument('--type_map',  type=str, help='New\
     type map.')

args = parser.parse_args()
def parse_lmp(input):
    # Load the dpdata folder
    data=dp.System(input,'lammps/lmp')
    with open(input,'r') as ifile:
        lines=ifile.read()
    masses=[(x.split()[0],x.split()[1]) for x in re.findall(r'Masses.*Atoms',lines,re.DOTALL)[0].split('\n')[2:-2]]
    print("Found Masses: ",masses)
    return(data,masses)

data,masses=parse_lmp(args.i)


#load the type map
with open(args.type_map,'r') as ifile:
    lines=ifile.readlines()
type_map=[(x.rstrip(),i+1,dp.periodic_table.Element(x.rstrip()).mass) for i,x in enumerate(lines)]
print("New Type Map: ",type_map)

old_to_new=[]
for i,mass_conf in enumerate(masses):
    for j,mass_type in enumerate(type_map):
        if float(mass_conf[1]) == mass_type[2]:
            old_to_new.append((i,j))
print("Conversions: ", old_to_new)

for pair in old_to_new:
    if pair[0] == pair[1]:
        pass
    else:
        for k,type in enumerate(data['atom_types']):
            if type == pair[0]:
                data['atom_types'][k]=pair[1]
            else:
                pass
new_at_count=[]
for i in range(len(masses),len(type_map)):
        data['atom_names'].append("Type_{}".format(type_map[i][1]))
for i in range(len(type_map)):
        new_at_count.append(np.count_nonzero(data['atom_types'] == i))
for i,count in enumerate(new_at_count):
    try:
        data['atom_numbs'][i]=count
    except:
        data['atom_numbs'].append(count)

    

print(data['atom_numbs'])
#1os.rename("conf.data", "conf.data.old_tmap")

data.to_lammps_lmp(args.o)

masses_new = [dp.periodic_table.Element(
        name[0]).mass for name in type_map]
text = " Masses\n\n"
for i, m in enumerate(masses_new):
    text += "{} {}\n".format(i+1, m)
text += '\n'
text += ' Atoms'
with open(args.o, 'r') as fp:
    content = fp.read()
    content = re.sub("Atoms", text, content)
with open(args.o, 'w') as fp:
    fp.write(content)

        
    
    
    
    
    


