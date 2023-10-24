import dpdata
import numpy as np
import argparse
import re

# Parse the command line arguments
parser = argparse.ArgumentParser(description='Convert a  dpdata folder\
     lammps data file')

parser.add_argument('--i', default='dpdata', type=str, help='Folder name \
    storing deepmd/npy files.')
parser.add_argument('--o', default='conf.data', type=str, help='Name of \
    lammps data file.')
parser.add_argument('--num', default=1, type=int, help='Number of\
     inputs to create.')

args = parser.parse_args()

# Load the dpdata folder
data = dpdata.System(args.i, fmt='deepmd/npy')

# Convert to lammps data file
nframes = len(data['cells'])
frames = [0]
if args.num > 1:
    frames = np.random.choice(np.array(range(nframes)),
                              size=args.num, replace=False)
for i in range(args.num):
    if args.num == 1:
        name = args.o
    else:
        name = args.o.replace('.data', '{}.data'.format(i))
    data.to_lammps_lmp(name, frame_idx=frames[i])

    masses = [dpdata.periodic_table.Element(
        name).mass for name in data['atom_names']]
    text = " Masses\n\n"
    for i, m in enumerate(masses):
        text += "{} {}\n".format(i+1, m)
    text += '\n'
    text += ' Atoms'
    with open(name, 'r') as fp:
        content = fp.read()
        content = re.sub("Atoms", text, content)
    with open(name, 'w') as fp:
        fp.write(content)
