import dpdata as dp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, default='./conf.data', help='input lammps data file')
parser.add_argument('--types', type=str, nargs='+', default=['O', 'H', 'C'], help='Elements of atom types')
args = parser.parse_args()

data=dp.System(args.input, fmt='lammps/lmp')
for i in range(len(data['atom_names'])):
    data['atom_names'][i]=args.types[i]
    

data.to('xyz', f'{args.input}.xyz')


