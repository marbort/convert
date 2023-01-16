import argparse
import numpy as np
import os

parser = argparse.ArgumentParser(
    description='Fix type_map.raw and type.raw, to use another type_map.raw which has a consistent type map.')
parser.add_argument(
    '--folder', help='Folder containing original type.raw and type_map.raw')
parser.add_argument('--new_type_map', help='Type map to use.')
args = parser.parse_args()

# Read type_map.raw
type_map_dict = {}
with open(os.path.join(args.folder, 'type_map.raw'), 'r') as f:
    for i, line in enumerate(f):
        type_map_dict[str(i)] = line.strip()

# Read type.raw
types = np.loadtxt(os.path.join(args.folder, 'type.raw'), dtype=int)


# Read new type_map.raw
type_map_dict_new = {}
with open(args.new_type_map, 'r') as f:
    for i, line in enumerate(f):
        type_map_dict_new[str(i)] = line.strip()

inv_type_map_dict_new = {v: k for k, v in type_map_dict_new.items()}
# Create new type.raw
types_new = np.zeros_like(types)
for i, t in enumerate(types):
    types_new[i] = int(inv_type_map_dict_new[type_map_dict[str(t)]])

# Write new type.raw and type_map.raw
np.savetxt(os.path.join(args.folder, 'type.raw'), types_new, fmt='%d')
np.savetxt(os.path.join(args.folder, 'type_map.raw'),
           np.array(list(type_map_dict_new.values())), fmt='%s')
