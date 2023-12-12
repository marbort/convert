#!/usr/bin/env python
import argparse
import sys
import dpdata
import ase
import os 

# Parse the command line arguments
parser = argparse.ArgumentParser(
    description='Convert DPData files to a single extended xyz formatted file')
parser.add_argument('--input', nargs='+', type=str, help='Input DPData files')
parser.add_argument('--output', type=str,
                    help='Output xyz file')
args = parser.parse_args()

if args.input is None:
    print('No input files specified')
    sys.exit(1)

def make_xyz_file(dpdata_files, prefix):
    """
    Convert DPData files to a single extended xyz formatted file

    Parameters
    ----------
    dpdata_files : list of str
        Input DPData files
    fname : str
        Output xyz file
    return : None
    """

    if isinstance(dpdata_files,list):
        for dpdata_file in dpdata_files:
            # Read the DPData type.raw file and extract parent folder
            path=os.path.dirname(dpdata_file)
            system = dpdata.LabeledSystem(path, fmt='deepmd/npy')
            # for each frame
            for ii in range(system.get_nframes()):
                # read the frame
                positions = system['coords'][ii]
                symbols = [system['atom_names'][type]
                        for type in system['atom_types']]
                cell = system['cells'][ii]
                pbc = True
                forces = system['forces'][ii]
                energies = system['energies'][ii]

                # make an ase.Atoms object
                atoms = ase.Atoms(
                    symbols=symbols, positions=positions, cell=cell, pbc=pbc)
                atoms.arrays['forces'] = forces
                atoms.info['energy'] = energies
                if prefix:
                    # write the frame to the extended xyz file
                    path_str=path.replace("/","_")
                    fname=f"{prefix}_{path_str}.xyz"
                    ase.io.write(os.path.join(path,fname), atoms, format='extxyz', append=True)
                else:
                    name=f"trj_{dpdata_file}.xyz"
                    ase.io.write(name, atoms, format='extxyz', append=True)
            print(fname)


# Define the main function
if __name__ == '__main__':
    make_xyz_file(args.input, args.output)
