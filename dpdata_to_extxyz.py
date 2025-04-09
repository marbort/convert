#!/usr/bin/env python
import argparse
import sys
import dpdata
import ase
from ase import io
import os 

# Parse the command line arguments
parser = argparse.ArgumentParser(
    description='Convert DPData files to a single extended xyz formatted file')
parser.add_argument('--input', nargs='+', type=str, help='Input DPData type.raw files')
parser.add_argument('--output', type=str,
                    help='Output xyz file')
parser.add_argument('--offset', type=int, default=1,help='Save every nth frame')
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
        Input DPData files type.raw
    fname : str
        Output xyz file
    return : None
    """

    if isinstance(dpdata_files,list):
        for dpdata_file in dpdata_files:
            # Read the DPData type.raw file and extract parent folder
            path=os.path.dirname(dpdata_file)
            try:
                system = dpdata.LabeledSystem(path, fmt='deepmd/npy')
            except:
                system = dpdata.System(path, fmt='deepmd/npy')
                print(f"Energies not found. Using non labeled system")
            # for each frame
            for ii in range(0,system.get_nframes(),args.offset):
                # read the frame
                positions = system['coords'][ii]
                symbols = [system['atom_names'][type]
                        for type in system['atom_types']]
                cell = system['cells'][ii]
                pbc = True
                try:
                    forces = system['forces'][ii]
                    energies = system['energies'][ii]
                except:
                    pass

                # make an ase.Atoms object
                atoms = ase.Atoms(
                    symbols=symbols, positions=positions, cell=cell, pbc=pbc)
                try:
                    atoms.arrays['forces'] = forces
                    atoms.info['energy'] = energies
                except:
                    pass
                if prefix:
                    # write the frame to the extended xyz file
                    path_str=path.replace("/","_")
                    fname=f"{prefix}_{path_str}.xyz"
                    ase.io.write(os.path.join(path,fname), atoms, format='extxyz', append=True,write_results=False)
                else:
                    name=f"trj_{dpdata_file}.xyz"
                    ase.io.write(name, atoms, format='extxyz', append=True,write_results=False)
            print(fname)


# Define the main function
if __name__ == '__main__':
    make_xyz_file(args.input, args.output)
