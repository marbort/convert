import argparse
import MDAnalysis as md

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", help="LAMMPS trajectory file", default="trajectory.lammpstrj")
    parser.add_argument("--output_pdb", help="Output PDB file for the first frame", default="first_frame.pdb")
    parser.add_argument("--output_dcd", help="Output DCD file for all frames", default="trajectory.dcd")
    parser.add_argument("--type_map", help="Type map file", default="type_map.raw")
    args = parser.parse_args()
    
    # Load the type map file
    type_map = {}
    with open(args.type_map) as f:
        for i, line in enumerate(f):
            type_map[i + 1] = line.strip()
    
    # Load the LAMMPS trajectory
    u = md.Universe(args.input_file, format='LAMMPSDUMP')
    
    # Assign atom names based on type map
#    for atom in u.atoms:
 #       atom.name = type_map[int(atom.type)-1]

    print(type_map)
    u.add_TopologyAttr('name',[type_map[int(atom.type)] for atom in u.atoms])
    
    # Write the first frame to a PDB file
    u.atoms.write(args.output_pdb)
    
    # Write all frames to a DCD file
    with md.Writer(args.output_dcd, n_atoms=u.atoms.n_atoms) as dcd:
        for ts in u.trajectory:
            dcd.write(u.atoms)

if __name__ == "__main__":
    main()
