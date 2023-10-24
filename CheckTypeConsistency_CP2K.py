import MDAnalysis
import dpdata
 
def crds_cp2k(file):
    start = False
    startb = False
    crds = []
    box = []
    types_map = []
    hartree_to_ev=1
    with open(file, 'r') as ifile:
            lines = ifile.readlines()
            
    for line in lines:
        # Extract coordinates from CP2K files
        if "&END COORD" in line:
            start = False
            # print(line)
        if start:
            # print(line)
            tmp = [float(line.split()[x]) for x in range(1, 4)]
            crds = crds+tmp
            # print(crds,len(crds))
        if "&COORD" in line:
            start = True

    return (crds)


def CheckConsistency(xyz_fn, dpdata_folder):
    #u = MDAnalysis.Universe(xyz_fn)
    start = False
    startb = False
    crds = []
    box = []
    types_map = []
    hartree_to_ev=1
    with open(file, 'r') as ifile:
            lines = ifile.readlines()
            
    for line in lines:
        # Extract coordinates from CP2K files
        if "&END COORD" in line:
            start = False
            # print(line)
        if start:
            # print(line)
            tmp = [float(line.split()[x]) for x in range(1, 4)]
            crds = crds+tmp
            # print(crds,len(crds))
        if "&COORD" in line:
            start = True

    atom_names = [x.split()[0] for x in crds]

    dpdata_system = dpdata.System(dpdata_folder, fmt='deepmd/npy')
    dpdata_atom_names = [dpdata_system.get_atom_names()[t]
                         for t in dpdata_system.get_atom_types()]

    # Check that each atom name of atom_names and dpdata_atom_names are the same
    if len(atom_names) != len(dpdata_atom_names):
        print('Number of atoms are different: {} vs {}'.format())
        return False
    for i in range(len(atom_names)):
        if atom_names[i] != dpdata_atom_names[i]:
            print('Atom names are different: {} vs {}'.format(
                atom_names[i], dpdata_atom_names[i]))
            return False

    return True


# Main function parsing command line arguments
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Check the consistency of xyz file and dpdata folder')
    parser.add_argument('--xyz_fn', help='xyz file')
    parser.add_argument('--dpdata_folder', help='dpdata folder')
    args = parser.parse_args()

    print('Checking {} and {}'.format(args.xyz_fn, args.dpdata_folder))

    is_consistent = CheckConsistency(args.xyz_fn, args.dpdata_folder)
    if is_consistent:
        print('Consistent')
    else:
        print('Inconsistent')
