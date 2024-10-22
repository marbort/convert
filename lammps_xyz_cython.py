
import lammpstrj_to_xyzbox_cython
import argparse


def main():
    # Parse the command line arguments
    parser = argparse.ArgumentParser(
        description='Convert DPData files to a single extended xyz formatted file')
    parser.add_argument('--input',  type=str, help='Input lammps trj')
    parser.add_argument('--output', type=str,
                        help='Output xyz file')
    parser.add_argument('--types_elm', dest='ty2elm', default=1,
                    type=str, help='types to element dictionary')
    args = parser.parse_args()
    
    
    lammpstrj_to_xyzbox_cython.lmp_to_xyzbox(args.input,args.output,args.ty2elm)


if __name__=="__main__":
    main()