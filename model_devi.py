def CalcModelDeviation(args):
    from deepmd.infer import calc_model_devi
    from deepmd.infer import DeepPot as DP
    import dpdata
    import numpy as np
    # Load the graphs
    graphs = [DP(g) for g in args.graphs]

    # Load the configurations, if first traj file crashes,

    # Load the configurations
    print('Loading the configurations...')
    start = 0
    while True:
        assert start < len(args.traj), 'No valid traj file found!'

        try:
            data = dpdata.System(args.traj[start], fmt='lammps/dump')
            data = data.sub_system(range(1, len(data['coords'])))
            start+=1
            break
        except:
            start += 1
            pass

    for i in range(start, len(args.traj)):
        try:
            dat = dpdata.System(args.traj[i], fmt='lammps/dump')

            # remove first frame
            dat = dat.sub_system(range(1, len(dat['coords'])))
            data.append(dat)
        except:
            pass

    # Perturb the configurations before computing the model deviations
    if args.perturb:
        print('Perturbing the configurations before computing the model deviations...')
        data = data.perturb(pert_num=2,
                            cell_pert_fraction=0.05,
                            atom_pert_distance=0.01,
                            atom_pert_style='normal')

    # Compute the model deviations
    print('Computing the model deviations...')
    natoms = int(len(data['coords'][0]))
    active_learning_configs = []
    with open('bad_configurations.dat', 'w') as fp:
        for i in np.array(range(len(data['coords'])))[::args.offset]:
            coord = data['coords'][i].reshape([1, natoms, 3])

            cell = data['cells'][i].reshape([1, -1])
            atype = data['atom_types']

            model_devi = calc_model_devi(coord, cell, atype, graphs)

            model_devi = model_devi[0][-3]
            if model_devi > args.f_min and model_devi < args.f_max:
                active_learning_configs.append(i)
            fp.write('{} {}\n'.format(i, model_devi))

    if len(active_learning_configs) == 0:
        print('No configurations found within the specified deviation range.')
        return

    # Pick the configurations within specified deviation range
    data = data.sub_system(active_learning_configs)

    # Pick the first args.max_num configurations
    data = data.sub_system(range(min(args.max_num, len(data['coords']))))

    # Save the picked configurations
    print('Saving the picked configurations...')
    data.to_deepmd_npy(args.out)


def main():
    import argparse
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='This script picks configurations\
                                    from a lammps trajectory file based on the\
                                        model deviations between the different models.')

    parser.add_argument('--offset', dest='offset', default=1,
                        type=int, help='Use only every nth frame.')
    parser.add_argument('--graphs', type=str, nargs='+',
                        help='List of deepmd graphs to compute force model deviations.')
    parser.add_argument('--traj', default='lmp.lammpstrj', type=str, nargs='+',
                        help='List of deepmd graphs to compute force model deviations.')
    parser.add_argument('--f_max', default='0.5', type=float,
                        help='Maximum allowed disagreement in forces between models.')
    parser.add_argument('--f_min', default='0.02', type=float,
                        help='Minimum allowed disagreement in forces between models.')
    parser.add_argument('--out', default='dpdata', type=str,
                        help='Folder name used for storing picked configurations.')
    parser.add_argument('--perturb', default=False, action='store_true',
                        help='Perturb the configurations before computing the model deviations.')
    parser.add_argument('--max_num', default=1000, type=int,
                        help='Maximum number of configurations to be picked.')
    args = parser.parse_args()
    print(args)
    CalcModelDeviation(args)


if __name__ == '__main__':
    main()
