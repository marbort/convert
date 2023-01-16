import json
import os
from re import T
import numpy as np
import argparse
import dpdata
boxl = ["A", "B", "C"]
name = "THF_iPrMgCl_CP2K"


parser = argparse.ArgumentParser(description='Plot data')
parser.add_argument('--i', dest='i', default='./dpdata',
                    help='input folder in dpdata format.')

parser.add_argument('--max_frames', type=int, dest='max_frames', default=1000000,
                    help='Maximum number of frames to convert. Is an alternative to --offset.')

parser.add_argument('--offset', dest='offset', default=1,
                    type=int, help='save only every nth frame')
parser.add_argument('--name', dest='name', default=1,
                    type=str, help='Project name')
parser.add_argument('--template', dest='template', default='cp2k_template.inp',
                    type=str, help='File path to template file.')

args = parser.parse_args()

data = dpdata.System(args.i, 'deepmd/npy')


with open(args.template, 'r') as ifile:
    lines = ifile.read()

    if args.max_frames < len(data['coords']):
        # choose args.max_frames from an array of length len(data['coords']), ignore offset
        # Make sure that first frame is ignored and that frames has exactly args.max_frames elements
        # Make sure that the frames are evenly spaced from the beginning to end
        frames = np.array(range(0, len(data['coords']), len(data['coords'])//args.max_frames))[1:]

        # remove superfluous frames randomly, use always the same seed to make sure that the same frames are removed
        np.random.seed(42)
        if len(frames) > args.max_frames:
            frames = np.random.choice(frames, args.max_frames, replace=False)

    else:
        frames = np.array(range(0, len(data['coords']), args.offset))

    print("Converting {} frames, out of {} frames available, to cp2k input.".format(len(frames),len(data['coords'])))

    for i in frames:
        crds=[]
        for j, item in enumerate(data['coords'][i]):
            X=item[0]
            Y=item[1]
            Z=item[2]
            crds.append("       "+data['atom_names'][data['atom_types'][j]] +
                        "   {:13.9f}   {:13.9f}   {:13.9f}".format(X, Y, Z))
            join='\n'.join(crds)
        struct=lines.replace("##coord##", join)
        box=[]
        for k, entry in enumerate(data['cells'][i]):
            X=entry[0]
            Y=entry[1]
            Z=entry[2]
            box.append(
                "       " + boxl[k]+"   {:13.9f}   {:13.9f}   {:13.9f}".format(X, Y, Z))
            bjoin='\n'.join(box)
        struct=struct.replace("##cell##", bjoin)
        struct=struct.replace(
            "##project##", "{}_{:06d}".format(args.name, i))

       # if os.path.isdir("inputs"):
            #print("Creating input for frame {}".format(i))
        if not os.path.isdir("inputs"):
            os.mkdir("inputs")
            #print("Creating input for frame {}".format(i))
        with open('inputs/cp2k_{}_{:06d}.inp'.format(args.name, i), 'w') as ofile:
            ofile.write(struct)
