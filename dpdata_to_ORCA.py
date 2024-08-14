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
parser.add_argument('--template', dest='template', default='ORCA_template.inp',
                    type=str, help='File path to template file.')

args = parser.parse_args()

data = dpdata.System(args.i, 'deepmd/npy')
box_0=data['cells'][0][0]+data['cells'][0][1]+data['cells'][0][2]

#wrap
for frame in range(len(data['coords'])):
    transl=data['coords'][frame][0]-(box_0)/2
    data['coords'][frame]-=transl
    for atom in range(len(data['coords'][frame])):
        for dim in range(len(data['coords'][frame][atom])):
            if data['coords'][frame][atom][dim] > box_0[dim]:
                data['coords'][frame][atom][dim]-=box_0[dim]
            if data['coords'][frame][atom][dim] < 0:
                data['coords'][frame][atom][dim]=box_0[dim]+data['coords'][frame][atom][dim]
                
print(data['coords'][0][0])
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

    print("Converting {} frames, out of {} frames available, to ORCA input.".format(len(frames),len(data['coords'])))

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
        


       # if os.path.isdir("inputs"):
            #print("Creating input for frame {}".format(i))
        if not os.path.isdir("inputs"):
            os.mkdir("inputs")
            #print("Creating input for frame {}".format(i))
        with open('inputs/orca_{}_{:06d}.inp'.format(args.name, i), 'w') as ofile:
            ofile.write(struct)
