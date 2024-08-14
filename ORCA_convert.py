import json
import os
from re import T
import numpy as np
import argparse
import shutil

parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--offset', dest='offset', default=1,
                    type=int, help='save only every nth frame')
parser.add_argument('--feval', dest='feval', action='store_true',
                    help='Extract forces in Force Eval format')
parser.add_argument('--types_map', dest='types_map', 
                    help='Use supplied type map. If none use atom ordering in CP2K file')




parser.add_argument('--box', dest='box', type=float,
                    help='XYZ Box Sizes (orthrombic only)', nargs=3)
args = parser.parse_args()

hartree_to_ev = 27.211399
# print(re.search(r'(?<=Element).*?(?=SUM OF)',lines,re.DOTALL)[0])

# This functions extracts coordinates and atom types from a CP2K trajectory
# positions xyz file


def extract_ORCA_coords(file, kind, tmap="none"):
    traj_ener = np.empty(0)
    types_map = {}
    types = []
    numbers = []
    natom = 0
    box_0 = [[20.080, 0., 0., 0., 20.080, 0., 0., 0., 20.080]]
    # with open("/run/media/marco/SAMSUNG/SHARED/SHARED_bkup/dottorato/scripts/PeriodicTableJSON.json",'r') as jfile:
    # jsonf=json.load(jfile)
    # atoms_dict={x['symbol']:x['number'] for x in jsonf['elements']}

    box = np.array(box_0)  # fixed box params taken from input
    box_traj = np.array(box_0)
    # read pos file for coordinates and atom types
    lines = file.readlines()

    natom = int(lines[0])
    atoms = [lines[x+2].split()[0] for x in range(natom)]
    nframes = len(lines)//(natom+2)
    print("Found {:d} frames".format(nframes))
    print("Extracting {:d} frames".format(nframes//args.offset))
    traj_ener = [np.append(traj_ener, float(lines[args.offset*k*(natom+2)+1].split()[8].rstrip())*hartree_to_ev)
                 for k in range(nframes//args.offset)]

    time = [float(lines[args.offset*k*(natom+2)+1].split()[5].strip(','))
            for k in range(max(nframes//args.offset,1))]
    print(time[0])

    # group val in an array for each frame
    val = [[float(lines[args.offset*k*(natom+2)+i+2].split()[l])
           for i in range(natom) for l in range(1, 4)]
           for k in range(nframes//args.offset)]

    # group val in an array for each element
    val_sch = [[[float(lines[args.offset*k*(natom+2)+i+2].split()[l])
                 for l in range(1, 4)] for i in range(natom)]
               for k in range(nframes//args.offset)]

    index = 0
    try: 
        with open(tmap, 'r') as ifile:
                        lines=ifile.readlines()
                        for line in lines:
                            types_map[line.split()[0]] = index
                            index = index+1
                        for i in atoms:
                            types.append(types_map[i])
    except:
        print("Type Map file not found or formatted properly")
        raise SystemExit(1)
    
    return (numbers, time, val, val_sch, traj_ener, types, types_map, nframes, natom)
           
"""
    else:    
        for i in atoms:
            if i in types_map:
                continue
            else:
                types_map[i] = index
                index = index+1
        for i in atoms:
            types.append(types_map[i])
            # numbers.append(atoms_dict[i])
""" 


   


# This functions extracts forces from a CP2K trajectory forces xyz file
def extract_cp2k_frc(file, time, natom):
    uconv = 51.42208619083232  # from Ha/Bohr to eV/Å
    print(time[0])
    print("Extracting Forces from {}".format(file))
    with open(file) as ifile:
        lines = ifile.readlines()
        nframes = len(lines)//(natom+2)

        # group val in an array for each frame
        frc = [[float(lines[k*(natom+2)+i+2].split()[l])*uconv
               for i in range(natom) for l in range(1, 4)]
               for k in range(nframes) if float(lines[k*(natom+2)+1].split()[5].rstrip(',')) in time]

        # group val in an array for each element
        frc_sch = [[[float(lines[k*(natom+2)+i+2].split()[l])*uconv
                     for l in range(1, 4)] for i in range(natom)]
                   for k in range(nframes) if float(lines[k*(natom+2)+1].split()[5].rstrip(',')) in time]

    return (frc, frc_sch)

# This functions takes a CP2K input file and a CP2K output file of
# a single point force evaluation run and extracts coordinates, forces
# and box informations



def feval(file, types_file):
    with open(file, 'r') as ifile:
        lines = ifile.readlines()
        start = False
        NormalTermination = False
        crds = []
        box = []
        types_map = []
        hartree_to_ev=27.2114
        uconv = 51.42208619083232
        frcs_found=False
        with open(types_file, 'r') as tfile:
            types_dict={}
            tlines=tfile.readlines()
            for i,line in enumerate(tlines):
                types_dict[line.split()[0]] = i
          

        for i,line in enumerate(lines):
            # Extract coordinates from CP2K files
            if "CARTESIAN COORDINATES (A.U.)" in line:
                start = False
                # print(line)
            if start:
                # print(line)
                try:
                    tmp = [float(lines[i].split()[x]) for x in range(1, 4)]
                except:
                    continue
                crds += tmp
                # print(crds,len(crds))
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                start = True
                
            if "FINAL SINGLE POINT ENERGY" in line:
                        nrg = float(line.split()[-1])*hartree_to_ev
            # Extract Forces from CP2K output file and convert it to dpdata units
            if "CARTESIAN GRADIENT" in line:
                frc_raw = [-float(
                    lines[k+3].split()[i]) for k in range(i, i+len(crds)//3) for i in range(3, 6)]
                frc = [-float(
                    lines[k+3].split()[i])*uconv for k in range(i, i+len(crds)//3) for i in range(3, 6)]
                #frc_raw=np.reshape(frc_raw,(len(crds),3))
                #frc=np.reshape(frc,(len(crds),3))
                types = [int(types_dict[lines[k+3].split()[1]]) for k in range(i, i+len(crds)//3)]
                #types_map_tmp = [str(lines2[k+3].split()[2])
                    #                for k in range(j, j+len(crds)//3)]
                types_map=[x for x in types_dict.keys()]
                #for x in types_map_tmp:
                    #   if x not in types_map:
                    #      types_map.append(x)
                frcs_found = True
            if "****ORCA TERMINATED NORMALLY****" in line:
                NormalTermination=True
        #print(len(crds))  
        return (crds, nrg, frc, types, types_map, frcs_found, box,frc_raw,NormalTermination)
      
                  

def savenpy(kind, val, traj_energy, types_map, types, dir, set, box):

    np.save(dir+set+"{}_numbers".format(kind), numbers)
    np.save(dir+set+"{}".format(kind), val)
    # np.save(dir+set+"{}_grouped".format(kind),val_sch)
    np.save(dir+set+"energy", traj_energy)

    if args.box:
        box_0 = [[args.box[0], 0., 0., 0., args.box[1], 0., 0., 0., args.box[2]]]
        box = box_0*len(traj_energy)
        np.save(dir+set+"{}".format("box"), box)
    else:
        if not box:
            print("Box not defined. May cause error in data conversion")
        np.save(dir+set+"{}".format("box"), box)

    with open(dir+"type_map.raw", 'w') as ofile:
        for i in types_map:
            ofile.write(i+"\n")

    with open(dir+"type.raw", 'w') as ofile:
        for i in types:
            ofile.write(str(i)+"\n")


def savenpy_frcs(val, dir, set):
    np.save(dir+set+"force", val)


def savefake_frcs(natoms, lgth, dir, set):
    val = np.zeros((lgth, natoms*3))
    np.save(dir+set+"force", val)


crds_traj = []
crds_sch_traj = []
frc_traj = []
frc_sch_traj = []
traj_energy_traj = []
numbers = []
time_crds = []
time_frcs = []
box_traj = args.box
frcs_found=False

#if args.feval:
for file in os.listdir("./"):
    if file.endswith(".out"):
        try:
            crds, nrg, frc, types, types_map, frcs_found, box,frc_raw,NormalTermination = feval(file,args.types_map)
            print(NormalTermination)
        except:
            print(f"Error in File {file}")
        if frcs_found:
            if NormalTermination:
                crds_traj.append(crds)
                frc_traj.append(frc)
                traj_energy_traj.append(nrg)
                box_traj.append(box) 


print(crds,frc_raw,len(frc_raw))


dir = "./dpdata/"
set = "/set.000/"
if os.path.isdir(dir):
    shutil.move(dir, "backup")
os.makedirs(dir+set)

kind = "coord"
try:
    savenpy(kind, crds_traj, traj_energy_traj,
            types_map, types, dir, set, box_traj)
except:
    print("Types Map not Found")
    exit()
    


if frcs_found:
    savenpy_frcs(frc_traj, dir, set)
else:
    print("No Forces Found. Creating fake array")
    savefake_frcs(len(types), len(crds_traj), dir, set)


# # Define a  main function that parses the command line arguments and calls the
# # main function of the script.
# def main():
#     parser = argparse.ArgumentParser(
#         description='Convert CP2K output to dpdata format')
#     parser.add_argument('-feval', action='store_true',
#                         help='If set, read from FEVAL output')
#     parser.add_argument('-box', nargs=3, type=float,
#                         help='Box dimensions in Angstrom')
#     args = parser.parse_args()
#     convert_cp2k_to_dpdata(args) í
