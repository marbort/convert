import json
import os
from re import T
import numpy as np
import argparse
import shutil
import MDAnalysis as mda

parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--offset', dest='offset', default=1,
                    type=int, help='save only every nth frame')
parser.add_argument('--feval', dest='feval', action='store_true',
                    help='Extract forces in Force Eval format')
parser.add_argument('--convert', dest='convert', action='store_true',
                    help='Convert dcd to xyz')
parser.add_argument('--topo', dest='topo', type=str,
                    help='Topology file for dcd conversion')
parser.add_argument('--enerfromforces', dest='ener_frc', action='store_true',
                    help='Extract energy from forces file')


parser.add_argument('--box', dest='box', type=float,
                    help='XYZ Box Sizes (orthrombic only)', nargs=3)
args = parser.parse_args()

hartree_to_ev = 27.211399
# print(re.search(r'(?<=Element).*?(?=SUM OF)',lines,re.DOTALL)[0])

# This functions extracts coordinates and atom types from a CP2K trajectory
# positions xyz file

def convert_dcd_xyz(topo,file):
    #convert a dcd file to xyz
    u = mda.Universe(topo,file)
    with mda.Writer(file.replace(".dcd", ".xyz"), n_atoms=u.atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)  # write the coordinates of the current frame to the output file
    return file.replace(".dcd", ".xyz")

    

def extract_cp2k_coords(file, kind,ener_from_forces):
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

    # group val in an array for each frame
    val = [[float(lines[args.offset*k*(natom+2)+i+2].split()[l])
           for i in range(natom) for l in range(1, 4)]
           for k in range(nframes//args.offset)]

    # group val in an array for each element
    val_sch = [[[float(lines[args.offset*k*(natom+2)+i+2].split()[l])
                 for l in range(1, 4)] for i in range(natom)]
               for k in range(nframes//args.offset)]

    index = 0
    for i in atoms:
        if i in types_map:
            continue
        else:
            types_map[i] = index
            index = index+1
    for i in atoms:
        types.append(types_map[i])
        # numbers.append(atoms_dict[i])
    
    if ener_from_forces:
            print(f"Ciao {types_map}")
            return(numbers,  val, val_sch, types, types_map, nframes, natom)
        
    else:
        print(f"Ciao {types_map}")
        traj_ener = [np.append(traj_ener, float(lines[args.offset*k*(natom+2)+1].split()[8].rstrip())*hartree_to_ev)
                 for k in range(nframes//args.offset)]

        time = [float(lines[args.offset*k*(natom+2)+1].split()[5].strip(','))
                for k in range(max(nframes//args.offset,1))]
        print(time[0])
        return (numbers, time, val, val_sch, traj_ener, types, types_map, nframes, natom)


# This functions extracts forces from a CP2K trajectory forces xyz file
def extract_cp2k_frc(file, time, natom,ener_from_forces):
    print(f"Ener from forces: {ener_from_forces}")
    uconv = 51.42208619083232  # from Ha/Bohr to eV/Å
    print("Extracting Forces from {}".format(file))
    with open(file) as ifile:
        lines = ifile.readlines()
        nframes = len(lines)//(natom+2)
    nframes = len(lines)//(natom+2)
    traj_ener = np.empty(0)
    if ener_from_forces:
        traj_ener = [np.append(traj_ener, float(lines[args.offset*k*(natom+2)+1].split()[8].rstrip())*hartree_to_ev)
                 for k in range(nframes//args.offset)]
        time = [float(lines[args.offset*k*(natom+2)+1].split()[5].strip(','))
                for k in range(max(nframes//args.offset,1))]
        print(f"Time0: {time[0]}")
    

    # group val in an array for each frame
    frc = [[float(lines[k*(natom+2)+i+2].split()[l])*uconv
            for i in range(natom) for l in range(1, 4)]
            for k in range(nframes) if float(lines[k*(natom+2)+1].split()[5].rstrip(',')) in time]

    # group val in an array for each element
    frc_sch = [[[float(lines[k*(natom+2)+i+2].split()[l])*uconv
                    for l in range(1, 4)] for i in range(natom)]
                for k in range(nframes) if float(lines[k*(natom+2)+1].split()[5].rstrip(',')) in time]
    if ener_from_forces:
        
        return (frc, frc_sch,time,traj_ener)
    else:
        return (frc, frc_sch)

# This functions takes a CP2K input file and a CP2K output file of
# a single point force evaluation run and extracts coordinates, forces
# and box informations



def feval(file):
    with open(file, 'r') as ifile:
        lines = ifile.readlines()
        start = False
        startb = False
        crds = []
        box = []
        types_map = []

        for line in lines:
            # Extract coordinates from CP2K files
            if "&END COORD" in line:
                start = False
                # print(line)
            if start:
                # print(line)
                tmp = [float(line.split()[x]) for x in range(1, 4)]
                crds = crds+tmp
                if line.split()[0] not in types_map:
                    types_map.append(line.split()[0])
                # print(crds,len(crds))
            if "&COORD" in line:
                start = True

            # Extract Box from CP2K input
            if "&END CELL" in line:
                startb = False
            if startb:
                tmp = [float(line.split()[x]) for x in range(1, 4)]
                box = box+tmp
            if "&CELL" in line:
                startb = True
        try:
            frcname = ifile.name.replace(".inp", ".out")
            uconv = 51.42208619083232
            with open(frcname, 'r') as frcfile:
                lines2 = frcfile.readlines()
                for j, line in enumerate(lines2):
                    # Extract Energy from CP2K output and convert it to eV
                    if "Total FORCE_EVAL" in line:
                        nrg = float(line.split()[-1])*hartree_to_ev
                    # Extract Forces from CP2K output file and convert it to dpdata units
                    if "ATOMIC FORCES in" in line:
                        frc = [float(
                            lines2[k+3].split()[i])*uconv for k in range(j, j+len(crds)//3) for i in range(3, 6)]
                        types = [int(lines2[k+3].split()[1]) -
                                 1 for k in range(j, j+len(crds)//3)]
                        """
                        
                        
                       
                        types_map_tmp = [str(lines2[k+3].split()[2])
                                         for k in range(j, j+len(crds)//3)]
                        for x in types_map_tmp:
                            if x not in types_map:
                                types_map.append(x)
                        """
                        frcs_found = True
            return (crds, nrg, frc, types, types_map, frcs_found, box)
        except:
            print("Forces file not found for {}".format(ifile.name))
            frcs_found = False
            exit()


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

    with open(os.path.join(dir,"type_map.raw"), 'w') as ofile:
        for i in types_map:
            ofile.write(i+"\n")

    with open(os.path.join(dir,"type.raw"), 'w') as ofile:
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
box_traj = []

if args.convert:
    for file in os.listdir("./"):
        if file.endswith(".dcd"):
            try:
                file = convert_dcd_xyz(args.topo,file)
            except:
                print("Error converting file {}".format(file))
                continue
            print("Converted file {}".format(file))

if args.feval:
    for file in os.listdir("./"):
        if file.endswith(".inp"):
            try:
                crds, nrg, frc, types, types_map, frcs_found, box = feval(file)
            except:
                continue
            if frcs_found:
                crds_traj.append(crds)
                frc_traj.append(frc)
                traj_energy_traj.append(nrg)
                box_traj.append(box)
    dir = "./dpdata/"
    set = "/set.000/"
    if os.path.isdir(dir):
        shutil.move(dir, "backup")
    os.makedirs(dir+set)

    kind = "coord"
    savenpy(kind, crds_traj, traj_energy_traj,
            types_map, types, dir, set, box_traj)


    if frcs_found:
        savenpy_frcs(frc_traj, dir, set)
    else:
        print("No Forces Found. Creating fake array")
        savefake_frcs(len(types), len(crds_traj), dir, set)


else:
    time=[]
    for file in os.listdir("./"):
        if file.endswith(".xyz"):
            crds_traj = []
            crds_sch_traj = []
            frc_traj = []
            frc_sch_traj = []
            traj_energy_traj = []
            numbers = []
            time_crds = []
            time_frcs = []
            box_traj = []

            with open(file, 'r') as ifile:
                if "pos" in ifile.name:
                    print("Extracting coordinates from {}".format(ifile.name))
                    kind = "crds"
                    frc_file = ifile.name.replace('pos', 'frc')
                    cell_file = ifile.name.replace('-pos', '').replace('.xyz', '.cell')
                    if args.ener_frc:
                        numbers, val, val_sch, types, types_map, nframes, natom = extract_cp2k_coords(
                        ifile, kind,args.ener_frc)
                    else:
                        numbers, time, val, val_sch, traj_ener, types, types_map, nframes, natom = extract_cp2k_coords(
                        ifile, kind,args.ener_frc)
                        
                    
                    try:
                        if args.ener_frc:
                                frc, frc_sch,time, traj_ener = extract_cp2k_frc(frc_file, time, natom,args.ener_frc)
                        else:
                                frc, frc_sch = extract_cp2k_frc(frc_file, time, natom,args.ener_frc)
                        frc_traj = frc_traj+frc
                        frc_sch_traj = frc_sch_traj+frc_sch
                        frcs_found = True
                    except:
                        print("No Forces Found. Skipping extraction")
                        frcs_found = False
                   
                    
                    crds_traj = crds_traj+val
                    crds_sch_traj = crds_sch_traj+val_sch
                    time_crds = time_crds+time
                    traj_energy_traj = traj_energy_traj+traj_ener
                   
                    
                    if not args.box:
                        
                        # Extract Cell from CP2K output file and convert it to dpdata units
                        # If does not exist, exit the program
                        if os.path.isfile(cell_file):
                            
                            data = np.loadtxt(cell_file)
                        else:
                            print("Cell file not found. Exiting")
                            exit()
                            
                        try:
                             
                            box_single=data[:,2:-1]
                            time_box=data[:,1]
                            # pick indices of box that match the time of the coordinates
                            box_single=box_single[np.in1d(time_box,time)]
                            box_traj+=box_single.tolist()

                        except:
                            print("No Cell Found. Skipping extraction")

                    dir = file.replace(".xyz", "_dpdata")
                    print(types_map)
                    set = "/set.000/"
                    if os.path.isdir(dir):
                        shutil.move(dir, f"backup_{dir}")
                    os.makedirs(dir+set)

                    kind = "coord"
                    print(f"Types map: {types_map}")
                    savenpy(kind, crds_traj, traj_energy_traj,
                            types_map, types, dir, set, box_traj)


                    if frcs_found:
                        savenpy_frcs(frc_traj, dir, set)
                    else:
                        print("No Forces Found. Creating fake array")
                        savefake_frcs(len(types), len(crds_traj), dir, set)

#print(dir)
#print(f"Types map: {types_map}")

                        


"""

dir = "./dpdata/"
set = "/set.000/"
if os.path.isdir(dir):
    shutil.move(dir, "backup")
os.makedirs(dir+set)

kind = "coord"
savenpy(kind, crds_traj, traj_energy_traj,
        types_map, types, dir, set, box_traj)


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
#     convert_cp2k_to_dpdata(args)
"""