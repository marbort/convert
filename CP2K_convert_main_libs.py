import os
import numpy as np
import argparse
import shutil


def feval(infile,outfile):
    hartree_to_ev = 27.211386245988
    with open(infile, 'r') as ifile:
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
            frcname = outfile
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
                        types_map_tmp = [str(lines2[k+3].split()[2])
                                         for k in range(j, j+len(crds)//3)]
                        for x in types_map_tmp:
                            if x not in types_map:
                                types_map.append(x)
                        frcs_found = True
            return (crds, nrg, frc, types, types_map, frcs_found, box)
        except:
            print("Forces file not found for {}".format(ifile.name))
            frcs_found = False
            exit()


def savenpy(kind, val, traj_energy, types_map, types, dir, set, box):

    # np.save(dir+set+"{}_numbers".format(kind), '')
    np.save(dir+set+"{}".format(kind), val)
    # np.save(dir+set+"{}_grouped".format(kind),val_sch)
    np.save(dir+set+"energy", traj_energy)

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


def main():
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--offset', dest='offset', default=1,
                        type=int, help='save only every nth frame')

    args = parser.parse_args()

    crds_traj = []
    crds_sch_traj = []
    frc_traj = []
    frc_sch_traj = []
    traj_energy_traj = []
    numbers = []
    time_crds = []
    time_frcs = []
    box_traj = []

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


if __name__ == "__main__":
    main()
