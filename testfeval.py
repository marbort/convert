import json
import os
from re import T
import numpy as np
import argparse
import shutil

def feval(file, types_file):
    print(types_file)
    with open(file, 'r') as ifile:
        lines = ifile.readlines()
        start = False
        startb = False
        crds = []
        box = []
        types_map = []
        hartree_to_ev=1
        with open(types_file, 'r') as tfile:
            types_dict={}
            tlines=tfile.readlines()
            for i,line in enumerate(tlines):
                types_dict[line.split()[0]] = i
            print(types_dict)

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
                        types = [int(types_dict[lines2[k+3].split()[2]]) for k in range(j, j+len(crds)//3)]
                        #types_map_tmp = [str(lines2[k+3].split()[2])
                         #                for k in range(j, j+len(crds)//3)]
                        types_map=[x for x in types_dict.keys()]
                        #for x in types_map_tmp:
                         #   if x not in types_map:
                          #      types_map.append(x)
                        frcs_found = True
                print(types)
            return (crds, nrg, frc, types, types_map, frcs_found, box)
        except:
            print("Forces file not found for {}".format(ifile.name))
            frcs_found = False
            exit()
            


crds, nrg, frc, types, types_map, frcs_found, box=feval('cp2k_It1_000023.inp','./type_map.raw')
print(types_map)
print(types)
    