import os
import argparse

def check_explosion(lammps_log):
    with open(lammps_log,'r') as ifile:
        lines=ifile.readlines()
    box_0=[]
    box=[]
    for i,line in enumerate(lines):
        spt=line.split()
        if "Lx" in spt:
            idx=spt.index("Lx")
            box_0=[float(lines[i+1].split()[idx]),float(lines[i+1].split()[idx+1]),float(lines[i+1].split()[idx+2])]
        if "Loop" in spt:
            box=[float(lines[i-1].split()[idx]),float(lines[i-1].split()[idx+1]),float(lines[i-1].split()[idx+2])]

    diff=[x-box_0[i] for i,x in enumerate(box)]
    for i,item in enumerate(diff):
        explosion=False
        if abs(item)/box_0[i] >= 0.5:
            explosion=True
    return(explosion)

def main():
    parser = argparse.ArgumentParser(
                    prog='CheckExplosions',
                    description='Check LAMMPS MD explosions',
                    epilog='Hope you are safe')

  
    parser.add_argument('--input', type=str, help='File path to LAMMPS log file.')

    args = parser.parse_args()
    
    check_explosion(args.input)


if __name__ == "__main__":
    main()

            
            
            