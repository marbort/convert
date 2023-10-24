import os
import argparse

def decompose_data(file,natom):
    last=0
    with open(file,'r') as ifile:
        lines=ifile.readlines()
    for i,line in enumerate(lines):
        if "time" in line:
            last=i
    reac=[x for x in lines[last+1:last+1+natom]]
    box=[x for x in lines[last+1+natom:]]
    with open("{}_react".format(file),'w') as ofile:
        for line in reac:
            ofile.write("{}".format(line))
    with open("{}_box".format(file),'w') as ofile:
        for line in box:
            ofile.write("{}".format(line))
    return(reac,box)

def compose_data(reac):
    box=reac.replace("react","box")
    with open(reac,'r') as ifile:
        lines_reac=ifile.readlines()
    with open(box,'r') as ifile:
        lines_box=ifile.readlines()
    natoms="{}\n".format(len(lines_reac)+len(lines_box)-2)
    header="i =     8800, time =     2200.000, E =     -8598.7226452004\n"
    with open(reac.replace('react','whole'),'w') as ofile:
        ofile.write(natoms)
        ofile.write(header)
        for line in lines_reac[2:]:
            ofile.write("{}".format(line))
        for line in lines_box:
            ofile.write("{}".format(line.lstrip()))
    
    
def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input' , dest='input',help='string to match for input. Default=pos',default='pos')
    parser.add_argument('--nat' , dest='nat', help='number of atoms of reactant',default=1, type=int)
    parser.add_argument('--compose' , dest='compose', action='store_true',help='merge react + box')
    args = parser.parse_args()
    
    if args.compose:
        for file in os.listdir('./'):
            if "reac" in file:
                compose_data(file)
    else:
        for file in os.listdir('./'):
            if args.input in file and file.endswith('.xyz'):
                decompose_data(file,args.nat)


if __name__=="__main__":
    main()
            
        