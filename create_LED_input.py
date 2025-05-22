import argparse

parser = argparse.ArgumentParser(description="Plot data")
parser.add_argument("--input", dest="input", type=str, help="input data")
parser.add_argument(
    "--led1", dest="led1", type=int, help="indexes (1-based) of atoms in fragment 1 for LED", default=None, nargs="+"
    )
parser.add_argument(
    "--led2", dest="led2", type=int, help="indexes (1-based) of atoms in fragment 2 for LED", default=None, nargs="+"
    )
args = parser.parse_args()

if args.led1:
    led1=[]
    led2=[]
    if not args.led2:
        raise ValueError("led2 is not defined")
    
    else:
        with open(args.input, 'r') as f:
            lines = f.readlines()
        if len(args.led1) + len(args.led2) != len(lines)-2:
            raise ValueError("led1 and led2 do not match the number of atoms in the system")
        for i in args.led1:
            if i in args.led2:
                raise ValueError(f"Atom {i} cannot be in both fragments!!")
            else:
                lines[i+1] = lines[i+1].replace(lines[i+1].split()[0], f"{lines[i+1].split()[0]}(1)")
                led1.append(lines[i+1])
        for i in args.led2:
            lines[i+1] = lines[i+1].replace(lines[i+1].split()[0], f"{lines[i+1].split()[0]}(2)")
            led2.append(lines[i+1])
    with open('orca.inp', 'w') as f:
        f.write("""
!DLPNO-CCSD(T) DEF2-SVP DEF2-SVP/C DEF2/J RIJCOSX VERYTIGHTSCF TIGHTPNO LED
%pal
nprocs 16
end


*XYZ 0 1
""")
        for line in lines[2:]:
            f.write(line)
        f.write("\n")
        f.write("*\n")
    with open('orca_F1.inp', 'w') as f:
        f.write("""
!DLPNO-CCSD(T) DEF2-SVP DEF2-SVP/C DEF2/J RIJCOSX VERYTIGHTSCF TIGHTPNO 
%pal
nprocs 16
end


*XYZ 0 1
""")
        for line in led1:
            f.write(line)
        f.write("\n")
        f.write("*\n")
    with open('orca_F2.inp', 'w') as f:
        f.write("""
!DLPNO-CCSD(T) DEF2-SVP DEF2-SVP/C DEF2/J RIJCOSX VERYTIGHTSCF TIGHTPNO 
%pal
nprocs 16
end


*XYZ 0 1\n
""")
        for line in led2:
            f.write(line)
        f.write("\n")
        f.write("*\n")