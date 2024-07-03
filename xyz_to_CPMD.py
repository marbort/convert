import sys
import numpy


"""
Example CPMD atom syntax

&ATOMS
*H_LDA.psp
 LMAX=S
  2
 0.371   0.000   0.000
-0.371   0.000   0.000
&END


"""

def extract_elements(input):
    angtobohr=1.89
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    elems={}
    tot_a=int(lines[0])
    for line in lines[2:]:
        if line.split()[0] in elems:
            elems[line.split()[0]].append([float(x)*angtobohr for x in line.split()[1:]])
        else:
                elems[line.split()[0]]=[]
                elems[line.split()[0]].append([float(x)*angtobohr for x in line.split()[1:]])
    tot_extracted=0
    for elm in elems:
        tot_extracted+=len(elems[elm])
    
    if tot_extracted != tot_a:
        print("ERROR")
        print(tot_a,tot_extracted)
    else:
        return(elems)

def write_atoms_CPMD(elems,pseudofile):
    pseudo={}
    with open(pseudofile,'r') as ifile:
        lines=ifile.readlines()
    for line in lines:
        split=line.split(' - ')
        pseudo[split[0].split()[0]]=[split[0].split()[1]]+[split[x] for x in range(1,len(split))]
        
    with open('atoms_CPMD.dat','w') as ofile:
        ofile.write('&ATOMS\n')
        for elm in elems:
            ofile.write(f"*{pseudo[elm][0]} {pseudo[elm][1]}\n")
            ofile.write(f"LMAX={pseudo[elm][2]}\n")
            ofile.write(f"{len(elems[elm])}\n")
            for crd in elems[elm]:
                ofile.write("  ".join([f"{x:.6f}" for x in crd])+"\n")
            ofile.write("\n")
        ofile.write("&END\n")
                
                
                



elems=extract_elements(sys.argv[1])
pseudofile=sys.argv[2]
write_atoms_CPMD(elems,pseudofile)