import sys
import os
import numpy as np
import json

def read_pdb(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
        print(lines[30][12:16])
        typ=[]
        number=[]
        name=[]
        alt_loc=[]
        resname=[]
        chain=  []
        resnumber=[]
        rescode=[]
        x=[]
        y=[]
        z=[]
        occupancy=[]
        tempfac=[]
        segment=[]
        element=[]
        charge=[]
        for i,line in enumerate(lines[1:-1]):
            try:
                typ.append(line[0:5])
                number.append(line[6:11])         #	Atom serial number	right	integer
                name.append(line[12:16])          #	Atom name	left*	character
                alt_loc.append(line[16])          #	Alternate location indicator		character
                resname.append(line[17:20])       #	Residue name	right	character
                chain.append(line[21])	      #  Chain identifier		character
                resnumber.append(line[22:26])     #	Residue sequence number	right	integer
                rescode.append(line[26])          #	Code for insertions of residues		character
                x.append(line[30:37])             #	X orthogonal Å coordinate	right	real (8.3)
                y.append(line[38:45])	          # Y orthogonal Å coordinate	right	real (8.3)
                z.append(line[46:55])             #	Z orthogonal Å coordinate	right	real (8.3)
                occupancy.append(line[54:59])      #	Occupancy	right	real (6.2)
                tempfac.append(line[60:65])	      # Temperature factor	right	real (6.2)
                segment.append(line[72:75])       #	Segment identifier¶	left	character
                element.append(line[76:77])	      # Element symbol	right	character
                charge.append(line[78:79])	      # Charge		character
            except:
                print(line,i)
    return(typ,number,name,alt_loc,resname,chain,resnumber,rescode,x,y,z,occupancy,tempfac,segment,
           element,charge)
    
    
def write_pqr(ff,typ,number,name,alt_loc,resname,chain,resnumber,rescode,x,y,z):
    with open(ff,'r') as fffile:
        data=json.load(fffile)
    print(data['Clm']['Clm'])
    pqrlines=[]
    for i,line in enumerate(typ):
        #try:
            pqrline=f"{line:6s}{number[i]:>5s}{name[i]:<3s}{resname[i]:>5s}{resnumber[i]:>5s}{float(x[i]):13.3f}{float(y[i]):9.3f}{float(z[i]):9.3f}\
{data[resname[i].strip()][name[i].strip()]['charge']:10.6f}{data[resname[i].strip().strip()][name[i].strip()]['vdW']:10.6f}"
            pqrlines.append(pqrline)
        #except:
        #   print(name[i])
        #   sys.exit()
    #print(pqrlines[-1])
    with open('output.pqr','w') as ofile:
        for line in pqrlines:
            ofile.write(f"{line}\n")
def check_pqr_charge(input):
    charges=[]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    for line in lines:
        charges.append(float(line[57:66]))
    return(charges,sum(charges))



def  main():
    typ,number,name,alt_loc,resname,chain,resnumber,rescode,x,y,z,occupancy,tempfac,segment,element,charge=read_pdb(sys.argv[1])
    print(f"CHL {resname.count('CHL')//21},Clm {resname.count('Clm')}, GCL {resname.count('GCL')/14} THF {resname.count('THF')/13}")
    write_pqr('/run/media/marco/T7 Shield/SHARED/GitHub/mio/convert/charges.json',
              typ,number,name,alt_loc,resname,chain,resnumber,rescode,x,y,z)
    charges,sumchrg=check_pqr_charge('output.pqr')
    print(sumchrg,str(charges[-1])[-1])



if __name__=="__main__":
    main()
    
        