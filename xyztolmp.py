import numpy as np
import sys
import dpdata as dp
import argparse


type_map={"Turbo":{"C":{"Index":0,"Mass":12.0107},"O":{"Index":1,"Mass":15.9994},"H":{"Index":2,"Mass":1.00794},
                   "Mg":{"Index":3,"Mass":24.305},"Cl":{"Index":4,"Mass":35.453},"Li":{"Index":5,"Mass":6.941}},
          "Normal":{"C":1,"O":2,"H":3,"Mg":4,"Cl":5},
          "Hauser":{"C":1,"O":2,"H":3,"Cl":4,"Li":5,"Mg":6, "N":7},
          "Turbo_Marinella":{"C":{"Index":0,"Mass":12.0107},"O":{"Index":1,"Mass":15.9994},"H":{"Index":2,"Mass":1.00794},
                   "Mg":{"Index":3,"Mass":24.305},"Cl":{"Index":4,"Mass":35.453},"Li":{"Index":5,"Mass":6.941},"Br":{"Index":6,"Mass":79.904}, "I":{"Index":7,"Mass":126.904}},
          "DES_new_AT":{
            "C":{"Index":0,"Mass":12.0107},
            "O":{"Index":1,"Mass":15.9994},
            "H":{"Index":2,"Mass":1.00794},
            "Mg":{"Index":3,"Mass":24.305},
            "Cl":{"Index":4,"Mass":35.453},
            "CR":{"Index":5,"Mass":12.0107},
            "N":{"Index":6,"Mass":14.006}}
            
}

#type_map={"C":1,"Cl":2,"H":3,"N":4,"O":5,"Mg":6}


"""OLD VERSION
elm,x,y,z=np.loadtxt(sys.argv[1],unpack=True,skiprows=2,dtype=str)
types=[type_map[sys.argv[2]][x] for x in elm]

with open('conf_data_temp','w') as ofile:
    for i,line in enumerate(types):
        ofile.write("{:8d} {:8d} {:10.6f} {:10.6f} {:10.6f}\n".format(i+1,line,float(x[i]),float(y[i]),float(z[i])))





print('ciao')
"""

def create_mass_section():
    mass_section="\nMasses\n\n"
    added=[]
    for i,item in type_map[args.type_set].items():
        if type(item) == dict:
            index=item["Index"]
            mass=item["Mass"]
        else:
            index=item
            mass=0.0
        if index not in added:
            mass_section+=f"{index+1} {mass:.4f} \n"
            added.append(index)
    mass_section+="\n"
    return mass_section

#NEW VERSION
parser = argparse.ArgumentParser(description='Convert XYZ to LAMMPS data file')
parser.add_argument('input', type=str, help='Input XYZ file')
parser.add_argument('type_set', type=str, help='Type set to use (Turbo, Normal, Hauser, Turbo_Marinella)')
parser.add_argument('--cell', type=float, default=None, help='cell size (low high)',nargs=2)
parser.add_argument('--triclinic', dest='triclinic',action='store_true',help='Set if the box is triclinic')
args = parser.parse_args()



data=dp.System(args.input, 'xyz')
for i,item in enumerate(data['atom_types']):
    data['atom_types'][i]=type_map[sys.argv[2]][data['atom_names'][int(item)]]["Index"]

    
    

data[-1].to('lammps/lmp','temp.data')
newlines=[]
with open(f"temp.data",'r') as ifile:
    lines=ifile.readlines()
    for line in lines:
        if 'yz' in line:
            if args.triclinic:
                newlines.append(line)
                newlines.append(create_mass_section())
            else:
                newlines.append(create_mass_section())
        elif 'atom types' in line:
            newlines.append(f"{len(type_map[args.type_set])} atom types\n")
        elif 'xlo xhi' in line and args.cell is not None:
            newlines.append(f"{args.cell[0]:15.10f} {args.cell[1]:15.10f} xlo xhi\n")
        elif 'ylo yhi' in line and args.cell is not None:
            newlines.append(f"{args.cell[0]:15.10f} {args.cell[1]:15.10f} ylo yhi\n")
        elif 'zlo zhi' in line and args.cell is not None:
            newlines.append(f"{args.cell[0]:15.10f} {args.cell[1]:15.10f} zlo zhi\n")
        else:
            newlines.append(line)

with open(f"{args.input.split('.')[0]}.data",'w') as ofile:
    for line in newlines:
        ofile.write(line)
