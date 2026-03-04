import MDAnalysis as mda
import itertools
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import json
import math

def calc_avg_bond(topol, traj, atom1,atom2, box, output,atom_style="default", mda_format="auto"):
    """
    Calculate the average bond between atom1, atom2.
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    atom1 (int): Index of the first atom in the bond.
    atom2 (int): Index of the second atom in the bond (central atom).
    output (str): Path to save the average bond.
    """
    
    u = mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")
    
    bonds = []
    for ts in u.trajectory:
        if len(u.atoms) < 3:
            continue
        a1 = u.atoms[atom1]
        a2 = u.atoms[atom2]
        
        bond = mda.lib.distances.calc_bonds(a2.position, a1.position, box=[box,box,box,90,90,90])
       
        bonds.append(bond)

    avg_bond = sum(bonds) / len(bonds) if bonds else 0.0
    std_bond = np.std(bonds) if bonds else 0.0
#    print(avg_bond, std_bond)
    #print(np.max(bonds), np.min(bonds), len(bonds))
    hist= np.histogram(bonds, bins=180, range=(0, np.pi))
    np.savetxt(f"{output}_hist.txt", np.column_stack((hist[1][:-1], hist[0])), header="bond (Angstroem) Frequency", fmt='%.6f %.0f') 
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(np.array(bonds), bins=np.max(bonds)//0.1, range=(0, np.max(bonds)),)
    ax.set_xlabel('bond (Angstrom)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Histogram of bonds between atoms {atom1}, {atom2}')
    plt.savefig(f"{output}_histogram.png")
    plt.close(fig)
    
                 
    with open(output, 'w') as f:
        f.write(f"Average bond: {avg_bond:.3f} Standard deviation: {std_bond:.3f}\n")
    
    print(f"Average bond calculated and saved to {output}.")
    
def calc_avg_bond_type(topol, traj, atom1, atom2, type, max_dist, output, box, atom_style="default", mda_format="auto"):
    #fig=plt.figure(figsize=(48,16),dpi=150,constrained_layout=True)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 60}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    print("##########")

    types_elm = {"1": "C", "2": "O", "3": "H", "4": "Mg", "5": "Cl", "6": "C"}
    
    u = mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")
    bonds = {}
    pairs_length = []
    frames_check=0
    centre= u.atoms[atom1]
    
    

    for ts in u.trajectory:
        extremes=u.select_atoms(f"type {type} and around {max_dist} index {atom1}")
        if len(extremes) < 1:
            print(f"Skipping frame {ts.frame} due to insufficient number of atoms of type {type} within {max_dist} of atom {atom1}.")
            frames_check += 1
        else:
            temp = []
            temp2 = []
            pairs = [(atom1, t.index) for t in extremes]
            pairs_length.append(len(extremes))
            #print(list(pairs))
            if len(extremes) not in list(bonds.keys()):
                bonds[len(extremes)] = {"data": {"type": [], "atoms": []}, "avg": {"type": 0.0, "atoms": 0.0}, "std": {"type": 0.0, "atoms": 0.0},"count": {"type": 0, "atoms": 0}}
            
        
            for t1,t2 in pairs:
 #               print(f"Calculating bond for frame {ts.frame} with {len(extremes)} atoms of type {type} within {max_dist} of atom {atom1}.")
                bond_t = mda.lib.distances.calc_bonds(u.atoms[t1].position, u.atoms[t2].position,[box,box,box,90,90,90])

  #              print(f"bond between {t1.index} and {t2.index} with centre {centre.index}: {bond_t*57.2958:.3f} degrees")
                #if bonds[len(extremes)][f"data_{j}"] == None:
                #    bonds[len(extremes)][f"data_{j}"] = []
                #    bonds[len(extremes)][f"data_{j}"].append(bond)
                #else:
                #    bonds[len(extremes)][f"data_{j}"].append(bond)
                if bond_t is not None:
                    bonds[len(extremes)]["data"]["type"].append(bond_t) 
    pairs_length = np.sort(np.unique(pairs_length))
    
    
    for i,key in enumerate(bonds.keys()):
        bonds[key]["avg"]["type"] = np.average(bonds[key]["data"]["type"]) 
        bonds[key]["avg"]["atoms"] = np.average(bonds[key]["data"]["atoms"])
        bonds[key]["std"]["type"] = np.std(bonds[key]["data"]["type"])
        bonds[key]["std"]["atoms"] = np.std(bonds[key]["data"]["atoms"])
        bonds[key]["count"]["type"] = len(bonds[key]["data"]["type"])
        bonds[key]["count"]["atoms"] = len(bonds[key]["data"]["atoms"])
        #assert  bonds[key]["count"]["type"]//math.comb(int(key),2)== bonds[key]["count"]["atoms"], f"WRONG NUMBER OF bonds CHECKED: type {bonds[key]['count']['type']} != atoms {bonds[key]['count']['atoms']}"
        frames_check += len(bonds[key]["data"])
        #for j in range(pairs_length[i]):
            #ax = plt.subplot(len(pairs_length[i]), len(bonds.keys()), j + 1)
            #ax.hist(np.array(bonds[key]["data"])*57.2958, bins=180, range=(0, 180), alpha=0.5)
            #ax.set_xlabel('bond (degrees)')
            #ax.set_ylabel('Frequency')
            #ax.set_title(f'Histogram of bonds with {key} atoms')
    #plt.tight_layout()
    #plt.savefig(f"{output}_histogram.png")
    #plt.close(fig)
    with open(f"{output}_bonds.json", 'w') as f:
        json.dump(bonds, f,indent=4)
            
            
    
    print(f"Total frames checked: {frames_check}")  
    print(f"Total frames {len(u.trajectory)}")
    
    with open(output, 'w') as f:
        f.write(f"Average bonds for atom {atom1} with type {type} within {max_dist}:\n")
    
        for key in sorted(bonds.keys()):
            f.write(f"Number of atoms: {key}, Average bond {types_elm[type]}_Mg_{types_elm[type]}: {bonds[key]['avg']['type']:.3f}, Standard deviation : {bonds[key]['std']['type']:.3f}, Count: {bonds[key]['count']['type']}\n")
            f.write(f"Number of atoms: {key}, Average bond A2_A1: {bonds[key]['avg']['atoms']:.3f}, Standard deviation : {bonds[key]['std']['atoms']:.3f}, Count: {bonds[key]['count']['atoms']}\n")
            f.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate average bond between three atoms or average bond with a specific type of atom.")
    parser.add_argument("topol", type=str, help="Path to the topology file.")
    parser.add_argument("traj", type=str, help="Path to the trajectory file.")
    parser.add_argument("atom1", type=int, help="Index of the first atom in the bond or the central atom for type calculation.")
    parser.add_argument("--atom2", type=int, nargs=1, default=None, help="Index of the second atom in the bond (optional for type calculation).")
   
    parser.add_argument("--type", type=str, default=None, help="Type of atom to calculate average bond with (optional for three-atom bond calculation).")
    parser.add_argument("--max_dist", type=float, default=5.0, help="Maximum distance to consider for type calculation (default: 5.0).")
    parser.add_argument("--output", type=str, default="average_bond.txt", help="Path to save the average bond (default: average_bond.txt).")
    parser.add_argument("--atom_style", type=str, default="default", help="MDAnalysis atom style (default: default).")
    parser.add_argument("--mda_format", type=str, default="auto", help="MDAnalysis format (default: auto).")
    parser.add_argument("--box", type=float, default=1, help="Box side length")
    args = parser.parse_args()
    
    if args.type is not None:
            calc_avg_bond_type(args.topol, args.traj, args.atom1, args.atom2,  args.type, args.max_dist, args.output, args.box, args.atom_style, args.mda_format)
    elif args.atom2 is not None:
        calc_avg_bond(args.topol, args.traj, args.atom1, args.atom2[0], args.box, args.output, args.atom_style, args.mda_format)
    else:
        print("Error: Please provide either --type with atom1 or  atom2  for two-atom bond calculation.")
        return
if __name__ == "__main__":
    main()
