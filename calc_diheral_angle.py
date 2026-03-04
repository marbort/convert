import MDAnalysis as mda
import itertools
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import json
import math

def calc_avg_dihed(topol, traj, atom1,atom2,atom3,atom4,cutoff, box, output,atom_style="default", mda_format="auto"):
    """
    Calculate the average dihed between atom1, atom2.
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    atom1 (int): Index of the first atom in the dihed.
    atom2 (int): Index of the second atom in the dihed (central atom).
    output (str): Path to save the average dihed.
    """
    
    u = mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")
    
    diheds = []
    for ts in u.trajectory:
        if len(u.atoms) < 3:
            continue
        a1 = u.atoms[atom1]
        a2 = u.atoms[atom2][0]
        a3 = u.atoms[atom3][0]
        a4 = u.atoms[atom4][0]
        #print(a1)
        #print(a2)
        #print(a3)
        #print(a4)
        dist=mda.lib.distances.calc_bonds(a2.position, a3.positionq, box=[box,box,box,90,90,90])
        
        if dist > cutoff:
            continue
        else:
        
            dihed = mda.lib.distances.calc_dihedrals(a1.position, a2.position, a3.position, a4.position, box=[box,box,box,90,90,90])
            diheds.append(dihed)

    avg_dihed = sum(diheds) / len(diheds) if diheds else 0.0
    std_dihed = np.std(diheds) if diheds else 0.0
#    print(avg_dihed, std_dihed)
    #print(np.max(diheds), np.min(diheds), len(diheds))
    hist= np.histogram(diheds, bins=360, range=(-np.pi, np.pi))
    np.savetxt(f"{output}_hist.txt", np.column_stack((hist[1][:-1], hist[0])), header="dihed (Angstroem) Frequency", fmt='%.6f %.0f') 
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(hist[1][:-1]*57.2958, hist[0]*57.2958, width=np.diff(hist[1])*57.2958, color='blue', edgecolor='black', alpha=0.7)
    #ax.hist(np.array(diheds), bins=np.max(diheds)//0.1, range=(0, np.max(diheds)),)
    ax.set_xlabel('dihed (Angstrom)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Histogram of diheds between atoms {atom1}, {atom2}, {atom3}, {atom4}')
    plt.savefig(f"{output}_histogram.png")
    plt.close(fig)
    
                 
    with open(output, 'w') as f:
        f.write(f"Average dihed: {avg_dihed:.3f} Standard deviation: {std_dihed:.3f}\n")
    
    print(f"Average dihed calculated and saved to {output}.")
    


def main():
    parser = argparse.ArgumentParser(description="Calculate average dihed between three atoms or average dihed with a specific type of atom.")
    parser.add_argument("topol", type=str, help="Path to the topology file.")
    parser.add_argument("traj", type=str, help="Path to the trajectory file.")
    parser.add_argument("--atom1", type=int, help="Index of the first atom in the dihed or the central atom for type calculation.")
    parser.add_argument("--atom2", type=int, nargs=1, default=None, help="Index of the second atom in the dihed.")
    parser.add_argument("--atom3", type=int, nargs=1, default=None, help="Index of the third atom in the dihed.")
    parser.add_argument("--atom4", type=int, nargs=1, default=None, help="Index of the fourth atom in the dihed.")
    parser.add_argument("--max_dist", type=float, default=5.0, help="Maximum distance to consider dihedral (default: 5.0).")
    parser.add_argument("--output", type=str, default="average_dihed.txt", help="Path to save the average dihed (default: average_dihed.txt).")
    parser.add_argument("--atom_style", type=str, default="default", help="MDAnalysis atom style (default: default).")
    parser.add_argument("--mda_format", type=str, default="auto", help="MDAnalysis format (default: auto).")
    parser.add_argument("--box", type=float, default=1, help="Box side length")
    args = parser.parse_args()
    
    if args.atom1 is not None and args.atom2 is not None and args.atom3 is not None and args.atom4 is not None:
            calc_avg_dihed(args.topol, args.traj, args.atom1, args.atom2,args.atom3,args.atom4,args.max_dist,args.box, args.output, args.atom_style, args.mda_format)
    else:
        print("Error: Please provide four atom indices for dihedral calculation.")
        return
if __name__ == "__main__":
    main()
