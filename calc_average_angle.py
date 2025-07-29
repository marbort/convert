import MDAnalysis as mda
import itertools
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import json

def calc_avg_angle(topol, traj, atom1,atom2,atom3, output,atom_style="default", mda_format="auto"):
    """
    Calculate the average angle between atom1, atom2, and atom3 (atom2 is the central atom in the angle).
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    atom1 (int): Index of the first atom in the angle.
    atom2 (int): Index of the second atom in the angle (central atom).
    atom3 (int): Index of the third atom in the angle.
    output (str): Path to save the average angle.
    """
    u = mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")
    
    angles = []
    for ts in u.trajectory:
        if len(u.atoms) < 3:
            continue
        a1 = u.atoms[atom1]
        a2 = u.atoms[atom2]
        a3 = u.atoms[atom3]
        angle = mda.lib.distances.calc_angles(a2.position, a1.position, a3.position)
       
        angles.append(angle)

    avg_angle = sum(angles) / len(angles) if angles else 0.0
    std_angle = np.std(angles) if angles else 0.0
    print(avg_angle, std_angle)
    #print(np.max(angles), np.min(angles), len(angles))
    hist= np.histogram(angles, bins=180, range=(0, np.pi))
    np.savetxt(f"{output}_hist.txt", np.column_stack((hist[1][:-1], hist[0])), header="Angle (radians) Frequency", fmt='%.6f %.0f') 
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(np.array(angles)*57.2958, bins=180, range=(0, np.pi),)
    ax.set_xlabel('Angle (degrees)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Histogram of angles between atoms {atom1}, {atom2}, and {atom3}')
    plt.savefig(f"{output}_histogram.png")
    plt.close(fig)
    
                 
    with open(output, 'w') as f:
        f.write(f"Average angle: {avg_angle*57.2958:.3f} Standard deviation: {std_angle*57.2958:.3f}\n")
    
    print(f"Average angle calculated and saved to {output}.")
    
def calc_avg_angle_type(topol, traj, atom1, type, max_dist, output, atom_style="default", mda_format="auto"):
    #fig=plt.figure(figsize=(48,16),dpi=150,constrained_layout=True)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 60}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    print("##########")
    
    u = mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")
    angles = {}
    pairs_length = []
    frames_check=0
    centre= u.atoms[atom1]

    for ts in u.trajectory:
        extremes=u.select_atoms(f"type {type} and around {max_dist} index {atom1}")
        if len(extremes) < 2:
            print(f"Skipping frame {ts.frame} due to insufficient number of atoms of type {type} within {max_dist} of atom {atom1}.")
            frames_check += 1
        else:
            print("Going this way")
            temp= []
            pairs = itertools.combinations(extremes, 2)
            print(list(pairs))
            if len(extremes) not in list(angles.keys()):
                angles[len(extremes)] = {"data": [], "avg": 0.0, "std": 0.0, "count": 0}
            if len(list(pairs)) not in pairs_length:
                pairs_length.append(len(list(pairs)))
        for a1,a2 in pairs:
            print("Another way")
            print(a1)
            angle = mda.lib.distances.calc_angles(a1.position, centre.position, a2.position)
            print(f"Angle between {a1.index} and {a2.index} with centre {centre.index}: {angle*57.2958:.3f} degrees")
            #if angles[len(extremes)][f"data_{j}"] == None:
            #    angles[len(extremes)][f"data_{j}"] = []
            #    angles[len(extremes)][f"data_{j}"].append(angle)
            #else:
            #    angles[len(extremes)][f"data_{j}"].append(angle)
            if angle is not None:
                temp.append(angle)
            print(f"TEMP: {temp}")
            angles[len(extremes)]["data"].append(np.average(temp)) if temp else None
    pairs_length = np.sort(pairs_length)
    with open(f"{output}_angles.json", 'w') as f:
        json.dump(angles, f)
    
    for i,key in enumerate(angles.keys()):
        try:
            angles[key]["avg"] = sum(angles[key]["data"]) / len(angles[key]["data"])
        except ZeroDivisionError:
            angles[key]["avg"] = 0.0
        angles[key]["std"] = np.std(angles[key]["data"])
        angles[key]["count"] = len(angles[key]["data"])
        frames_check += len(angles[key]["data"])
        #for j in range(pairs_length[i]):
            #ax = plt.subplot(len(pairs_length[i]), len(angles.keys()), j + 1)
            #ax.hist(np.array(angles[key]["data"])*57.2958, bins=180, range=(0, 180), alpha=0.5)
            #ax.set_xlabel('Angle (degrees)')
            #ax.set_ylabel('Frequency')
            #ax.set_title(f'Histogram of angles with {key} atoms')
    #plt.tight_layout()
    #plt.savefig(f"{output}_histogram.png")
    #plt.close(fig)
            
            
            
    #assert frames_check == len(u.trajectory), f"WRONG NUMBER OF FRAMES CHECKED: total frames {len(u.trajectory)} != frames checked {frames_check}"
    print(f"Total frames checked: {frames_check}")  
    print(f"Total frames {len(u.trajectory)}")
    
    with open(output, 'w') as f:
        f.write(f"Average angles for atom {atom1} with type {type} within {max_dist}:\n")
        for key in sorted(angles.keys()):
            f.write(f"Number of atoms: {key}, Average angle: {angles[key]['avg']*57.2958:.3f}, Standard deviation: {angles[key]['std']*57.2958:.3f}, Count: {angles[key]['count']}\n")
            
    
    
def main():
    parser = argparse.ArgumentParser(description="Calculate average angle between three atoms or average angle with a specific type of atom.")
    parser.add_argument("topol", type=str, help="Path to the topology file.")
    parser.add_argument("traj", type=str, help="Path to the trajectory file.")
    parser.add_argument("atom1", type=int, help="Index of the first atom in the angle or the central atom for type calculation.")
    parser.add_argument("--atom2", type=int, nargs=1, default=None, help="Index of the second atom in the angle (optional for type calculation).")
    parser.add_argument("--atom3", type=int, nargs=1, default=None, help="Index of the third atom in the angle (optional for type calculation).")
    parser.add_argument("--type", type=str, default=None, help="Type of atom to calculate average angle with (optional for three-atom angle calculation).")
    parser.add_argument("--max_dist", type=float, default=5.0, help="Maximum distance to consider for type calculation (default: 5.0).")
    parser.add_argument("--output", type=str, default="average_angle.txt", help="Path to save the average angle (default: average_angle.txt).")
    parser.add_argument("--atom_style", type=str, default="default", help="MDAnalysis atom style (default: default).")
    parser.add_argument("--mda_format", type=str, default="auto", help="MDAnalysis format (default: auto).")
    args = parser.parse_args()
    
    if args.type is not None:
        if args.atom2 is not None or args.atom3 is not None:
            print("Warning: --type option is used, atom2 and atom3 will be ignored.")
        print("CIAO")
        calc_avg_angle_type(args.topol, args.traj, args.atom1, args.type, args.max_dist, args.output, args.atom_style, args.mda_format)
    elif args.atom2 is not None and args.atom3 is not None:
        calc_avg_angle(args.topol, args.traj, args.atom1, args.atom2[0], args.atom3[0], args.output, args.atom_style, args.mda_format)
    else:
        print("Error: Please provide either --type with atom1 or both atom2 and atom3 for three-atom angle calculation.")
        return
if __name__ == "__main__":
    main()