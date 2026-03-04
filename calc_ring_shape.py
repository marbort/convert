import MDAnalysis as mda
import MDAnalysis.transformations as mda_transform_translate
from MDAnalysis.transformations import TransformationBase
import MDAnalysis.transformations.translate as pluto
import MDAnalysis.transformations.wrap as pippo
import itertools
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import json
import sys
import inspect
from functools import partial

class center_in_box(TransformationBase):
    """
    Translates the coordinates of a given :class:`~MDAnalysis.coordinates.timestep.Timestep`
    instance so that the center of geometry/mass of the given :class:`~MDAnalysis.core.groups.AtomGroup`
    is centered on the unit cell. The unit cell dimensions are taken from the input Timestep object.
    If a point is given, the center of the atomgroup will be translated to this point instead.

    Example
    -------

    .. code-block:: python

        ag = u.residues[1].atoms
        ts = MDAnalysis.transformations.center(ag,center='mass')(ts)

    Parameters
    ----------
    ag: AtomGroup
        atom group to be centered on the unit cell.
    center: str, optional
        used to choose the method of centering on the given atom group. Can be 'geometry'
        or 'mass'
    point: array-like, optional
        overrides the unit cell center - the coordinates of the Timestep are translated so
        that the center of mass/geometry of the given AtomGroup is aligned to this position
        instead. Defined as an array of size 3.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup.

    Returns
    -------
    :class:`~MDAnalysis.coordinates.timestep.Timestep` object


    .. versionchanged:: 2.0.0
        The transformation was changed from a function/closure to a class
        with ``__call__``.
    .. versionchanged:: 2.0.0
       The transformation was changed to inherit from the base class for
       limiting threads and checking if it can be used in parallel analysis.
    """

    def __init__(
        self,
        ag,
        center="geometry",
        point=None,
        wrap=False,
        max_threads=None,
        parallelizable=True,
    ):
        super().__init__(
            max_threads=max_threads, parallelizable=parallelizable
        )

        self.ag = ag
        self.center = center
        self.point = point
        self.wrap = wrap

        pbc_arg = self.wrap
        if self.point:
            self.point = np.asarray(self.point, np.float32)
            if self.point.shape != (3,) and self.point.shape != (1, 3):
                raise ValueError("{} is not a valid point".format(self.point))
        try:
            if self.center == "geometry":
                self.center_method = partial(
                    self.ag.center_of_geometry, wrap=pbc_arg
                )
            elif self.center == "mass":
                self.center_method = partial(
                    self.ag.center_of_mass, wrap=pbc_arg
                )
            else:
                raise ValueError(f"{self.center} is valid for center")
        except AttributeError:
            if self.center == "mass":
                errmsg = f"{self.ag} is not an AtomGroup object with masses"
                raise AttributeError(errmsg) from None
            else:
                raise ValueError(
                    f"{self.ag} is not an AtomGroup object"
                ) from None

    def _transform(self, ts):
        if self.point is None:
            if ts.dimensions is None:
                raise ValueError("Box is None")
            boxcenter = np.sum(ts.triclinic_dimensions, axis=0) / 2
        else:
            boxcenter = self.point

        ag_center = self.center_method()

        vector = boxcenter - ag_center
        ts.positions += vector

        return ts






def get_universe(topol, traj, atom_style="default", mda_format="auto"):
    """
    Create a Universe object from the topology and trajectory files.
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    
    Returns:
    MDAnalysis Universe object.
    """
    return mda.Universe(topol, traj, atom_style=f"{atom_style}", format=f"{mda_format}")

def calc_ring_shape(universe, output, a1,a2,a3,a4, r1, r2, max_dist,box=None):
    """
    Calculate the shape of a ring defined by four atoms and their distances.
    
    Parameters:
    universe (MDAnalysis Universe): The universe containing the atoms.
    output (str): Path to save the ring shape data.
    a1, a2, a3, a4 (int): Indices of the four atoms defining the ring.
    r1, r2 (float): Distances for the ring shape calculation.
    max_dist (float): Maximum distance for filtering.
    
    Returns:
    None
    """
    angles = []
    frames=[]
    tors_R1 = []
    tors_R2 = []
    r1=universe.atoms[r1]
    r2=universe.atoms[r2]
    center=universe.select_atoms(f"index {a1}")
    
    
    for ts in universe.trajectory:   
        cog=center.center_of_geometry(compound='residues')
        shift = universe.dimensions[:3] / 2 - cog
        universe.atoms.translate(shift)
        universe.atoms.wrap(compound='atoms')
        frame= ts.frame
        frames.append(frame)
        #ts = center_in_box(center,center='geometry')(ts)
        ring=True
        if len(universe.atoms) < 4:
            continue
        positions = []
        positions.append(universe.atoms[a1].position)
        positions.append(universe.atoms[a2].position)
        positions.append(universe.atoms[a3].position)
        positions.append(universe.atoms[a4].position)
        
        
        for i,pos in enumerate(positions):
            if i == len(positions) - 1:
                next_pos = positions[0]
            else:
                next_pos = positions[i + 1]
            
            dist=mda.lib.distances.calc_bonds(pos, next_pos,box)
            if dist > max_dist:
                ring = False
        if ring:
            
            g1=universe.select_atoms(f"index {a2} or index {a4}")
            middle_point= g1.center_of_geometry(compound='group')
            angle=mda.lib.distances.calc_angles(positions[0], middle_point, positions[2],box)
            if angle < 0.8:
                print(f"Small angle detected: {angle*57.2958:.3f} degrees at frame {frame}")
            angles.append(angle)
            tors_R1.append(mda.lib.distances.calc_dihedrals(r1.position, positions[0], positions[1], positions[2],box))
            tors_R2.append(mda.lib.distances.calc_dihedrals(r2.position, positions[2], positions[3], positions[0],box))
            
    
    with open(f"{output}.json", 'w') as f:
        json.dump({
            "angles": angles,
            "tors_R1": tors_R1,
            "tors_R2": tors_R2
        }, f)
    print(f"Min angle: {min(angles) * 57.2958:.3f} degrees at frame {frames[angles.index(min(angles))]}")
    with open(output, 'w') as f:
        f.write(f"Average angle: {np.mean(angles) * 57.2958:.3f} degrees\n")
        f.write(f"Standard deviation of angles: {np.std(angles) * 57.2958:.3f} degrees\n")
        f.write(f"Average torsion R1: {np.mean(tors_R1) * 57.2958:.3f} degrees\n")
        f.write(f"Standard deviation of torsion R1: {np.std(tors_R1) * 57.2958:.3f} degrees\n")
        f.write(f"Average torsion R2: {np.mean(tors_R2) * 57.2958:.3f} degrees\n")
        f.write(f"Standard deviation of torsion R2: {np.std(tors_R2) * 57.2958:.3f} degrees\n")
    return angles, tors_R1, tors_R2

def plot_histograms(angles, tors_R1, tors_R2, output):
    #Plotting the results
    fig = plt.figure(figsize=(16, 10), dpi=150)
    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    
    hist_angles,bins_angles = np.histogram(np.array(angles) * 57.2958, bins=180,density=True)
    hist_tors_R1, bins_tors_R1 = np.histogram(np.array(tors_R1) * 57.2958, bins=180,density=True)
    hist_tors_R2, bins_tors_R2 = np.histogram(np.array(tors_R2) * 57.2958, bins=180,density=True)
    
    ax1.bar(bins_angles[:-1], hist_angles, width=np.diff(bins_angles), color='blue', alpha=0.7, edgecolor='black')
    ax2.bar(bins_tors_R1[:-1], hist_tors_R1, width=np.diff(bins_tors_R1), color='green', alpha=0.7, edgecolor='black')
    ax3.bar(bins_tors_R2[:-1], hist_tors_R2, width=np.diff(bins_tors_R2), color='orange', alpha=0.7, edgecolor='black')
    
    ax1.axvline(np.mean(angles) * 57.2958, color='red', linestyle='dashed', linewidth=1)
    ax2.axvline(np.mean(tors_R1) * 57.2958, color='red', linestyle='dashed', linewidth=1)
    ax3.axvline(np.mean(tors_R2) * 57.2958, color='red', linestyle='dashed', linewidth=1)
    
    ax1.set_xlabel('Angle (degrees)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Histogram of Angles')
    ax1.set_xlim(110, 180)
    ax1.set_ylim(0, 0.05)
    
    ax2.set_xlabel('Torsion R1 (degrees)')
    ax2.set_ylabel('Frequency')
    ax2.set_title('Histogram of Torsion R1')
    ax2.set_xlim(-180, 180)
    ax2.set_ylim(0, 0.01)
    
    ax3.set_xlabel('Torsion R2 (degrees)')
    ax3.set_ylabel('Frequency')
    ax3.set_title('Histogram of Torsion R2')
    ax3.set_xlim(-180, 180)
    ax3.set_ylim(0, 0.01)
    plt.tight_layout()
    plt.savefig(f"{output}_histograms.png", dpi=150)
    plt.close(fig)
    

def main():
    parser = argparse.ArgumentParser(description='Calculate ring shape and angles')
    parser.add_argument('--topol', type=str, required=True, help='Path to the topology file')
    parser.add_argument('--traj', type=str, required=True, help='Path to the trajectory file')
    parser.add_argument('--output', type=str, required=True, help='Output file for the ring shape data')
    parser.add_argument('--a1', type=int, required=True, help='Index of atom 1')
    parser.add_argument('--a2', type=int, required=True, help='Index of atom 2')
    parser.add_argument('--a3', type=int, required=True, help='Index of atom 3')
    parser.add_argument('--a4', type=int, required=True, help='Index of atom 4')
    parser.add_argument('--r1', type=int, required=True, help='Index of reference atom 1')
    parser.add_argument('--r2', type=int, required=True, help='Index of reference atom 2')
    parser.add_argument('--max_dist', type=float, default=5.0, help='Maximum distance for filtering (default: 5.0)')
    parser.add_argument('--atom_style', type=str, default="default", help='Atom style for MDAnalysis (default: "default")')
    parser.add_argument('--mda_format', type=str, default="auto", help='MDAnalysis format (default: "auto")')
    parser.add_argument('--box', type=float, nargs=6, default=None, help='Box dimensions and angles (optional)')
    parser.add_argument('--plotonly', type=str, default=None, help='Plot only from existing data file (optional)')
    
    args = parser.parse_args()
    
    if args.plotonly:
        with open(f"{args.plotonly}", 'r') as f:
            data = json.load(f)
        angles = data["angles"]
        tors_R1 = data["tors_R1"]
        tors_R2 = data["tors_R2"]
        plot_histograms(angles, tors_R1, tors_R2, args.output)
        print(f"Histograms saved to {args.output}_histograms.png")
    else:
        universe = get_universe(args.topol, args.traj,args.atom_style, args.mda_format)
        angles, tors_R1, tors_R2 = calc_ring_shape(universe, args.output, args.a1, args.a2, args.a3, args.a4, args.r1, args.r2, args.max_dist,args.box)
        plot_histograms(angles, tors_R1, tors_R2, args.output)
        print(f"Ring shape data saved to {args.output}.json")
        print(f"Histograms saved to {args.output}_histograms.png")

if __name__ == "__main__":
    main()
    
    
            
            