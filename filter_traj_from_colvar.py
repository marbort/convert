import MDAnalysis as mda 
import numpy as np



def get_frames_to_extract(colvar,timestep):
    data=np.loadtxt(colvar,skiprows=1,unpack=True)
    frames_to_extract = [x / timestep for x in data[0]]
    return frames_to_extract

def get_filter_traj_from_colvar(frames_to_extract, topol, traj, output_traj):
    """
    Extract frames from a trajectory based on a colvar file.
    
    Parameters:
    frames_to_extract (list): List of frame indices to extract.
    topol (str): Path to the topologyfile
    traj (str): Path to the trajectory file.
    output_traj (str): Path to save the filtered trajectory.
    
    Returns:
    None
    """
    u = mda.Universe(topol,traj,atom_style='id type x y z',format='LAMMPSDUMP')
    with mda.Writer(output_traj, u.atoms.n_atoms) as W:
        for ts in frames_to_extract:
            ts = int(ts)  # Ensure ts is an integer index
            if ts < len(u.trajectory):  # Check if the frame index is within bounds
                u.trajectory[ts]
                # Write the current frame to the output trajectory
                # Note: mda.trajectory.TrajectoryWriter expects the frame index to be an integer
                # and the atoms to be written in the format specified by the topology.
                # Here, we write the atoms of the current frame to the output trajectory.
                # If the frame index is valid, write the atoms to the output trajectory.
                W.write(u.atoms)

def main(): 
    import argparse
    parser = argparse.ArgumentParser(description='Filter trajectory based on colvar file')
    parser.add_argument('--colvar', type=str, help='Colvar file to filter', required=True)
    parser.add_argument('--topol', type=str, help='Topology file', required=True)
    parser.add_argument('--traj', type=str, help='Trajectory file', required=True)
    parser.add_argument('--output_traj', type=str, help='Output trajectory file', required=True)
    parser.add_argument('--timestep', type=float, help='Timestep in ps', default=1.0)

    args = parser.parse_args()

    frames_to_extract = get_frames_to_extract(args.colvar, args.timestep)
    get_filter_traj_from_colvar(frames_to_extract, args.topol, args.traj, args.output_traj)
    
if __name__ == "__main__":
    main()
