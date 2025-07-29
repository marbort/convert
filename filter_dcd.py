import  MDAnalysis as mda
import argparse

def get_trajectory(topol, traj):
    """
    Load the trajectory using MDAnalysis.
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    
    Returns:
    Universe: MDAnalysis Universe object containing the trajectory.
    """
    return mda.Universe(topol, traj)

def filter_dcd(topol, traj, output,begin=0, end=-1):
    """
    Filter the DCD trajectory to only include frames from begin to end.
    
    Parameters:
    topol (str): Path to the topology file.
    traj (str): Path to the trajectory file.
    output (str): Path to save the filtered trajectory.
    """
    u = get_trajectory(topol, traj)
    
    with mda.Writer(output, u.atoms.n_atoms) as W:
        print(f"Total frames in trajectory: {len(u.trajectory)}")
        print(f"Filtering frames from {begin} to {end}...")
        for ts in u.trajectory[begin:end]:
            if u.atoms[0].position[2] > 0:  # Check z-coordinate of the first atom
                W.write(u)

def main():
    parser = argparse.ArgumentParser(description="Filter DCD trajectory from begin to end frame.")
    parser.add_argument("topol", type=str, help="Path to the topology file.")
    parser.add_argument("traj", type=str, help="Path to the trajectory file.")
    parser.add_argument("output", type=str, help="Path to save the filtered trajectory.")
    parser.add_argument("--begin", type=int, default=0, help="Start frame index (inclusive).")
    parser.add_argument("--end", type=int, default=-1, help="End frame index (exclusive). Use -1 for all frames.")

    args = parser.parse_args()
    
    filter_dcd(args.topol, args.traj, args.output, args.begin, args.end)

if __name__ == "__main__":
    main()

