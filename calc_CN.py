import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
import sys

# Load the RDF data from a GROMACS .xvg file
def load_rdf(xvg_file):
    data = np.loadtxt(xvg_file, comments=['#', '@'])  # Ignore header lines
    r = data[:, 0]  # Radial distances (r)
    g_r = data[:, 1]  # RDF values (g(r))
    return r, g_r

def get_rho(num_particles, volume):
    volume = volume **3   # Convert volume from A^3 to nm^3 
    return num_particles / volume

# Function to calculate the coordination number
def calculate_coordination_number(r, g_r, rho, r_min, r_max):
    # We integrate g(r) * r^2 over the specified range [r_min, r_max]
    # First, find the indices corresponding to r_min and r_max
    peak=g_r.max()
    peak_idx = np.argmax(g_r)
    peak_distance = r[peak_idx]
    r_min_idx = np.searchsorted(r, r_min)
    r_max_idx = np.searchsorted(r, r_max)
    
    # Slice the RDF data within the desired range
    r_sliced = r[r_min_idx:r_max_idx]
    g_r_sliced = g_r[r_min_idx:r_max_idx]
    
    # Calculate the integral using Simpson's rule for numerical integration
    integral = simps(g_r_sliced * r_sliced**2, r_sliced)
    
    # The coordination number CN
    CN = 4 * np.pi * rho * integral
    return CN,peak_distance,peak

# Example usage:
xvg_file = sys.argv[1]  # Your GROMACS RDF file
r, g_r = load_rdf(xvg_file)

# Number density (rho) in units of particles per nm^3, for example
#rho = 1.0  # Example: 1 particle per nm^3
rho=get_rho(float(sys.argv[2]),float(sys.argv[3]))

# Define the range for integration (in nm)
r_min = 0.0  # Example minimum radius (nm)
r_max = sys.argv[4]  # Example maximum radius (nm)

# Calculate the coordination number

CN,peak_distance,peak = calculate_coordination_number(r, g_r, rho, r_min, r_max)

print(f"Peak Distance: {peak_distance}")
print(f"Peak Value: {peak}")
print(f"Coordination Number: {CN}")

