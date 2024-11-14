import math
import sys
import numpy as np

def mean_charge_dipole_interaction(charge_e, dipole_D, distance_A, temperature_K):
    """
    Calculate the mean interaction energy for a freely rotating dipole and a point charge.
    
    Parameters:
    - charge_e: Charge in units of elementary charge (e.g., e for an electron).
    - dipole_D: Dipole moment in Debye (D).
    - distance_A: Distance between charge and dipole in Angstroms (Å).
    - temperature_K: Temperature in Kelvin.
    
    Returns:
    - Mean interaction energy in kJ/mol.
    """
    # Constants
    e_charge = 1.602e-19  # Elementary charge in C
    D_to_Cm = 3.33564e-30  # Conversion factor from Debye to C·m
    epsilon_0 = 8.854e-12  # Permittivity of free space in C^2/(N·m^2)
    k_B = 1.3806e-23  # Boltzmann's constant in J/K
    NA = 6.022e23  # Avogadro's number
    
    # Convert inputs to SI units
    q = charge_e * e_charge  # Convert charge from e to C
    mu = dipole_D * D_to_Cm  # Convert dipole from Debye to C·m
    r = distance_A * 1e-10  # Convert distance from Å to m
    
    # Calculate mean interaction energy in J
    U_mean_J = - (q * mu**2) / (4 * math.pi * epsilon_0 * r**3 * k_B * temperature_K)
    U_mean_J_alt = - (q**2 * mu**2) / ((4 * math.pi * epsilon_0 * r**2)**2 * 3 * k_B * temperature_K)
    
    # Convert mean interaction energy to kJ/mol
    U_mean_kJ_per_mol = U_mean_J * 1e-3 * NA
    U_mean_alt_kJ_per_mol = U_mean_J_alt * 1e-3 * NA
    
    return(U_mean_kJ_per_mol,U_mean_alt_kJ_per_mol)



with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()
    
    
for line in lines:
    mol=line.split()[0]
    # Example usage
    charge_e = float(line.split()[1])  # charge in units of e (e.g., 1 for Cl-)
    dipole_D = float(line.split()[2])  # dipole moment in Debye
    distance_A = float(line.split()[3])  # distance in Å
    temperature_K = float(line.split()[4])  # temperature in K
    mean_interaction_energy = mean_charge_dipole_interaction(charge_e, dipole_D, distance_A, temperature_K)
    print(f"The mean interaction energy for {mol} is {mean_interaction_energy[0]:.2f} kJ/mol")
    print(f"The mean alt interaction energy for {mol} is {mean_interaction_energy[1]:.2f} kJ/mol")
    print("######")
    




