import numpy as np

#Calculate the electrostatic potential energy of a atomic size charged sphere in a solvent
def charged_sphere_potential(radius, charge, permittivity):
    """
    Calculate the electrostatic potential energy of a charged sphere in a solvent.
    
    Parameters:
    radius (float): Radius of the sphere in Angstrom.
    charge (float): Charge of the sphere in electrons.
    permittivity (float): Relative permittivity of the solvent.
    
    Returns:
    float: Electrostatic potential energy in joules.
    """
    k = 8.9875517873681764e9  # Coulomb's constant in N m²/C²
    angstrom_to_meter = 1e-10  # Conversion factor from Angstrom to meters
    radius_meters = radius * angstrom_to_meter
    charge_coulombs = charge * 1.602176634e-19  # Conversion factor from electrons to coulombs
    e_0 = 8.854187817e-12  # Permittivity of free space in F/m
    # Calculate the potential energy using the formula U = k * q^2 / (r * ε)
    # where U is the potential energy, k is Coulomb's constant, q is charge, r is radius, and ε is permittivity.
    if radius_meters == 0:
        raise ValueError("Radius must be greater than zero to avoid division by zero.")
    if permittivity == 0:
        raise ValueError("Permittivity must be greater than zero to avoid division by zero.")
    
    return (k * charge_coulombs**2) / (radius_meters * permittivity*2)

def main():
    # Example parameters
    radius = 13  # Radius in Angstroms
    charge = 1  # Charge in electrons
    permittivity = 7.58  # Permittivity of free space in F/m

    potential_energy = charged_sphere_potential(radius, charge, permittivity)
    print(f"Electrostatic potential energy of a charged sphere with radius {radius} Å and charge {charge} e in a solvent with relative permittivity {permittivity}:")
    print(f"Electrostatic potential energy: {potential_energy:.2e} J")
    #convert to kcal/mol for convenience
    potential_energy_kcal = potential_energy / 4.184e3 * 6.022e23  # Convert J to kcal/mol
    print(f"Electrostatic potential energy: {potential_energy_kcal:.2f} kcal/mol")

    
if __name__ == "__main__":
    main()