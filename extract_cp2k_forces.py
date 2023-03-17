import os

dir_path = os.getcwd()

for filename in os.listdir(dir_path):
    if filename.endswith(".out"):
        with open(os.path.join(dir_path, filename), 'r') as f:
            lines = f.readlines()

        # Find the start of the atomic forces section
        start = 0
        while "ATOMIC FORCES" not in lines[start]:
            start += 1

        # Find the end of the atomic forces section
        end = start
        while "SUM OF ATOMIC FORCES" not in lines[end]:
            end += 1

        # Extract the atomic forces
        forces = []
        for i in range(start+3, end):
            try:
                atom_index = int(lines[i].split()[0])
                kind = int(lines[i].split()[1])
                element = lines[i].split()[2]
                force = [float(x) for x in lines[i].split()[3:6]]
                total_norm = sum([x**2 for x in force])**0.5
                forces.append((atom_index, kind, element, force, total_norm))
            except ValueError:
                pass

        # Find the atom with the highest norm
        max_norm_atom = max(forces, key=lambda atom: atom[4])
        max_norm_line = f"{max_norm_atom[0]:>6d}   {max_norm_atom[1]:>4d}   {max_norm_atom[2]:>2s}{''.join([f'{x:16.8f}' for x in max_norm_atom[3]])}{max_norm_atom[4]:16.8f}{max_norm_atom[4]*51.421:16.8f}\n"

        # Write the atomic forces to a new file
        output_filename = filename.split(".")[0] + "_forces.log"
        with open(output_filename, "w") as f:
            f.write("  ATOM     KIND     ELEMENT       FX             FY              FZ          |FORCE( a.u)|   |FORCE (eV/Å)|\n")
            f.write("-----------------------------------------------------------------------------------------------------------\n")
            for atom in forces:
                atom_index, kind, element, force, total_norm = atom
                f.write(f"{atom_index:>6d}   {kind:>6d}   {element:>6s}{''.join([f'{x:16.8f}' for x in force])}{total_norm:16.8f}{total_norm*51.421:16.8f}\n")
            f.write(f"Atom with highest norm:")
            f.write(f"{max_norm_line}")
