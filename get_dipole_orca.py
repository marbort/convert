import sys 
with open(sys.argv[1], 'r') as f:
     lines=f.readlines()

for line in lines:
    if "Total Dipole" in line:
       dipole=line.split()[-3:]
    if "electric origin" in line:
        origin=line.split()[-3:]
    if "Magnitude (a.u.)" in line:
        magnitude=line.split()[-1]
   
    

dipole_angstrom=[float(i)*0.5291 for i in dipole]
origin_angstrom=[float(i)*0.5291 for i in origin]
end_angstrom=[float(origin_angstrom[i])+float(dipole_angstrom[i]) for i in range(3)]
print(f"Dipole moment (e x Å): {dipole}")
print(f"Dipole origin (Å): {origin}")
print(f"Dipole magnitude (a.u.): {magnitude}")

print(f"VMD command to draw dipole:")   
print(f"vmd_draw_arrow top {{{origin_angstrom[0]:.3f} {origin_angstrom[1]:.3f} {origin_angstrom[2]:.3f}}} \
{{{end_angstrom[0]:.3f} {end_angstrom[1]:.3f} {end_angstrom[2]:.3f}}}")
print(f"VMD command to draw dipole vector z-component:")
print(f"vmd_draw_arrow top {{{origin_angstrom[0]:.3f} {origin_angstrom[1]:.3f} {origin_angstrom[2]:.3f}}} \
{{{origin_angstrom[0]:.3f} {origin_angstrom[1]:.3f} {end_angstrom[2]:.3f}}}")
print(f"Cosine of angle with z-axis: {dipole_angstrom[2]/float(magnitude)/0.5291:.4f}")
