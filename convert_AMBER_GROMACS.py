import parmed as pmd
import sys

# convert AMBER topology[1] and crds[2] to GROMACS, CHARMM formats 
amber = pmd.load_file(sys.argv[1], sys.argv[2])
gmx=sys.argv[1].split('.')[0]


# Save a GROMACS topology and GRO files
amber.save('{}_GRO.top'.format(gmx))
amber.save('{}.gro'.format(gmx))
print(gmx)
