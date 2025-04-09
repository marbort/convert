import MDAnalysis as mda





def rescale_charges(infile, factor, outfile):
    u = mda.Universe(infile)
    charges = u.atoms.charges
    "rescale charges by factor"
    charges=charges*factor
    u.atoms.charges=charges
    with mda.coordinates.MOL2.MOL2Writer(outfile) as writer:
        writer.write(u)
    return u

for i in [0.95,0.90,0.85,0.80]:
    rescale_charges("AcPh.mol2", i, f"AcPh_rescaled_{int(i*100)}.mol2")
    

    