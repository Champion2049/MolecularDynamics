import MDAnalysis as mda

# Load trajectory and topology
u = mda.Universe("md_0_1.tpr", "md_0_1.xtc")

# Write out DCD file
with mda.Writer("md_0_1.dcd", u.atoms.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u.atoms)
