import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np

# Load the PDB file
pdb = 'ini.pdb' # Path of PDB file
xtc = 'try.xtc' # Path of XTC file
u = mda.Universe(pdb, xtc)

# Select heavy atoms (excluding hydrogens)
protein_atoms = u.select_atoms("protein and not name H*")
dna_atoms = u.select_atoms("nucleic and not name H*")

with open('protein_dna_atom_min_dist.dat', 'w') as f:
    for ts in u.trajectory:
        # Calculate distance array between protein and DNA atoms
        dist_array = distances.distance_array(protein_atoms.positions, dna_atoms.positions)
        
        # Find the minimum distance for this frame
        min_distance = dist_array.min()
        
        # Write frame number and minimum distance to file
        f.write(f"{ts.frame} {min_distance:.4f}\n")
