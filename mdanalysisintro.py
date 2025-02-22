import mdtraj as md
import numpy as np
import os

pdf_file = '7qdi.pdb'
if not os.path.exists(pdf_file):
    print('File not found: {}'.format(pdf_file))
    exit()

atom_count = 100
traj = md.load(pdf_file)
for i in range (atom_count):
    atom1_index = i
    for j in range (i+1, atom_count):
        atom2_index = j
        distances = md.compute_distances(traj, [[atom1_index, atom2_index]])
        print(f"Distance between atom {(atom1_index)+1} and atom {(atom2_index)+1} is {distances[0][0]:.3f} Angstroms")
        md.compute_angles(traj, [[atom1_index, atom2_index, atom2_index+1]])
        print(f"Angle between atom {(atom1_index)+1}, atom {(atom2_index)+1}, and atom {(atom2_index)+1} is {distances[0][0]:.3f} degrees")
bond_distance_threshold = 2.0
bond_counts = np.zeros(atom_count)
for i in range(atom_count):
    for j in range(i + 1, atom_count):
        distance = np.linalg.norm(traj.xyz[0][i] - traj.xyz[0][j])
        if distance < bond_distance_threshold:
            bond_counts[i] += 1
            bond_counts[j] += 1
for i in range(atom_count):
    print(f"Atom {i+1} ({traj.atom(i).name}) has {bond_counts[i]} bonds.")