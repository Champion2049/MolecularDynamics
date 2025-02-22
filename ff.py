import numpy as np
import mdtraj as md
import MDAnalysis as mda

# Load the PDB file
pdb_file = "4f51.pdb"
traj = md.load(pdb_file)
u = mda.Universe(pdb_file)

# Extract atom properties
atoms = u.atoms
atom_names = atoms.names

# Default partial charges (adjust if needed)
default_charges = {
    "C": 0.1, "O": -0.3, "N": -0.2, "H": 0.1, "Cl": -0.3,
}
charges = np.array([default_charges.get(atom[0], 0.0) for atom in atom_names])

# Extract topology and connectivity
topology = traj.topology
bond_indices = [[bond.atom1.index, bond.atom2.index] for bond in topology.bonds]

# Compute bond lengths
bond_lengths = md.compute_distances(traj, bond_indices)[0]

# Handle zero or very small bond lengths (avoid division by zero)
bond_lengths[bond_lengths < 0.1] = 0.1

# Manually find out the angles
angle_indices = []
for i, bond1 in enumerate(bond_indices):
    for bond2 in bond_indices[i + 1:]:
        shared_atom = set(bond1) & set(bond2)
        if shared_atom:
            a, b = bond1
            c, d = bond2
            common = shared_atom.pop()
            angle = [a, common, d] if common == b else [c, common, b]
            if angle not in angle_indices:
                angle_indices.append(angle)

# Compute bond angles and replace NaN values
bond_angles = np.degrees(md.compute_angles(traj, angle_indices)[0])
bond_angles = np.nan_to_num(bond_angles, nan=109.5)  # Replace NaN with default angle

# Compute dihedrals (torsions)
dihedral_indices = []
for i, angle1 in enumerate(angle_indices):
    for angle2 in angle_indices[i + 1:]:
        shared_atoms = set(angle1) & set(angle2)
        if len(shared_atoms) == 2:
            a, b, c = angle1
            d = angle2[2] if angle2[0] == b else angle2[0]
            dihedral = [a, b, c, d]
            if dihedral not in dihedral_indices:
                dihedral_indices.append(dihedral)

# Compute torsions
torsions = np.degrees(md.compute_dihedrals(traj, dihedral_indices)[0])

# Extract force field parameters (random values as placeholders)
force_constants = {f"{atoms[b[0]].name}-{atoms[b[1]].name}": np.random.uniform(100, 500) for b in bond_indices}

# Define energy functions
def bond_energy(l, l0, k_b):
    return 0.5 * k_b * (l - l0) ** 2

def angle_energy(theta, theta0, k_theta):
    return 0.5 * k_theta * (theta - theta0) ** 2

def torsion_energy(omega, V_n, n, gamma):
    return (V_n / 2) * (1 + np.cos(np.radians(n * omega) - gamma))

def vdw_energy(r, epsilon, Rij):
    return epsilon * ((Rij / r) * 12 - 2 * (Rij / r) * 6)

def electrostatic_energy(q_i, q_j, r, epsilon_0):
    r = max(r, 0.1)  # Minimum cutoff to prevent division by zero
    return (q_i * q_j) / (4 * np.pi * epsilon_0 * r)

# Compute energy components
E_bond = sum(bond_energy(l, 1.54, force_constants.get(f"{atoms[b[0]].name}-{atoms[b[1]].name}", 300)) for b, l in zip(bond_indices, bond_lengths))
E_angle = sum(angle_energy(theta, 109.5, 120) for theta in bond_angles)
E_torsion = sum(torsion_energy(omega, 12, 3, 0) for omega in torsions)
E_vdw = sum(vdw_energy(r, 0.3, 0.35) for r in bond_lengths)
E_electrostatic = sum(electrostatic_energy(charges[b[0]], charges[b[1]], l, 8.854e-12) for b, l in zip(bond_indices, bond_lengths))

# Total energy
E_total = E_bond + E_angle + E_torsion + E_vdw + E_electrostatic

# Print extracted parameters
print("\nExtracted Force Field Parameters")
print(" Bonds:", force_constants)
print(f"Angles: {len(angle_indices)} detected")
print(f"Torsions: {len(dihedral_indices)} detected")

# Print energy components
print("\nEnergy Components")
print(f"Bond Energy: {E_bond:.4f} kJ/mol")
print(f"Angle Energy: {E_angle:.4f} kJ/mol")
print(f"Torsion Energy: {E_torsion:.4f} kJ/mol")
print(f"Van der Waals Energy: {E_vdw:.4f} kJ/mol")
print(f"Electrostatic Energy: {E_electrostatic:.4f} kJ/mol")

print(f"\nTotal Potential Energy: {E_total:.4f} kJ/mol")