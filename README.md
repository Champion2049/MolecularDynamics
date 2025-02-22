# Molecular Dynamics Analysis of 4F51

![Molecular Dynamics](https://upload.wikimedia.org/wikipedia/commons/3/3d/Molecular_dynamics_simulation.gif)

## Overview
This project utilizes **MDAnalysis** and **GROMACS** to perform molecular dynamics simulations and analyze the trajectory of the molecule with the **PDB ID: 4F51**. The study focuses on evaluating the molecule's properties based on its trajectory data.

## Features
- **Molecular Dynamics Simulation:** Uses GROMACS to generate trajectory data.
- **Trajectory Analysis:** Leverages MDAnalysis to extract and analyze key molecular properties.
- **Visualization:** Provides graphical representations of simulation results.

## Installation
Ensure you have the following dependencies installed:

```bash
pip install MDAnalysis matplotlib numpy pandas
```

For GROMACS, install it via:

```bash
sudo apt install gromacs  # Linux
# or
brew install gromacs  # macOS
```

## Usage
### 1. Setting Up the System
Prepare the molecular system and configure GROMACS for simulation:

```bash
gmx pdb2gmx -f 4f51.pdb -o processed.gro -water spce
gmx grompp -f md.mdp -c processed.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md
```

### 2. Analyzing the Trajectory
Run the Python script to analyze the trajectory:

```bash
python analyze_trajectory.py
```

### 3. Visualization
Generate plots for analysis:

```python
import MDAnalysis as mda
import matplotlib.pyplot as plt

u = mda.Universe("md.tpr", "md.xtc")

# Example: Calculate and plot RMSD
rmsd = []
for ts in u.trajectory:
    rmsd.append(ts.time)

plt.plot(rmsd)
plt.xlabel("Time (ps)")
plt.ylabel("RMSD")
plt.show()
```

## Results
- **Trajectory of the Molecule:** Calculates the Molecules trajectory to get simulations on VMD or PyMOL
- **Root Mean Square Deviation (RMSD):** Measures the structural stability over time.
- **Root Mean Square Fluctuation (RMSF):** Calculation of Individual Residue Flexibility
- **Radius of Gyration (Rg):** Evaluates compactness of the molecule.
- **Hydrogen Bond Analysis:** Examines molecular interactions.

## Contributing
Feel free to submit issues or pull requests for enhancements.

## References
- [MDAnalysis Documentation](https://www.mdanalysis.org/)
- [GROMACS Manual](http://www.gromacs.org/)
- [PDB ID 4F51](https://www.rcsb.org/structure/4F51)

---
For more details, visit the [GitHub Repository](https://github.com/Champion2049/MolecularDynamics).

