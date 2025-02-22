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
### 1. Analyzing the Trajectory
Run the Python script to analyze the trajectory:

```bash
python3 trajectory.py
```

### 2. Visualization
Generate plots for analysis:

```bash
python3 viewxvg.py
python3 rmsd.py
```

### 3. Simulating the Trajectory
Using PyMOL or VMD, load the md_0_1.gro file and on top of that file load md_0_1_noPBC.xtc file and click on play to run the simulation of the molecule.

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

