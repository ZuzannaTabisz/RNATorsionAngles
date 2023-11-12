# RNATorsionAngles

RNA Angle Extraction from PDB File

Overview
This Python code extracts RNA dihedral angles from a PDB file. The angles include alpha, beta, gamma, delta, epsilon, zeta, chi_purines, and chi_pyrimidines. The extracted angles are then presented in a 2D matrix.

Requirements

Python 3.x
Biopython library (pip install biopython)
NumPy library (pip install numpy)

Usage
bash
Copy code
python script.py <pdb_file_path>
Replace <pdb_file_path> with the path to the PDB file you want to analyze.

Output
File named angle_matrix.txt. The matrix as well as all the calculated dihedral angles are printed in console.
