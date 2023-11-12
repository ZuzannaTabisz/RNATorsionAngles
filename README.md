# RNATorsionAngles

RNA Angle Extraction from PDB File

**Overview**

This Python code extracts RNA dihedral angles from a PDB file. The angles include alpha, beta, gamma, delta, epsilon, zeta, chi_purines, and chi_pyrimidines. The extracted angles are then presented in a 2D matrix.

**Requirements**

Python 3.x
Biopython library (pip install biopython)
NumPy library (pip install numpy)

**Usage**

 ```bash
    python z2rna.py <pdb_file_path>
```
Replace <pdb_file_path> with the path to the PDB file you want to analyze.


**Output**

File named angle_matrix.txt, with the first row having angle symbols. The matrix as well as all the calculated dihedral angles are printed in console.

Example:
α	β	γ	δ	ε	ζ	χ
-	-128.1	67.8	82.9	-155.6	-68.6	-167.8
-67.4	-178.4	53.8	83.4	-145.1	-76.8	-163.8
-74.5	169.7	59.5	80.7	-148.3	-80.0	-161.9
-64.4	162.2	60.7	82.2	-157.4	-68.7	-168.7
-74.7	-176.5	53.4	84.9	-137.5	-81.7	-162.9
-48.8	157.6	55.3	81.3	-151.0	-77.0	-160.0
-59.5	-178.7	62.5	137.3	-105.9	-52.0	-133.1
-83.8	-145.6	55.4	78.6	-142.8	-118.6	-161.5
-69.7	-141.7	52.3	147.8	-106.2	-77.3	-70.5
177.8	147.2	60.1	89.3	-126.2	-88.7	169.6
-56.1	167.9	48.2	87.2	-150.5	-69.9	-160.9
-67.8	172.9	51.8	80.7	-158.5	-65.2	-158.3
166.6	-169.9	178.6	82.5	-153.1	-97.4	-168.3
83.4	-158.3	-114.6	92.0	-125.5	-57.3	-170.7
-55.1	162.5	51.9	79.8	-136.3	-143.9	-164.5
-6.1	91.2	76.8	96.8	-61.8	-131.2	-85.8
27.8	107.7	174.1	94.8	178.0	76.2	-142.5
45.4	-159.4	59.0	150.6	-95.2	-179.1	-99.5
-71.4	-178.9	53.8	153.8	-91.6	-83.7	-80.3
-81.3	-150.7	47.8	89.9	-122.3	-54.1	177.8
-75.6	148.6	-176.6	78.2	-168.9	-75.6	-160.2
158.8	153.5	179.3	82.0	-145.0	-80.4	-175.5
-53.3	174.8	52.5	82.3	-155.3	-66.4	-158.0
-68.8	178.2	46.8	83.6	-144.3	-72.8	-160.7
-65.1	168.9	53.9	83.3	-145.1	-68.4	-160.3
-53.8	170.8	47.7	86.0	-136.3	-76.9	-163.4
-53.0	166.9	43.6	83.4	-148.5	-73.4	-168.2
-72.4	178.3	49.3	80.1	-152.1	-67.0	-160.6
-66.6	174.0	55.6	81.4	-155.5	-78.3	-165.9
-54.0	165.9	56.9	83.6	-144.7	-62.3	-171.7
-69.9	177.8	52.3	83.7	-137.0	-75.5	-156.7
-52.7	161.4	49.3	80.1	-145.9	-71.2	-149.9
-67.7	-177.0	47.0	82.1	-148.0	-53.7	-148.2
171.1	148.1	52.5	83.4	-132.5	-71.8	-171.2
-47.7	163.7	40.2	80.9	-143.7	-59.5	-154.4
-52.4	165.7	51.3	72.2	-160.4	-85.2	-158.4
-57.5	163.0	47.8	81.1	-148.1	-67.0	-168.8
-61.8	-180.0	46.9	82.5	-136.8	-76.4	-169.4
-47.7	160.4	53.3	79.3	-140.1	-68.6	-166.9
-67.4	172.0	56.2	83.2	-154.2	-74.9	-162.6
-68.2	-179.4	52.4	78.9	-137.3	-84.7	-169.0
-47.9	158.7	55.6	79.8	-160.3	-70.3	-169.0
-67.0	-178.3	55.6	81.6	-154.9	-76.4	-160.2
-59.7	162.1	60.0	85.3	-142.8	-57.2	-159.4
-71.9	-176.9	51.0	87.6	-135.1	-78.7	-149.3
-56.8	-146.5	48.4	141.6	-102.7	-137.9	-65.8
62.4	-164.0	44.4	146.1	-93.7	-78.0	-112.0
-73.5	-174.3	161.5	145.6	-143.5	75.6	-140.1
50.7	168.5	42.2	84.3	-145.0	-82.1	-173.6
-51.7	177.2	42.1	80.4	-150.6	-67.8	-165.3
-63.9	176.8	52.8	79.4	-150.4	-71.3	-156.6
-64.7	173.6	48.5	80.3	-156.5	-69.4	-164.0
-56.9	171.5	56.2	83.9	-159.4	-64.9	-169.2
-79.7	-172.8	57.7	77.6	-128.6	-70.7	-161.5
-49.7	168.8	44.1	76.6	-140.8	-69.9	-152.4
166.4	171.8	53.3	83.4	-132.7	-70.6	-161.5
-65.7	167.1	57.5	81.7	-145.2	-67.6	-159.3
-60.8	-146.1	71.8	156.7	-78.3	-169.3	-86.3
72.6	-158.8	63.7	84.6	-148.8	-53.7	-165.6
-72.2	179.5	66.0	148.3	-97.1	-66.4	-117.8
-84.3	179.8	38.2	83.0	-152.3	-74.5	-166.7
-60.1	179.6	46.9	80.5	-145.6	-74.1	-158.7
-62.0	167.3	50.9	80.7	-152.3	-70.7	-152.6
-66.9	180.0	44.1	75.8	-147.5	-76.5	-161.8
-44.0	164.2	49.9	79.8	-152.0	-73.3	-172.8
-57.9	178.5	52.0	81.7	-151.0	-73.5	-164.9
-62.0	164.1	54.2	83.2	-152.2	-78.3	-162.8
-59.8	175.3	47.3	82.2	-152.9	-65.4	-160.1
-63.8	168.1	55.1	79.1	-155.4	-85.6	-161.4
-61.7	164.6	53.1	79.0	-158.5	-64.5	-152.0
-78.4	173.6	60.3	80.3	-149.6	-68.4	-162.8
-73.2	176.2	62.1	83.0	-152.3	-67.9	-161.6
-63.3	177.7	50.4	81.6	-148.2	-66.1	-167.4
-66.9	-174.9	50.7	85.9	-145.0	-58.8	-153.1
-52.3	175.7	42.3	85.6	-131.9	163.9	-151.7
-71.0	130.2	164.6	160.9	-	-	138.5
