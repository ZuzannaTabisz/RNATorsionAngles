import io
import sys
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import calc_dihedral
from math import degrees

def extract_angles(structure):
    angles = {
        "alpha": ("O3'", "P", "O5'", "C5'"),
        "beta": ("P", "O5'", "C5'", "C4'"),
        "gamma": ("O5'", "C5'", "C4'", "C3'"),
        "delta": ("C5'", "C4'", "C3'", "O3'"),
        "epsilon": ("C4'", "C3'", "O3'", "P"),
        "zeta": ("C3'", "O3'", "P", "O5'"),
        "chi_purines": ("O4'", "C1'", "N9", "C4"),
        "chi_pyrimidines": ("O4'", "C1'", "N1", "C2"),
    }

    all_residue_angles = {}

    for model in structure:
        for chain in model:
            for i, residue in enumerate(chain):
                angles_values = {}
                #From the dictionary above
                for angle_name, atoms in angles.items():
                    try:
                        #All the atoms that should be conisered form previous or next residue
                        if angle_name == "alpha" and i == 0:
                            raise e
                        elif angle_name == "alpha" and i > 0:
                            previous_residue = chain.get_list()[i - 1]
                            atoms_objects = [previous_residue[atoms[0]]] + [residue[atoms[i]] for i in range(1, 4)]
                            # print(atoms_objects,"alpha")
                        elif angle_name == "epsilon":
                            next_residue = chain.get_list()[i + 1]
                            atoms_objects =  [residue[atoms[i]] for i in range(0, 3)] + [next_residue[atoms[3]]]
                            # print(atoms_objects,"epsilon")
                        elif angle_name == "zeta":
                            next_residue = chain.get_list()[i + 1]
                            atoms_objects = [residue[atoms[i]] for i in range(0, 2)] + [next_residue[atoms[i]] for i in range(2,4)]
                            # print(atoms_objects,"zeta")
                        else:
                            atoms_objects = [residue[atoms[i]] for i in range(4)]
                        vectors = [atom.get_vector() for atom in atoms_objects]
                        # print(vectors)
                        #Calc_dihedral returns radians - converting it to degrees
                        angle = calc_dihedral(vectors[0], vectors[1], vectors[2], vectors[3])
                        angle_degrees = round(degrees(angle), 1)

                        angles_values[angle_name] = angle_degrees
                    except Exception as e:
                        angles_values[angle_name] = "-"
                # Saving only the residues that have at least one angle
                if any(angle != "-" for angle in angles_values.values()):
                    all_residue_angles[residue.id[1]] = angles_values

    return all_residue_angles

def save_angles_to_2Dmatrix(all_residue_angles, structure):

    total_residues = sum(1 for _ in structure.get_residues())

    #matrix is only as big as all_residue angles
    matrix = np.empty((len(all_residue_angles), 7), dtype=object)

    #Still we go through all the residues, because we need to check whether there are purines or pyrmidines
    for i, residue in enumerate(structure.get_residues()):

        index = i
        #Here we break before the residues that we dont need
        if index == len(all_residue_angles):
            break
        angles_values = all_residue_angles.get(residue.id[1], {})

        matrix[index, 0] = angles_values.get("alpha", "-")
        matrix[index, 1] = angles_values.get("beta", "-")
        matrix[index, 2] = angles_values.get("gamma", "-")
        matrix[index, 3] = angles_values.get("delta", "-")
        matrix[index, 4] = angles_values.get("epsilon", "-")
        matrix[index, 5] = angles_values.get("zeta", "-")

        #Adjusting to pyrimidines and purines as shown at the website : http://rnafrabase.ibch.poznan.pl/?act=search
        if any(letter in residue.get_resname() for letter in ["A", "G"]):
            matrix[index, 6] = angles_values.get("chi_purines", "-")
        elif any(letter in residue.get_resname() for letter in ["C", "U"]):
            matrix[index, 6] = angles_values.get("chi_pyrimidines", "-")
        else:
            matrix[index, 6] = "-"

    return matrix


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file_path>")
        sys.exit(1)

    pdb_file_path = sys.argv[1]

    #We parse only ATOM and HETATM residues as in the results form : http://rnafrabase.ibch.poznan.pl/?act=search
    atom_lines = [line.strip() for line in open(pdb_file_path) if line.startswith("ATOM") or line.startswith("HETATM")]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("residue", io.StringIO('\n'.join(atom_lines)))

    residue_angles = extract_angles(structure)


    print("res no.\tα\tβ\tγ\tδ\tε\tζ\tχ_purines\tχ_pyrimidines")


    for residue_num, angles_values in residue_angles.items():
        print(f"{residue_num}\t"
            f"{angles_values.get('alpha', '-')}\t"
            f"{angles_values.get('beta', '-')}\t"
            f"{angles_values.get('gamma', '-')}\t"
            f"{angles_values.get('delta', '-')}\t"
            f"{angles_values.get('epsilon', '-')}\t"
            f"{angles_values.get('zeta', '-')}\t"
            f"{angles_values.get('chi_purines', '-')}\t"
            f"{angles_values.get('chi_pyrimidines', '-')}")

    #Saving the matrix
    matrix = save_angles_to_2Dmatrix(residue_angles, structure)


    print(matrix)

    output_file_path = "angle_matrix.txt"
    np.savetxt(output_file_path, matrix, fmt='%s', delimiter='\t', header="α\tβ\tγ\tδ\tε\tζ\tχ", comments='', encoding='utf-8')


if __name__ == "__main__":
    main()
