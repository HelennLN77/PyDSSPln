import pymol
from pymol import cmd
import os

def add_hydrogens_to_pdb(input_pdbs: list, output_directory: str):
    """
    Add hydrogen atoms to multiple PDB files using PyMOL and save the results to a specified output directory.

    Parameters:
    input_pdbs (list): List of paths to input PDB files.
    output_directory (str): Directory to save the output PDB files with added hydrogens.
    """
    pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet mode

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for input_pdb in input_pdbs:
        try:
            if not os.path.exists(input_pdb):
                print(f"Error: File {input_pdb} does not exist. Skipping.")
                continue

            # Load the PDB file into PyMOL
            cmd.load(input_pdb)

            # Add hydrogens
            cmd.h_add()

            # Extract the filename without the path and extension
            base_filename = os.path.basename(input_pdb).replace(".pdb", "_hydro.pdb")
            output_pdb = os.path.join(output_directory, base_filename)

            # Save the new PDB file with added hydrogens
            cmd.save(output_pdb)

            # Clear the PyMOL session for the next file
            cmd.reinitialize()

            print(f"Hydrogens added to {input_pdb}. Output saved to {output_pdb}")

        except Exception as e:
            print(f"An error occurred while processing {input_pdb}: {e}")

    pymol.cmd.quit()


# Example usage
input_files = ["Documents/M2BI/Programmation_avancee/projet_court_DSSP/1qk1.pdb", "Documents/M2BI/Programmation_avancee/projet_court_DSSP/1shr.pdb","Documents/M2BI/Programmation_avancee/projet_court_DSSP/7kqy.pdb"]
output_dir = "output_pdbs"
add_hydrogens_to_pdb(input_files, output_dir)
