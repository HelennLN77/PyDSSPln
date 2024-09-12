# PyDSSPln
# Python implementation of DSSP algorithm for secondary structure in proteins
# Date: 2024-09-12
# Author: Bounsay Hélène

# Importing libraries
import os
import numpy as np
import multiprocessing as mp
import logging
from Bio.PDB import PDBParser, DSSP
from sklearn.metrics import confusion_matrix, accuracy_score

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class HBond:
    """
    Represents a hydrogen bond between two residues.
    
    Attributes:
        co_i (dict): Coordinates of the donor residue atoms.
        nh_j (dict): Coordinates of the acceptor residue atoms.
        cutoff (float): Energy cutoff to determine if a hydrogen bond exists which is dertermine here at -0.5.
        energy (float): Calculated energy of the hydrogen bond.
    """

    def __init__(self, co_i, nh_j, cutoff=-0.5):
        """
        Initializes an HBond instance.
        
        Args:
            co_i (dict): Coordinates of the donor residue atoms.
            nh_j (dict): Coordinates of the acceptor residue atoms.
            cutoff (float): Energy cutoff for hydrogen bond detection.
        """
        self.co_i = co_i
        self.nh_j = nh_j
        self.cutoff = cutoff
        self.energy = self.calculate_energy()

    def calculate_energy(self):
        """
        Calculates the energy of the hydrogen bond.
        
        Returns:
            float: The energy of the hydrogen bond.
            q1, q2, f are constants used in the calculation of the energy of the hydrogen bond given by the article of Kabsch and Sander.
        """
        q1 = 0.42  # Partial charge on C,O
        q2 = 0.20  # Partial charge on N,H
        f = 332.0  # Dimensional factor 

        try:
            r_on = np.linalg.norm(self.co_i['O'] - self.nh_j['N'])
            r_ch = np.linalg.norm(self.co_i['C'] - self.nh_j['H'])
            r_oh = np.linalg.norm(self.co_i['O'] - self.nh_j['H'])
            r_cn = np.linalg.norm(self.co_i['C'] - self.nh_j['N'])
        except KeyError as e:
            logging.warning(f"Missing atom coordinate in HBond calculation: {e}")
            return float('inf')

        E = q1 * q2 * (1/r_on + 1/r_ch - 1/r_oh - 1/r_cn) * f
        return E
    
    def is_HBond(self):
        """
        Determines if the calculated energy indicates a hydrogen bond.
        
        Returns:
            bool: True if the energy is below the cutoff, indicating a hydrogen bond.
        """
        return self.energy < self.cutoff 

class SecondaryStructure:
    """
    Analyzes the secondary structure of a protein based on its residues.
    
    Attributes:
        residues (list): List of residue dictionaries with atom coordinates.
        hbonds (list): List of detected hydrogen bonds.
    """

    def __init__(self, residues):
        """
        Initializes a SecondaryStructure instance.
        
        Args:
            residues (list): List of residue dictionaries with atom coordinates.
        """
        self.residues = residues
        self.hbonds = self.calculate_hbonds()

    def calculate_hbonds(self):
        """
        Calculates all possible hydrogen bonds between residues.
        
        Returns:
            list: A list of tuples representing hydrogen bonds.
        """
        h_bonds = []
        for i, res_i in enumerate(self.residues):
            for j, res_j in enumerate(self.residues):
                if abs(i - j) >= 3:
                    h_bond = HBond(res_i['CO'], res_j['NH'])
                    if h_bond.is_HBond():
                        h_bonds.append((i, j, h_bond))
        logging.debug(f"Detected {len(h_bonds)} hydrogen bonds")
        return h_bonds

    def detect_turns(self, n):
        """
        Detects turns of a given size in the protein structure.
        
        Args:
            n (int): The size of the turn.
        
        Returns:
            list: A list of detected turns of the specified size.
        """
        turns = []
        for i, j, hbond in self.hbonds:
            if j == i + n:
                turns.append((i, j, n))
        logging.debug(f"Detected {len(turns)} turns of size {n}")
        return turns
    
    def detect_alpha_helix(self):
        """
        Detects alpha helices of sizes 3, 4, and 5 based on turns.
        
        Returns:
            dict: A dictionary of detected alpha helices categorized by size.
        """
        alpha_helix = {'3-helix': [], '4-helix': [], '5-helix': []}
        
        for n in range(3, 6):
            turns = self.detect_turns(n)
            for i, j, _ in turns:
                if (i + 1, i + n - 1, n) in turns:
                    alpha_helix[f'{n}-helix'].append((i, j))
        
        logging.debug(f"Detected helices: {alpha_helix}")
        return alpha_helix
    
    def detect_bridges(self):
        """
        Detects parallel and antiparallel bridges in the protein structure.
        
        Returns:
            dict: A dictionary of detected bridges categorized by type.
        """
        bridges = {'parallel': [], 'antiparallel': []}
    
        for i, j, hbond_ij in self.hbonds:
            if (i - 1, j) in [(hbond[0], hbond[1]) for hbond in self.hbonds] and \
               (j, i + 1) in [(hbond[0], hbond[1]) for hbond in self.hbonds]:
                bridges['parallel'].append((i, j))
            
            if (j - 1, i) in [(hbond[0], hbond[1]) for hbond in self.hbonds] and \
               (i, j + 1) in [(hbond[0], hbond[1]) for hbond in self.hbonds]:
                bridges['antiparallel'].append((i, j))
        logging.debug(f"Detected bridges: {bridges}")
        return bridges

    def detect_beta_sheet(self):
        """
        Detects beta sheets in the protein structure based on detected bridges.
        
        Returns:
            list: A list of detected beta sheets and ladders.
        """
        beta_sheet = []
        bridges = self.detect_bridges()

        for bridge_type, bridge_list in bridges.items():
            ladder = []
            for i, (res_i, res_j) in enumerate(bridge_list):
                if i == 0 or res_i != bridge_list[i - 1][0] + 1:
                    if ladder:
                        beta_sheet.append(ladder)
                    ladder = [(res_i, res_j, bridge_type)]
                else:
                    ladder.append((res_i, res_j, bridge_type))

            if ladder:
                beta_sheet.append(ladder)
        logging.debug(f"Detected beta sheets: {beta_sheet}")
        return beta_sheet

class DSSP:
    """
    Analyzes protein secondary structure using the DSSP algorithm.
    
    Attributes:
        pdb_file (str): Path to the PDB file.
        residues (list): List of residue dictionaries with atom coordinates.
        secondary_structure (SecondaryStructure): An instance of SecondaryStructure.
    """

    def __init__(self, pdb_file):
        """
        Initializes a DSSP instance.
        
        Args:
            pdb_file (str): Path to the PDB file.
        """
        self.residues = self.parse_pdb(pdb_file)
        self.secondary_structure = SecondaryStructure(self.residues)

    def parse_pdb(self, pdb_file):
        """
        Parses the PDB file to extract residue coordinates.
        
        Args:
            pdb_file (str): Path to the PDB file.
        
        Returns:
            list: A list of residue dictionaries with atom coordinates.
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_file)
        residues = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() in ['CYS', 'SER', 'THR', 'ARG', 'GLY']:
                        try:
                            co = residue['CA'].get_coord()
                            nh = residue['N'].get_coord()
                            co_o = residue['O'].get_coord() if 'O' in residue else np.zeros(3)
                            nh_h = residue['H'].get_coord() if 'H' in residue else np.zeros(3)
                            residues.append({'CO': {'C': co, 'O': co_o},
                                            'NH': {'N': nh, 'H': nh_h}})
                        except KeyError as e:
                            logging.warning(f"Missing atom in residue {residue.get_id()}: {e}")
                            continue
        return residues
    
    def get_secondary_structure(self):
        """
        Prints the detected secondary structure elements including helices and beta sheets.
        """
        helix = self.secondary_structure.detect_alpha_helix()
        sheets = self.secondary_structure.detect_beta_sheet()

        print("Detected helices: ")
        for helix_type, helix in helix.items():
            print(f"{helix_type}: {helix}")

        print("Detected beta sheets and ladders: ")
        for sheet in sheets:
            print(sheet)

def process_file(pdb_file, output_file):
    """
    Processes a single PDB file to analyze its secondary structure and writes results to an output file.
    
    Args:
        pdb_file (str): Path to the PDB file.
        output_file (str): Path to the output file where results will be written.
    """
    logging.info(f"Processing file: {pdb_file}")
    try:
        dssp = DSSP(pdb_file)
        result = []

        helices = dssp.secondary_structure.detect_alpha_helix()
        beta_sheets = dssp.secondary_structure.detect_beta_sheet()
        hydrogen_bonds = dssp.secondary_structure.calculate_hbonds()
        turns = dssp.secondary_structure.detect_turns(3)  

        result.append(f"Results for file: {pdb_file}\n")

        if hydrogen_bonds:
            result.append(f"Detected hydrogen bonds: {len(hydrogen_bonds)}\n")
            for bond in hydrogen_bonds:
                result.append(f"Bond between residues {bond[0]} and {bond[1]}\n")
        else:
            result.append("No hydrogen bonds detected.\n")

        if turns:
            result.append("Detected turns:\n")
            for size, turn_list in turns.items():
                result.append(f"Turns of size {size}: {len(turn_list)}\n")
                for turn in turn_list:
                    result.append(f"Turn between residues {turn[0]} and {turn[1]}\n")
        else:
            result.append("No turns detected.\n")

        result.append("Detected helices:\n")
        if helices:
            for helix_type, helix_list in helices.items():
                result.append(f"{helix_type}: {helix_list}\n")
        else:
            result.append("No helices detected.\n")

        result.append("Detected beta sheets:\n")
        if beta_sheets:
            for sheet in beta_sheets:
                result.append(f"{sheet}\n")
        else:
            result.append("No beta sheets detected.\n")

        with open(output_file, 'a') as f:
            f.writelines(result)

        logging.info(f"Finished processing: {pdb_file}")
    except Exception as e:
        logging.error(f"Error processing file {pdb_file}: {e}")

def main():
    """
    Main function to process all PDB files in a directory using multiprocessing.
    Generates an output file with the results of the secondary structure analysis.
    """
    pdb_directory = '/home/etudiant/Documents/M2BI/Programmation_avancee/projet_court_DSSP/data/'
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
    output_file = "/home/etudiant/Documents/M2BI/Programmation_avancee/projet_court_DSSP/res/output2.txt"

    file_pairs = [(os.path.join(pdb_directory, pdb_file), output_file) for pdb_file in pdb_files]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.starmap(process_file, file_pairs)

if __name__ == "__main__":
    main()
