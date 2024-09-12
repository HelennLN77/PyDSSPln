# PyDSSPln - DSSP implementation project
This project is an implementation of the DSSP ( Dictionary of Secondary Structure of Proteins) algorithm. Its role is to parse an input Protein Data Bank ( PDB) file to determine the secondary structure of proteins (helices, sheets) The projects aims to imlemente the DSSP algorithm and output the results in .mkdssp format. It can be used to study the protein structures and compare the accuracy, sensitivity, specificity with the offical DSSP toom from the BioPython library.

## Features
Hydrogen Bond Calculation: Calculates and detects hydrogen bonds between residues.
Helix Detection: Identifies alpha helices of different sizes (3, 4, and 5 residues).
Beta Sheet Detection: Detects beta sheets and their type (parallel or antiparallel).
Turn Detection: Identifies turns of size 3 residues.

## Requirements
Python 3.x
NumPy
Biopython
scikit-learn (for metrics)

# Setup
To install all the packages, you need to perform the folowwing steps :

### Clone the repository 
```bash
git clone https://github.com/HelennLN77/PyDSSPln.git
cd PyDSSPln
```
### Install conda 
### Create a conda environment 
```bash
conda env create -f pydssp_ln.yml
```
### Activate  the conda environment 
```bash
conda activate pydssp_ln
```
## Run the script 
```bash
python PyDSSPln.py
```
## Output 
The results will be saved to the specified output file. The output will include:

Detected hydrogen bonds
Detected alpha helices
Detected beta sheets and their types
Detected turns
