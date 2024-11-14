# Import necessary modules
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
import xml.etree.ElementTree as ET
import numpy as np
import glob

# Script to convert VASP outputs(MTP fake vasprun.xml) to a single DFSET file for ALAMODE
# author: Jiongzhi Zheng
# email: jiongzhi.zheng@dartmouth.edu
# institution: Dartmouth College
# date: 2024-11-13
# version: 1.0

# usage: python3 mtp_to_alamode_dfset.py
# This script reads the POSCAR file and a series of vasprun_*.xml files in the directory and writes the data to a single DFSET file.
# The DFSET file contains the Cartesian coordinates of the atoms, the forces on the atoms, 
# and the displacements of the atoms from the equilibrium positions.

# The displacements are adjusted to be within the -0.5 to 0.5 range.
# The units are converted from Angstrom to Bohr for positions and from eV/Angstrom to Ryd/bohr for forces.


# Define the path to the directory containing the files
path = "/Users/jiongzhizheng/Desktop/harm_vasp_outputs/"

# Read the POSCAR file and initialize Phonopy object
unitcell = read_vasp(path + "POSCAR")
# change the supercell size to your supercell size
phonon = Phonopy(unitcell, [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

# Get the supercell and its lattice vectors
supercell = phonon.get_supercell()
lattice = supercell.get_cell()
lattice_inv = np.linalg.inv(lattice)

positions_t = supercell.get_scaled_positions()
positions_t = np.dot(positions_t, lattice)

# Function to extract data from a specific varray
def extract_varray_data(root, varray_name):
    data = []
    for v in root.findall(f".//varray[@name='{varray_name}']/v"):
        values = list(map(float, v.text.split()))
        data.append(values)
    return data

# Prepare to write the combined data to a single DFSET file
with open(path + "DFSET", "w") as f:

    # Loop over each vasprun_*.xml file in the directory
    for vasprun_file in glob.glob(path + "vasprun_*.xml"):
        # Parse the XML file
        tree = ET.parse(vasprun_file)
        root = tree.getroot()

        # Extract positions and forces
        positions = np.array(extract_varray_data(root, "positions"))
        forces = np.array(extract_varray_data(root, "forces"))

        # Convert positions to fractional coordinates
        positions_frac = np.dot(positions, lattice_inv)
        disps = positions_frac - phonon.supercell.get_scaled_positions()

        # Adjust displacements to be within -0.5 to 0.5 range
        disps = np.where(disps > 0.5, disps - 1.0, disps)
        disps = np.where(disps < -0.5, disps + 1.0, disps)

        # Convert fractional displacements to Cartesian coordinates
        positions_cart = np.dot(disps, lattice)

        # Convert units: positions from Angstrom to Bohr and forces from eV/Angstrom to Ryd/bohr
        positions_cart *= 1.88973  # Angstrom to Bohr
        forces *= 0.038935     # eV/Angstrom to Ryd/bohr

        # Write the data to DFSET with a header line for each file
        f.write(f"# {vasprun_file}\n")
        for i in range(len(positions_cart)):
            f.write(f"{positions_cart[i][0]:.16f} {positions_cart[i][1]:.16f} {positions_cart[i][2]:.16f}     "
                    f"{forces[i][0]:.16f} {forces[i][1]:.16f} {forces[i][2]:.16f}\n")
