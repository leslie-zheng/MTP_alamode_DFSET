# MTP_alamode_DFSET

Script to convert VASP outputs(MTP fake vasprun.xml) to a single DFSET file for ALAMODE
author: Jiongzhi Zheng
email: jiongzhi.zheng@dartmouth.edu
institution: Dartmouth College
date: 2024-11-13
version: 1.0

usage: python3 mtp_to_alamode_dfset.py
This script reads the POSCAR file and a series of vasprun_*.xml files in the directory and writes the data to a single DFSET file.
The DFSET file contains the Cartesian coordinates of the atoms, the forces on the atoms, 
and the displacements of the atoms from the equilibrium positions.

The displacements are adjusted to be within the -0.5 to 0.5 range.
The units are converted from Angstrom to Bohr for positions and from eV/Angstrom to Ryd/bohr for forces.
