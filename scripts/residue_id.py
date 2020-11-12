#!/usr/bin/env python3

import os
import sys
import time
import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser(
         description="Generate control and osc. strength files"
         )
parser.add_argument("-f", "--file", default='structures/LHCII',
                    help="File with residues we want to renumber")
parser.add_argument("-t", "--target", default='LHCII',
                    help="The reference PDB file with correctly numbered residues")
args = parser.parse_args()

arr = np.genfromtxt(args.file, encoding='utf-8', dtype=None, skip_footer=1)
target = np.genfromtxt(args.target, encoding='utf-8', dtype=None, skip_footer=2)
residues = [np.unique([str("{}{}").format(line[3], str(line[4])) for line in file]) for file in [arr, target]]
assert len(residues[0]) == len(residues[1]), "Number of residues not equal!"
print("Residues to renumber: ", residues[0])
print("Target residues: ", residues[1])

squared_distances = np.zeros((len(residues[0]), len(residues[1])))
atom_counts = np.zeros((2, len(residues[0])))

t0 = time.time_ns()
for arr_line in arr:
    arr_residue = str("{}{}").format(arr_line[3], str(arr_line[4]))
    arr_pigment = str("{}").format(arr_line[3])
    arr_atom_code = str("{}").format(arr_line[2])
    arr_index = np.where(residues[0] == arr_residue)
    atom_counts[0, arr_index] += 1
    for target_line in target:
        target_pigment = str("{}").format(target_line[3])
        target_residue = str("{}{}").format(target_line[3], str(target_line[4]))
        target_atom_code = str("{}").format(arr_line[2])
        target_index = np.where(residues[1] == target_residue)
        atom_counts[1, target_index] += 1
        if (arr_pigment == target_pigment):
            if (arr_atom_code == target_atom_code):
                sqdist = ((float(arr_line[5]) - float(target_line[5]))**2 + (float(arr_line[6]) - float(target_line[6]))**2 + (float(arr_line[7]) - float(target_line[7]))**2)
                squared_distances[arr_index, target_index] += sqdist


t1 = time.time_ns()
print("Total time taken (s): {:6.3f}".format(float(t1 - t0) / 1000000000))

for i in range(len(residues[0])):
    print("File residues and atom counts: ", residues[0][i], atom_counts[0][i])

for i in range(len(residues[1])):
    print("Target residues and atom counts: ", residues[1][i], atom_counts[1][i])
print(squared_distances)

np.savetxt("file_residues.txt", residues[0])
np.savetxt("target_residues.txt", residues[1])
np.savetxt("squared_distances.txt", squared_distances)
