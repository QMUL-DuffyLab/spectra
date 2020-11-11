#!/usr/bin/env python3

import os
import sys
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

arr = np.genfromtxt(args.file, dtype=None, skip_header=5, skip_footer=2)
target = np.genfromtxt(args.target, dtype=None, skip_header=4, skip_footer=2)
initial_residues = []
target_residues = []
for line in arr:
    initial_residues.append(str("{}{}".format(line[3], str(line[4]))))

for line in target:
    target_residues.append(str("{}{}".format(line[3], str(line[4]))))

initial_residues = np.unique(np.array(initial_residues))
target_residues = np.unique(np.array(target_residues))
# initial_residues = np.unique(np.array(["".join(str(line[3]) + str(line[4]) for line in arr)]))
print(initial_residues)
print(target_residues)

# for item in filelist:
#     # same deal - load in the PDB file
#     arr = np.genfromtxt(item, encoding='utf-8', dtype=None, skip_footer=2)
#     # we'll append each row of x, y, z, partial charge as we go
#     temp_list = []
#     # PDB file contains ligand code - use
#     # that to check which tresp data to use
#     output_dir = args.input_dir + '/' + pigments[str(arr[0][4])]
#     os.makedirs(output_dir, exist_ok=True)
#     frame = str(item)[int((str(item)).find(".")) + 1:]
#     output_file = output_dir + '/frame{}.csv'.format(frame)
#     print("Creating file {}".format(output_file))

#     ligand_code = str(arr[0][3])
#     if ligand_code == "CLA":
#         tresp_dict = cla_dict
#     elif ligand_code == "CHL":
#         tresp_dict = chl_dict
#     elif ligand_code == "LUT":
#         tresp_dict = lut_dict
#     else:
#         continue
#         # raise ValueError("Invalid ligand code: {}".format(ligand_code))

#     for pdb_line in arr:
#         atom_code = pdb_line[2]
#         # could create a lookup table for this
#         # but i doubt it's worth the effort
#         if atom_code in tresp_dict.keys():
#             # the 0.001 is because the values are reported as e * 10^3
#             temp_list.append([pdb_line[5], pdb_line[6], pdb_line[7],
#                               0.001 * float(tresp_dict[atom_code])])
#         else:
#             # any atom not reported in tresp file has zero charge
#             temp_list.append([pdb_line[6], pdb_line[7], pdb_line[8], 0.0])

#     output = np.array(temp_list)
#     # this fmt call is to make all the data regular, because
#     # otherwise it's an absolute nightmare to read in fortran
#     np.savetxt(output_file, output, fmt='%016.8e')
