#!/usr/bin/env python3

import os
import sys
import time
import glob
import pprint
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
         description="Generate control and osc. strength files"
         )
parser.add_argument("-d", "--input_dir", default='structures/LHCII_PDB/LHCII/LHCII_trimer/vangelis',
                    help="Directory with residues we want to renumber")
parser.add_argument("-t", "--target_dir", default='structures/LHCII_PDB/LHCII/LHCII_trimer/target',
                    help="The reference PDB directory with correctly numbered residues")
parser.add_argument("-p", "--pigment", default='CHL',
                    help="The pigment we're looking at")
args = parser.parse_args()

file_list = glob.glob(args.input_dir + '/*')
target_list = glob.glob(args.target_dir + '/*')

assert len(file_list) == len(target_list), "Number of residues not equal!"

sqdist = 0.0
atom_count = 0
rmsd_dict = {}
residues = []
targets = []

for f in file_list:
    f_name = os.path.basename(f)
    f_pigment = f_name[0:3]
    if (f_pigment != args.pigment):
        continue

    f_residue = f_name[0:-4]
    residues.append(f_residue)
    rmsd_dict[f_residue] = {}
    temp_dict = {}
    for g in target_list:
        g_name = os.path.basename(g)
        g_pigment = g_name[0:3]
        if (g_pigment != args.pigment):
            continue

        g_residue = g_name[0:-4]
        targets.append(g_residue)
        arr = np.genfromtxt(f, encoding='utf-8', dtype=None, skip_footer=1)
        target = np.genfromtxt(g, encoding='utf-8', dtype=None, skip_footer=1)
        print("Calculating RMSD for residues {} and {}: ".format(f_residue, g_residue), end="")
        sqdist = 0.0
        atom_count = 0
        for arr_line in arr:
            arr_atom_code = str("{}").format(arr_line[2])
            for target_line in target:
                target_atom_code = str("{}").format(target_line[2])
                if (arr_atom_code == target_atom_code):
                    atom_count += 1
                    r = ((float(arr_line[5]) - float(target_line[5]))**2 + (float(arr_line[6]) - float(target_line[6]))**2 + (float(arr_line[7]) - float(target_line[7]))**2)
                    sqdist += r

        sqdist /= atom_count
        rmsd = np.sqrt(sqdist)
        print(rmsd, "atom count = {}".format(atom_count))
        temp_dict[g_residue] = (rmsd, atom_count)

    rmsd_dict[f_residue] = temp_dict

targets = np.unique(targets)
f = open("{}_rmsd_data.csv".format(args.pigment), "w")
f.write(str(rmsd_dict))
f.close()

# pprint.pprint(rmsd_dict)

for key in rmsd_dict.keys():
    d = rmsd_dict[key]
    # key_max = max(d.keys(), key=(lambda k: d[k][0]))
    key_min = min(d.keys(), key=(lambda k: d[k][0]))
    print("Residue {} has minimum RMSD {} with target residue {}".format(key, d[key_min], key_min))

# for i in range(len(residues)):
#     for j in range(len(targets)):
#         rmsd_array[i, j] = rmsd_dict[residues[i]][targets[j]]

# print(rmsd_array)
# fig, ax = plt.subplots()
# im = ax.imshow(d)
# plt.imshow(d, norm=norm)
# plt.colorbar(im, label=label, norm=norm)
# plt.xlabel(xlabel)
# plt.ylabel(ylabel)
# plt.savefig(output)


# for i in range(len(residues[0])):
    # print("File residues and atom counts: ", residues[0][i], atom_counts[0][i])

# for i in range(len(residues[1])):
    # print("Target residues and atom counts: ", residues[1][i], atom_counts[1][i])
# print(squared_distances)

# np.savetxt("file_residues.txt", residues[0])
# np.savetxt("target_residues.txt", residues[1])
# np.savetxt("squared_distances.txt", squared_distances)
