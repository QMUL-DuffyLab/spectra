#!/usr/bin/env python3

import os
import glob
import pprint
import pickle
import argparse
import numpy as np

parser = argparse.ArgumentParser(
         description="Generate control and osc. strength files"
         )
parser.add_argument("-d", "--input_dir", default='structures/LHCII_PDB/LHCII/LHCII_trimer/vangelis',
                    help="Directory with residues we want to renumber")
parser.add_argument("-t", "--target_dir", default='structures/LHCII_PDB/LHCII/LHCII_trimer/target',
                    help="The reference PDB directory with correctly numbered residues")
args = parser.parse_args()

file_list = glob.glob(args.input_dir + '/*')
target_list = glob.glob(args.target_dir + '/*')

assert len(file_list) == len(target_list), "Number of residues not equal!"

sqdist = 0.0
atom_count = 0
rmsd_dict = {}
residues = []
targets = []
pigments = ['CLA', 'CHL', 'LUT', 'NEX', 'XAT']
# atoms we want to "fix" to identify pairs of residues
# leave blank to sum over every common atom
atoms = ['MG', 'NA', 'NB', 'NC', 'ND']

for pigment in pigments:
    for f in target_list:
        f_name = os.path.basename(f)
        if (f_name[0:3] != pigment):
            continue

        f_residue = f_name[0:-4]
        residues.append(f_residue)
        # rmsd_dict[f_residue] = {}
        temp_dict = {}
        for g in file_list:
            g_name = os.path.basename(g)
            if (g_name[0:3] != pigment):
                continue

            g_residue = g_name[0:-4]
            targets.append(g_residue)
            arr = np.genfromtxt(f, encoding='utf-8', dtype=None)
            target = np.genfromtxt(g, encoding='utf-8', dtype=None)
            # print("Calculating RMSD for residues {} and {}: ".format(f_residue, g_residue), end="")
            sqdist = 0.0
            atom_count = 0
            for arr_line in arr:
                arr_atom_code = (str("{}").format(arr_line[2])).replace(" ", "")
                if (arr_atom_code in atoms or len(atoms) == 0):
                    for target_line in target:
                        target_atom_code = (str("{}").format(target_line[2])).replace(" ", "")
                        if (arr_atom_code == target_atom_code):
                            atom_count += 1
                            try:
                                r = ((float(arr_line[5]) - float(target_line[5]))**2 + (float(arr_line[6]) - float(target_line[6]))**2 + (float(arr_line[7]) - float(target_line[7]))**2)
                                sqdist += r
                            except ValueError:
                                print("Value error in RMSD calc - relevant info:")
                                print(arr_line)
                                print(target_line)

            if (atom_count > 0):
                sqdist /= atom_count
                rmsd = np.sqrt(sqdist)
                # print("{:12.8f},".format(rmsd), "atom count = {:3d}".format(atom_count))
                temp_dict[g_residue] = (rmsd, atom_count)
                rmsd_dict[f_residue] = temp_dict


targets = np.unique(targets)
assert len(residues) == len(targets), "Number of residues doesn't equal number of targets???"

rmsd_array = np.zeros((len(residues), len(targets)))
for i, r in enumerate(residues):
    if (r in rmsd_dict.keys()):
        d = rmsd_dict[r]
    else:
        d = {}
    # list(targets) because targets is an np.array() and we want to use it as one later
    for j, t in enumerate(list(targets)):
        if t in d.keys():
            rmsd_array[i, j] = d[t][0]
        else:
            rmsd_array[i, j] = 0.0

with open("rmsd.pkl", "wb") as handle:
    pickle.dump(rmsd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

np.savetxt("rmsd_array.csv", rmsd_array)

counts = np.zeros(len(targets), dtype=int)
pairs = {}
for key in sorted(rmsd_dict):
    d = rmsd_dict[key]
    key_min = min(d, key=d.get)
    counts[np.where(targets == key_min)] += 1
    # there's probably a way to format the bits of the tuple along with the strings but I tried a few different things and none worked. so format the tuple separately
    formatted_tuple_string = "{0:12.8f}, {1:3d}".format(*d[key_min])
    pairs[key] = key_min
    print("Residue {} has minimum RMSD {} with target residue {}".format(key, formatted_tuple_string, key_min))

with open("pairs.pkl", "wb") as handle:
    pickle.dump(pairs, handle, protocol=pickle.HIGHEST_PROTOCOL)

pprint.pprint(pairs)
# this gives a list of tuples - the first in each tuple is the target residue (the correctly numbered one), and the second is the incorrectly numbered residue which has a minimum RMSD with it.
duplicates = [(value, key) for val in targets[np.where(counts > 1)] for key, value in pairs.items() if value == val]
zeroes = targets[np.where(counts == 0)]
print("Targets not chosen: ", zeroes)

# now we turn that into a dict - we can look at each of the correctly numbered residues and see which of the incorrectly numbered ones have been assigned to it
inverse_dict = {}
for (key, value) in duplicates:
    if not key in inverse_dict:
        inverse_dict[key] = [value]
    else:
        inverse_dict[key].append(value)

pprint.pprint(inverse_dict)
# for (key, value) in inverse_dict.items():
#     for res in value:
        # print("Residue:", res)
        # pprint.pprint(sorted(rmsd_dict[res].items(), key=lambda item: item[1]))

# now the ones that weren't chosen
for res in zeroes:
    print("Target: ", res)
    dists = []
    for key in rmsd_dict.keys():
        if (res in rmsd_dict[key].keys()):
            dists.append((key, rmsd_dict[key][res][0]))

    pprint.pprint(sorted(dists, key=lambda item: item[1]))
