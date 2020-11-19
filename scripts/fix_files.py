#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import numpy as np

parser = argparse.ArgumentParser(
         description="Generate control and osc. strength files"
         )
parser.add_argument("-i", "--input_dir", default='structures/LHCII',
                    help="Relative path to input directory.")
parser.add_argument("-p", "--protein", default='LHCII',
                    help="the protein we're looking at")
args = parser.parse_args()

# assumes the PDB files we're fixing have filenames in the form
# R0.frame where R0 is a sequence of digits which represent the
# residue number from the MD trajectory and .frame is the frame number.
# this is how I output the files from cpptraj so should be a safe assumption
filelist = glob.glob(args.input_dir + '/[0-9]*')
print("{} files found.\n".format(len(filelist)))

# want to keep the atom strings so don't loadtxt!
cla_tresp = np.genfromtxt("in/cla_tresp.dat", comments='#',
                          encoding='utf-8', dtype=None)
chl_tresp = np.genfromtxt("in/chl_tresp.dat", comments='#',
                          encoding='utf-8', dtype=None)
lut_tresp = np.genfromtxt("in/lut_tresp.dat", comments='#',
                          encoding='utf-8', dtype=None)
cla_tresp = [list(line) for line in zip(*cla_tresp)]
chl_tresp = [list(line) for line in zip(*chl_tresp)]
lut_tresp = [list(line) for line in zip(*lut_tresp)]

# ugly but it's just translating the residue numbers from the MD
# into the equivalent ones from the Liu structure
if (args.protein) == 'LHCII':
    pigments = {
        '1433': 'CHL601',
        '1434': 'CHL605',
        '1435': 'CHL606',
        '1436': 'CHL607',
        '1437': 'CHL608',
        '1438': 'CHL609',
        '1439': 'CLA602',
        '1440': 'CLA603',
        '1441': 'CLA604',
        '1442': 'CLA610',
        '1443': 'CLA611',
        '1444': 'CLA612',
        '1445': 'CLA613',
        '1446': 'CLA614',
        '1447': 'LUT620',
        '1448': 'LUT621',
        '1449': 'NEX623',
        '1450': 'XAT622',
        }
elif (args.protein) == 'NLLZ':
    pigments = {
        '1433': 'CHL601',
        '1434': 'CHL605',
        '1435': 'CHL606',
        '1436': 'CHL607',
        '1437': 'CHL608',
        '1438': 'CHL609',
        '1439': 'CLA602',
        '1440': 'CLA603',
        '1441': 'CLA604',
        '1442': 'CLA610',
        '1443': 'CLA611',
        '1444': 'CLA612',
        '1445': 'CLA613',
        '1446': 'CLA614',
        '1447': 'LUT620',
        '1448': 'LUT621',
        '1449': 'NEX623',
        '1450': 'ZEA622',
        }
elif (args.protein) == 'VANGELIS':
    pigments = {}
else:
    pigments = {}

'''
   there's probably an easier way of doing this but
   basically we turn the tresp data into a dictionary:
   then we can check each line of the PDB file and get the
   partial charge for the corresponding atom. It's easier
   to do this way because the tresp file doesn't list atoms
   for which the charge is zero, and also the atoms aren't
   listed in the same order.
'''
cla_dict = dict(zip(cla_tresp[0], cla_tresp[1]))
chl_dict = dict(zip(chl_tresp[0], chl_tresp[1]))
lut_dict = dict(zip(lut_tresp[0], lut_tresp[1]))

for item in filelist:
    # same deal - load in the PDB file
    arr = np.genfromtxt(item, encoding='utf-8', dtype=None, skip_footer=2)
    # we'll append each row of x, y, z, partial charge as we go
    temp_list = []
    # PDB file contains ligand code - use
    # that to check which tresp data to use
    if (args.protein == 'VANGELIS'):
        res = "{}{}".format(arr[0][3], str(arr[0][4]))
        if (pigments[res][3] == '1'):
            output_dir = args.input_dir + '/B/' + pigments[res]
        elif (pigments[res][3] == '5'):
            output_dir = args.input_dir + '/F/' + pigments[res]
        elif (pigments[res][3] == '6'):
            output_dir = args.input_dir + '/G/' + pigments[res]
    else:
        output_dir = args.input_dir + '/' + pigments[str(arr[0][4])]

    os.makedirs(output_dir, exist_ok=True)
    frame = str(item)[int((str(item)).find(".")) + 1:]
    output_file = output_dir + '/frame{}.csv'.format(frame)
    print("Creating file {}".format(output_file))

    ligand_code = str(arr[0][3])
    if ligand_code == "CLA":
        tresp_dict = cla_dict
    elif ligand_code == "CHL":
        tresp_dict = chl_dict
    elif ligand_code == "LUT":
        tresp_dict = lut_dict
    else:
        continue
        # raise ValueError("Invalid ligand code: {}".format(ligand_code))

    for pdb_line in arr:
        atom_code = pdb_line[2]
        # could create a lookup table for this
        # but i doubt it's worth the effort
        if atom_code in tresp_dict.keys():
            # the 0.001 is because the values are reported as e * 10^3
            temp_list.append([pdb_line[5], pdb_line[6], pdb_line[7],
                              0.001 * float(tresp_dict[atom_code])])
        else:
            # any atom not reported in tresp file has zero charge
            temp_list.append([pdb_line[6], pdb_line[7], pdb_line[8], 0.0])

    output = np.array(temp_list)
    # this fmt call is to make all the data regular, because
    # otherwise it's an absolute nightmare to read in fortran
    np.savetxt(output_file, output, fmt='%016.8e')
