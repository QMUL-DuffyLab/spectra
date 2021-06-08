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
# filelist = glob.glob(args.input_dir + '/[0-9]*')
# filelist = glob.glob(args.input_dir + '/*')
filelist = [f.path for f in os.scandir(args.input_dir) if f.is_file()]
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

ligands = {
        'CLA' : cla_dict,
        'CHL' : chl_dict,
        'LUT' : lut_dict,
        }

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
    pigments = {
        '262' : 'C/CHL601', '266' : 'E/CHL601', '270' : 'H/CHL601',
        '232' : 'C/CLA602', '238' : 'E/CLA602', '244' : 'H/CLA602',
        '233' : 'C/CLA603', '239' : 'E/CLA603', '245' : 'H/CLA603',
        '250' : 'C/CLA604', '252' : 'E/CLA604', '254' : 'H/CLA604',
        '256' : 'C/CHL605', '258' : 'E/CHL605', '260' : 'H/CHL605',
        '257' : 'C/CHL606', '259' : 'E/CHL606', '261' : 'H/CHL606',
        '263' : 'C/CHL607', '267' : 'E/CHL607', '271' : 'H/CHL607',
        '264' : 'C/CHL608', '268' : 'E/CHL608', '272' : 'H/CHL608',
        '265' : 'C/CHL609', '269' : 'E/CHL609', '273' : 'H/CHL609',
        '234' : 'C/CLA610', '240' : 'E/CLA610', '246' : 'H/CLA610',
        '235' : 'C/CLA611', '241' : 'E/CLA611', '247' : 'H/CLA611',
        '236' : 'C/CLA612', '242' : 'E/CLA612', '248' : 'H/CLA612',
        '237' : 'C/CLA613', '243' : 'E/CLA613', '249' : 'H/CLA613',
        '251' : 'C/CLA614', '253' : 'E/CLA614', '255' : 'H/CLA614',
        '274' : 'C/LUT620', '276' : 'E/LUT620', '278' : 'H/LUT620',
        '275' : 'C/LUT621', '277' : 'E/LUT621', '279' : 'H/LUT621',
        '283' : 'C/XAT622', '284' : 'E/XAT622', '285' : 'H/XAT622',
        '280' : 'C/NEX623', '281' : 'E/NEX623', '282' : 'H/NEX623',
        }
else:
    pigments = {}


'''
if we could guarantee the residues are all sequentially numbered
then i guess we wouldn't need to generate a new dictionary; but
this covers any case, since the indices are only needed to keep
track of each residue separately, they don't mean anything
'''
res_to_index = {}
output_basenames = {}
for index, key in enumerate(pigments.keys()):
    res_to_index[key] = index
    # not fully done yet - will need frame numbers later
    output_basenames[key] = args.input_dir + '/' + pigments[key]
    os.makedirs(output_basenames[key], exist_ok=True)

print(res_to_index)
print(output_basenames)

for item in filelist:
    # same deal - load in the PDB file
    f = open(item, mode='r')
    text = f.read()
    f.close()
    # final newline will create an empty list at the end, hence [:-1]
    arr = [x.split() for x in text.split("\n")[:-1]]

    # append each row of x, y, z, partial charge as we go
    temp_list = [[] for i in range(len(pigments))]
    footer_list = [[] for i in range(len(pigments))]
    output_files = [None] * len(pigments)

    dot = int((str(item)).find(".")) + 1
    suffix = str(item)[dot:]
    # pymol outputs frames in the form _xxxx.pdb
    # where xxxx is the frame number; cpptraj outputs them
    # in the form res_number.frame; this should find both
    if "pdb" in suffix:
        # note that this assumes there are no other underscores
        # after the frame number, but i think that should be safe
        frame = str(item)[int((str(item)).rfind("_")) + 1:dot - 1]
    else:
        frame = suffix

    # pymol left pads with 0's; get rid of those if they're there
    frame = str(int(frame))

    for pdb_line in arr:
        if pdb_line[0] != "ATOM":
            '''
            the continue here is just papering over a crack really.
            i was including the CONECT info in the input files for
            the code, because they're useful for the GCN later. but
            adding the correct CONECT info in the case where multiple
            residues are all in one PDB requires us to keep track of
            the atom numbers for each residue, check every line and then
            add the CONECT info for the correct atoms at the end?
            it can be done but i'm not sure it's worth doing,
            especially because i cut out the atom numbers
            '''
            continue
            # index = res_to_index[str(pdb_line[4])]
            # footer_list[index].append(pdb_line)
        else:
            index = res_to_index[str(pdb_line[4])]
            output_files[index] = output_basenames[str(pdb_line[4])]\
                    + '/frame{}.csv'.format(frame) 
            atom_code = pdb_line[2]
            try:
                tresp_dict = ligands[str(pdb_line[3])]
            except KeyError:
                continue

            # PDB file contains ligand code - use
            # that to check which tresp data to use
            if atom_code in tresp_dict.keys():
                # the 0.001 is because the values are reported as e * 10^3
                temp_list[index].append(\
                        [pdb_line[2], pdb_line[5], \
                         pdb_line[6], pdb_line[7], \
                         0.001 * float(tresp_dict[atom_code])])
            else:
                # any atom not reported in tresp file has zero charge
                temp_list[index].append(\
                        [pdb_line[2], pdb_line[5],\
                         pdb_line[6], pdb_line[7], 0.0])

    # i was using numpy for this but it really doesn't like assigning the
    # structured array and keeping the correct dtypes, for some reason, and
    # i can't be bothered to figure out why, so do this the old-fashioned way
    for i in range(len(pigments)):
        with open(output_files[i], 'w') as f:
            print("Creating file {}".format(output_files[i]))
            for line in temp_list[i]:
                f.write("{:5s} {:016.8e} {:016.8e} {:016.8e} {:016.8e}\n".format(line[0], float(line[1]), float(line[2]), float(line[3]), float(line[4])))

            for line in footer_list[i]:
                f.write("{}\n".format(" ".join(line)))
