#!/usr/bin/env python3

import os
import glob
import numpy as np

filelist = glob.glob('LHCII' + '/**/*.*', recursive=True)
# want to keep the atom strings so don't loadtxt!
cla_tresp = np.genfromtxt("in/cla_s01_no_h.dat", dtype=None)
chl_tresp = np.genfromtxt("in/chl_s01_no_h.dat", dtype=None)
cla_tresp = [list(line) for line in zip(*cla_tresp)]
chl_tresp = [list(line) for line in zip(*chl_tresp)]

'''
   there's probably an easier way of doing this but
   basically we turn the tresp data into a dictionary:
   then we can check each line of the PDB file and get the
   partial charge for the corresponding atom. It's easier 
   to do this way because the tresp file doesn't list atoms
   for which the charge is zero, and also the atoms aren't
   listed in the same order.
'''
cla_dict = dict(zip(cla_tresp[1], cla_tresp[2]))
chl_dict = dict(zip(chl_tresp[1], chl_tresp[2]))

for item in filelist:
    # same deal - load in the PDB file
    arr = np.genfromtxt(item, dtype=None, skip_footer=2)
    # we'll append each row of x, y, z, partial charge as we go
    temp_list = []
    print(item)
    # PDB file contains ligand code - use 
    # that to check which tresp data to use
    ligand_code = str(arr[0][3], 'utf-8')
    if ligand_code == "CLA":
        tresp_dict = cla_dict
    elif ligand_code == "CHL":
        tresp_dict = chl_dict
    else:
        raise ValueError("Invalid ligand code: {}".format(ligand_code))

    for pdb_line in arr:
        atom_code = pdb_line[2]
        # could create a lookup table for this
        # but i doubt it's worth the effort
        if atom_code in tresp_dict.keys():
            # the 0.001 is because the values are reported as e * 10^3
            temp_list.append([pdb_line[5], pdb_line[6], 
                pdb_line[7], 0.001 * float(tresp_dict[atom_code])])
        else:
            # any atom not reported in tresp file has zero charge
            temp_list.append([pdb_line[5], pdb_line[6], pdb_line[7], 0.0])

    output = np.array(temp_list)
    test_output = "{}_test".format(str(item))
    # this fmt call is to make all the data regular, because
    # otherwise it's an absolute nightmare to read in fortran
    np.savetxt(test_output, output, fmt='%016.8e')
