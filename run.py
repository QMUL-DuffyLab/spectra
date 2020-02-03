#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Eventually - should take keyword arguments e.g. S1-Qy, S2-Qy etc. etc. and generate control file.
Also, should generate the diagonal of the matrix using oscillator strengths for each given pigment - can do this by reading strings from control file but I suppose we won't need to if we're generating the control file here.
"""

import os
import sys
import itertools
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-s", "--states", nargs="*", default=['Qy', 'S1', 'Qx', 'S2'],
        help="List of states (default is Qx, Qy, S1, S2) to add to the control file.")

args = parser.parse_args()

control_file = "J_control.txt"
osc_file = "osc.txt"
osc_strengths = {
        'Chla_Qy': 14900.0, 'Chla_Qx': 5000.0,
        'Chlc_Qy': 10000.0, 'Chlc_Qx': 2500.0, 
        'Fuco_S2': 20000.0, 'Diad_S2': 20000.0,
        'Fuco_S1': 14950.0, 'Diad_S1': 14950.0
                 }

'''
these lists are ugly; they come from naming conventions for pdb/tresp files.
it'd be nice if I didn't have to hardcode them but I don't think there's
any way of knowing the numbers for the PDB files beforehand?
'''
chl_list = ['Chla401', 'Chla402', 'Chla404', 'Chla405', 'Chla406', 'Chla407', 'Chla409', 'Chlc2_403', 'Chlc1_408']
car_list = ['Fuco_301', 'Fuco_302', 'Fuco_303', 'Fuco_304', 'Fuco_305', 'Fuco_306', 'Fuco_307', 'Diadino_308']
chl_pdb_list = ["coords/{}.pdb".format(p) for p in chl_list]
car_pdb_list = ["coords/{}.pdb".format(p) for p in car_list]
# chl_tresp_list = ['Chla', 'Chla', 'Chla', 'Chla', 'Chla', 'Chla', 'Chla', 'Chlc2_anion', 'Chlc1_anion']
chl_tresp_list = ["{}_anion".format(p[0:-4]) if "Chlc" in p else p[0:-4] for p in chl_list]
car_tresp_list = [p[0:-4] for p in car_list]

pdb_list = []
tresp_list = []
state = []
for s in args.states:
    if 'Q' in s:
        pdb_list.append(chl_pdb_list)
        tresp_list.append(["tresp/{0}_{1}_TrESP.txt".format(tr, s) for tr in chl_tresp_list])
        state.append([s for p in chl_pdb_list])
    elif 'S' in s:
        pdb_list.append(car_pdb_list)
        tresp_list.append(["tresp/{0}_{1}_TrESP.txt".format(tr, s) for tr in car_tresp_list])
        state.append([s for p in car_pdb_list])


pdb_list = [y for x in pdb_list for y in x]
tresp_list = [y for x in tresp_list for y in x]
state = [y for x in state for y in x]
print(pdb_list)
print(tresp_list)
f = open(control_file, "w")
g = open(osc_file, "w")
for i, var in enumerate(pdb_list):
    print("{} {}".format(pdb_list[i], tresp_list[i]), file=f)
    print(osc_strengths["{}_{}".format(var[7:11], state[i])], file=g)

f.close()
g.close()
