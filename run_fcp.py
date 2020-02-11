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
ei_file = "ei.txt"
lambda_file = "lambda.txt"
gnt_file = "gnt.txt"
lifetimes_file = "lifetimes.txt"
lineshape_dir = "/Users/cgray/code/lineshape"

'''
comments on following excitation energies:
    chl a ones are taken from valkunas 2013 jcpl
    chl c taken from peak locations listed on wikipedia lol
    fuco and diadino :shrug:, although enriquez 2010 gives
    diadino S1 at about 14170
'''
ei = {
        'Chla_Qy': 14900.0, 'Chla_Qx': 16300.0,
        'Chlc_Qy': 15970.0, 'Chlc_Qx': 17330.0,
        'Fuco_S2': 20000.0, 'Diad_S2': 20000.0,
        'Fuco_S1': 14900.0, 'Diad_S1': 14900.0
                 }

# reorganisation energies - these are specific to the parameters used
# in the spectral density functions
# reorgs = {
#         'Chla_Qy': 14900.0, 'Chla_Qx': 16300.0,
#         'Chlc_Qy': 15970.0, 'Chlc_Qx': 17330.0,
#         'Fuco_S2': 20000.0, 'Diad_S2': 20000.0,
#         'Fuco_S1': 14900.0, 'Diad_S1': 14900.0
#                 }

# excited state lifetimes to put into the exponentials later
# these are in nanoseconds!
lifetimes = {
        'Chla_Qy': 4, 'Chla_Qx': 4,
        'Chlc_Qy': 4, 'Chlc_Qx': 4,
        'Fuco_S2': 0.01, 'Diad_S2': 0.01,
        'Fuco_S1': 0.01, 'Diad_S1': 0.01
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
chl_tresp_list = ["{}_anion".format(p[0:-4]) if "Chlc" in p else p[0:-3] for p in chl_list]
car_tresp_list = [p[0:-4] for p in car_list]

chl_lig_list = ["CLA" if "Chla" in p else "CLC" for p in chl_list]
car_lig_list = ["A86" if "Fuco" in p else "DD6" for p in car_list]
chl_lambda_list = ["{}/out/{}_lambda.dat".format(lineshape_dir, p) for p in chl_lig_list]
car_lambda_list = ["{}/out/{}_lambda.dat".format(lineshape_dir, p) for p in car_lig_list]
chl_gnt_list = ["{}/out/{}_gt.dat".format(lineshape_dir, p) for p in chl_lig_list]
car_gnt_list = ["{}/out/{}_gt.dat".format(lineshape_dir, p) for p in car_lig_list]

# the lambdas are just one number - read them in here, makes it easier in the fortran
for i, l in enumerate(chl_lambda_list):
    with open(l) as f:
        # ugly but ¯\_(ツ)_/¯
        chl_lambda_list[i] = float((f.readline()).split()[0])

for i, l in enumerate(car_lambda_list):
    with open(l) as f:
        car_lambda_list[i] = float((f.readline()).split()[0])

pdb_list = []
tresp_list = []
lambda_list = []
gnt_list = []
state = []
for s in args.states:
    if 'Q' in s:
        pdb_list.append(chl_pdb_list)
        lambda_list.append(chl_lambda_list)
        gnt_list.append(chl_gnt_list)
        tresp_list.append(["tresp/{0}_{1}_TrESP.txt".format(tr, s) for tr in chl_tresp_list])
        state.append([s for p in chl_pdb_list])
    elif 'S' in s:
        pdb_list.append(car_pdb_list)
        lambda_list.append(car_lambda_list)
        gnt_list.append(car_gnt_list)
        tresp_list.append(["tresp/{0}_{1}_TrESP.txt".format(tr, s) for tr in car_tresp_list])
        state.append([s for p in car_pdb_list])


# need to flatten the lists for the loop below
# there must be a way of doing these all at once
flatten = lambda l : [y for x in l for y in x]
pdb_list = flatten(pdb_list)
lambda_list = flatten(lambda_list)
gnt_list = flatten(gnt_list)
tresp_list = flatten(tresp_list)
state = flatten(state)

f = open(control_file, "w")
g = open(ei_file, "w")
h = open(lambda_file, "w")
j = open(gnt_file, "w")
k = open(lifetimes_file, "w")
for i, var in enumerate(pdb_list):
    print("{} {}".format(pdb_list[i], tresp_list[i]), file=f)
    print(ei["{}_{}".format(var[7:11], state[i])], file=g)
    print(lifetimes["{}_{}".format(var[7:11], state[i])], file=k)
    print("{}".format(lambda_list[i]), file=h)
    print("{}".format(gnt_list[i]), file=j)

f.close()
g.close()
h.close()
j.close()
k.close()
