#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Generate parameter and job files for each combination of parameters
passed to this script, and run them locally or on HPC clusters.
"""

import os
import sys
import itertools
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Run script for FCP lineshape code")
# non-parameter args
parser.add_argument("-r", "--run", default=1, type=int, help="Defaults to 1; set to 0 to prevent running of program")
parser.add_argument("-l", "--ligand", help="Name of ligand")

args = parser.parse_args()

# separate into parameter and non-parameter dicts:
# makes it easier to iterate over

# sometimes i add branch info to the executable name
prog_name = os.popen("grep 'EXEC = ' makefile | awk '{print $3}'").read().rstrip('\n')
label = 'out/{}'.format(args.ligand)
p['Aw_file'] = os.path.join(prefix, label + "Aw.dat")
p['Fw_file'] = os.path.join(prefix, label + "Fw.dat")
p['At_file'] = os.path.join(prefix, label + "At.dat")
p['Ft_file'] = os.path.join(prefix, label + "Ft.dat")
p['lambda_file'] = os.path.join(prefix, label + "lambda.dat")
f = open(temp_file, "w")
for key in p.keys():
    print("{} = {}".format(key, p[key]), file=f)

f.close()
run_cmd = "{0} ./{1} {2} & ".format(
         prog_name, ligand_file, npd['protocol_file'])
if npd['run'] == 1:
    os.system(run_cmd)
else:
    print(run_cmd)
