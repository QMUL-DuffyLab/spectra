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

# (optional) parameter arguments
# default=[] ensures the dict comprehension below works
# even if none of these arguments are included
parser.add_argument("--nx", nargs='*', default=[],
                    help="List of nx values to use")

args = parser.parse_args()

# separate into parameter and non-parameter dicts:
# this makes it easier to iterate over
non_params = ['run', 'clusters', 'no_slots', 'analyse',
        'input_file', 'protocol_file', 'output_dir']
npd = dict((k, v) for k, v in vars(args).items() if k in non_params)
pd  = dict((k, v) for k, v in vars(args).items() if k not in non_params)

# read in the input file
p = {}
with open(npd['input_file']) as f:
    for line in f:
        if line[0] == '#':
            continue
        (key, val) = line.split(' = ')
        p[key] = val.rstrip('\n')

# relative path from cwd
prefix = os.path.join(os.getcwd(), npd['output_dir'])

# concatenate with the params in input_file and take unique values
# the p[key] has to be enclosed in []; have to add two lists
d = {key:list(set((pd[key] + [p[key]]))) for key in pd.keys()}
# d.keys() isn't indexable but i want to index it below
d_keys = list(d.keys())
# every permutation, one from each list of optional parameters
perm = list(itertools.product(*d.values()))

# sometimes i add branch info to the executable name
prog_name = os.popen("grep 'EXEC = ' makefile | awk '{print $3}'").read().rstrip('\n')

for i in range(len(perm)):
    # label output files by parameters
    label = ''
    for j in range(len(perm[i])):
        p[d_keys[j]] = perm[i][j]
        label = label + "{}_{}_".format(d_keys[j], perm[i][j])
    
    p['mag'] = os.path.join(prefix, label + "ac_mag.out")
    p['stats'] = os.path.join(prefix, label + "stats.out")

    temp_file = os.path.join(os.getcwd(),
                os.path.dirname(npd['input_file']), label + "params")
    f = open(temp_file, "w")
    for key in p.keys():
        print("{} = {}".format(key, p[key]), file=f)

    f.close()
    if npd['clusters'] == 1:
        job_file = os.path.join(os.getcwd(),
                   os.path.dirname(npd['input_file']), label[:-1] + ".job")
        jf = open(job_file, "w")
        job_script = print("""#!/bin/bash -f
# =================================
#$ -M callum.gray.10@ucl.ac.uk
#$ -m bes
#$ -V
#$ -j y
#$ -cwd
#$ -N '{0}'
#$ -S /bin/bash
#$ -l vf={1}
#$ -pe ompi {2}
#
IPWD=`pwd`
echo "Got ${{NSLOTS}} slots in ${{IPWD}}."
mpirun --mca btl ^openib --mca mtl ^psm --n ${{NSLOTS}} ./{3} {4} {5}
exit 0
        """.format(label[:-1],ram,npd['no_slots'],prog_name,
                   temp_file,npd['protocol_file']), file=jf)
        jf.close()
        run_cmd = "qsub -q {0} {1}".format(queue, job_file)
        if npd['run'] == 1:
            os.system(run_cmd)
        else:
            print(run_cmd)

    else:
        run_cmd = "mpirun -np {0} ./{1} {2} {3} | tee {2}.out & ".format(npd['no_slots'],
                 prog_name, temp_file, npd['protocol_file'])
        if npd['run'] == 1:
            os.system(run_cmd)
        else:
            print(run_cmd)
