#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Eventually - should take keyword arguments e.g. S1-Qy, S2-Qy etc. etc. and generate control file.
Also, should generate the diagonal of the matrix using oscillator strengths for each given pigment - can do this by reading strings from control file but I suppose we won't need to if we're generating the control file here.
"""

import os
import numpy as np # this is extremely wasteful to read in one float but i want to make sure everything else works first
import argparse
import subprocess
import pigment_data

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='LHCII',
        help="Relative path to input directory.")
parser.add_argument("-o", "--output_dir", default='out',
        help="Relative path to output directory.")

args = parser.parse_args()

def get_pigments(input_dir):
    pigment_dirs = []
    # programatically get pigment name/numbers
    for item in os.scandir(input_dir):
        if item.is_dir() and "_CSV" in item.name:
            pigment_dirs.append(item.name)

    return pigment_dirs


def construct_input_files(pigment_dirs, direc, snapshot_number):
    # fortran won't create the directory; do it here
    output_path = "{}/{}".format(direc, snapshot_number)
    os.makedirs(output_path, exist_ok=True)
    input_file = "in/pigments.{}".format(snapshot_number) 
    energy_file = "in/ei.txt"
    lifetimes_file = "in/lifetimes.txt"
    lambda_file = "in/lambda.txt"
    gnt_file = "in/gnt.txt"
    lineshape_file = "in/lineshapes.{}".format(snapshot_number) 
    f = open(input_file, "w")
    g = open(energy_file, "w")
    h = open(lifetimes_file, "w")
    j = open(lambda_file, "w")
    k = open(gnt_file, "w")
    l = open(lineshape_file, "w")
    for p in pigment_dirs:
        reorg = np.loadtxt("lineshape/out/{}_lambda.dat".format(p[0:3]))[0]
        gn = "lineshape/out/{}_gn.dat".format(p[0:3])
        lineshape = "lineshape/in/{}.def".format(p[0:3])
        print("{}/{}/frame{}.csv".format(input_dir, p, snapshot_number), file=f)
        print(pigment_data.pigment_data[p[0:3]]["S1"]["energy"], file=g)
        print(pigment_data.pigment_data[p[0:3]]["S1"]["lifetime"], file=h)
        print(reorg, file=j)
        print(gn, file=k)
        print(lineshape, file=l)

    f.close()
    g.close()
    h.close()
    j.close()
    return (input_file, output_path)

input_dir = os.path.join(os.getcwd(), args.input_dir)
output_dir = os.path.join(os.getcwd(), args.output_dir)
snapshot_number = 1 # replace this with for loop to iterate obv

for i in range(1000):
    pigment_dirs = get_pigments(input_dir)
    input_file, output_path = construct_input_files(pigment_dirs, output_dir, i + 1)
    print("Frame {} complete.".format(output_path))
    os.system("./coupling_calc {} {}".format(input_file, output_path))
