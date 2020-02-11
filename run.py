#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# tested on python 3.7
"""
Eventually - should take keyword arguments e.g. S1-Qy, S2-Qy etc. etc. and generate control file.
Also, should generate the diagonal of the matrix using oscillator strengths for each given pigment - can do this by reading strings from control file but I suppose we won't need to if we're generating the control file here.
"""

import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description="Generate control and osc. strength files")
parser.add_argument("-i", "--input_dir", default='LHCII',
        help="Relative path to input directory.")
parser.add_argument("-o", "--output_dir", default='hamiltonians',
        help="Relative path to output directory.")

args = parser.parse_args()

def get_pigments(input_dir):
    pigment_dirs = []
    # programatically get pigment name/numbers
    for item in os.scandir(input_dir):
        if item.is_dir() and "_CSV" in item.name():
            pigment_dirs.append(item.name())

    return pigment_dirs


def construct_input_file(pigment_dirs, output_dir, snapshot_number):
    # fortran won't create the directory; do it here
    os.makedirs("{}/{}".format(output_dir, snapshot_number))
    input_file = "pigments.{}".format(snapshot_number) 
    f = open(input_file, "w")
    for p in pigment_dirs:
        print("{}/frame.{}".format(p, snapshot_number), file=f)

    f.close()
    # f = open("parameters.{}".format(snapshot_number), "w")
    # print('''
    # pigment_file = {0}
    # output_dir = {1}
    # '''.format(), file=f)
    return input_file

input_dir = os.path.join(os.cwd(), args.input_dir)
output_dir = os.path.join(os.cwd(), args.output_dir)
snapshot_number = 1 # replace this with for loop to iterate obv

pigment_dirs = get_pigments(input_dir)
input_file = construct_input_file(pigment_dirs, output_dir, snapshot_number)

subprocess.popen(exe, input_file, output_dir)
